"""One-off utility to repair DATE-OBS in previously rendered EUV map FITS files."""

import argparse
import base64
from dataclasses import dataclass
import gzip
import json
import math
import os
from pathlib import Path
import re
import sys
import warnings

import astropy.units as u
from astropy.io import fits

from suncet_instrument_simulator.make_radiance_maps import model_frame_observation_time


RADIANCE_MAP_PATTERN = re.compile(r'radiance_maps_(\d+)\.fits')
BACKUP_SUFFIX = '.timestamp-headers.json.gz'


@dataclass(frozen=True)
class RepairResult:
    path: Path
    expected_date_obs: str
    previous_dates: tuple
    needs_update: bool


def expected_date_obs(path, cadence_seconds):
    """Calculate the expected timestamp from the model path and file index."""
    if not math.isfinite(cadence_seconds) or cadence_seconds <= 0:
        raise ValueError('Cadence must be finite and greater than zero')
    match = RADIANCE_MAP_PATTERN.fullmatch(path.name)
    if match is None:
        raise ValueError('Could not determine frame index from {}'.format(path))
    frame_index = int(match.group(1))
    return model_frame_observation_time(
        str(path), frame_index, cadence_seconds * u.second)


def _encode(block):
    return base64.b64encode(block).decode('ascii')


def _decode(block):
    return base64.b64decode(block, validate=True)


def _read_backup(path):
    with gzip.open(path, 'rt', encoding='utf-8') as stream:
        return json.load(stream)


def _check_current_headers(stream, manifest):
    """Allow an original, completed, or interrupted repair; reject other edits."""
    if os.fstat(stream.fileno()).st_size != manifest['file_size']:
        raise ValueError('File size differs from its timestamp backup')
    for block in manifest['headers']:
        original, repaired = _decode(block['original']), _decode(block['repaired'])
        if len(original) != len(repaired) or len(original) % 2880:
            raise ValueError('Invalid header lengths in timestamp backup')
        stream.seek(block['offset'])
        if stream.read(len(original)) not in (original, repaired):
            raise ValueError('Header differs from both original and repaired backup')


def _write_headers(stream, manifest, version):
    for block in manifest['headers']:
        stream.seek(block['offset'])
        stream.write(_decode(block[version]))
    stream.flush()
    os.fsync(stream.fileno())


def restore_backup(backup_path):
    """Restore original header blocks after checking no other edits occurred."""
    manifest = _read_backup(backup_path)
    path = Path(manifest['path'])
    with path.open('r+b') as stream:
        _check_current_headers(stream, manifest)
        _write_headers(stream, manifest, 'original')
    return path


def repair_file(path, cadence_seconds=10.0, write=False):
    """Repair fixed-size header blocks only, with a compressed recovery backup.

    Image bytes and HDU offsets stay unchanged. Existing DATASUM values allow
    CHECKSUM updates without reading image arrays. If CHECKSUM exists without
    DATASUM, Astropy calculates and adds the data checksum from the original
    file first (provided the additional card fits in the same header block).
    """
    path = Path(path).resolve()
    expected = expected_date_obs(path, cadence_seconds)
    blocks = []
    with fits.open(path, mode='readonly', memmap=True) as hdul:
        previous_dates = tuple(hdu.header.get('DATE-OBS') for hdu in hdul)
        if write and any(value != expected for value in previous_dates):
            with path.open('rb') as stream:
                file_size = os.fstat(stream.fileno()).st_size
                for hdu in hdul:
                    info = hdu.fileinfo()
                    offset = info['hdrLoc']
                    length = info['datLoc'] - offset
                    stream.seek(offset)
                    original = stream.read(length)
                    if hdu.header.get('DATE-OBS') == expected:
                        continue
                    data_checksum = None
                    if 'CHECKSUM' in hdu.header:
                        if 'DATASUM' in hdu.header:
                            data_checksum = int(hdu.header['DATASUM'])
                            if not 0 <= data_checksum <= 0xffffffff:
                                raise ValueError('Invalid DATASUM in {}'.format(path))
                        else:
                            data_checksum = hdu._calculate_datasum()
                        # Uses only header bytes and the unchanged data checksum.
                        if hdu._calculate_checksum(data_checksum) != hdu.header['CHECKSUM']:
                            raise ValueError('Existing CHECKSUM is inconsistent in {}'.format(path))
                        if 'DATASUM' not in hdu.header:
                            hdu.header['DATASUM'] = str(data_checksum)
                    hdu.header['DATE-OBS'] = expected
                    if data_checksum is not None:
                        hdu.header['CHECKSUM'] = hdu._calculate_checksum(data_checksum)
                    repaired = hdu.header.tostring().encode('ascii')
                    if len(repaired) != length:
                        raise ValueError('DATE-OBS edit would change an HDU header size')
                    blocks.append({'offset': offset, 'original': _encode(original),
                                   'repaired': _encode(repaired)})

    needs_update = any(value != expected for value in previous_dates)
    if write and needs_update:
        manifest = {'path': str(path), 'file_size': file_size, 'headers': blocks}
        backup_path = path.with_name(path.name + BACKUP_SUFFIX)
        if backup_path.exists():
            saved = _read_backup(backup_path)
            if saved['path'] != manifest['path'] or saved['file_size'] != file_size:
                raise ValueError('Existing backup belongs to a different file')
            saved_blocks = {block['offset']: block for block in saved['headers']}
            if any(block['offset'] not in saved_blocks or
                   block['repaired'] != saved_blocks[block['offset']]['repaired']
                   for block in blocks):
                raise ValueError('Existing backup describes a different timestamp repair')
            manifest = saved
        else:
            # Exclusive creation preserves the first repair's original headers.
            with backup_path.open('xb') as backup_stream:
                with gzip.GzipFile(fileobj=backup_stream, mode='wb') as compressed:
                    compressed.write(json.dumps(manifest).encode('utf-8'))
                backup_stream.flush()
                os.fsync(backup_stream.fileno())
        with path.open('r+b') as stream:
            _check_current_headers(stream, manifest)
            try:
                _write_headers(stream, manifest, 'repaired')
            except BaseException:
                # Roll back ordinary failures; the backup also permits recovery
                # after process termination between individual header writes.
                _write_headers(stream, manifest, 'original')
                raise

    return RepairResult(path, expected, previous_dates, needs_update)


def find_radiance_maps(roots):
    """Find radiance-map FITS files under one or more files/directories."""
    files = set()
    for root in roots:
        root = Path(root).expanduser()
        if root.is_file():
            if not RADIANCE_MAP_PATTERN.fullmatch(root.name):
                raise ValueError('Not a standard radiance-map filename: {}'.format(root))
            files.add(root.resolve())
        elif root.is_dir():
            for path in root.rglob('radiance_maps_*.fits'):
                if RADIANCE_MAP_PATTERN.fullmatch(path.name):
                    files.add(path.resolve())
                else:
                    warnings.warn('Ignoring nonstandard radiance-map filename: {}'.format(path))
        else:
            raise FileNotFoundError(root)
    return sorted(files)


def _parser():
    parser = argparse.ArgumentParser(
        description=(
            'Check or repair DATE-OBS in rendered EUV radiance maps. The default '
            'is a dry run; pass --write to update fixed-size headers in place. '
            'Original headers are saved alongside each file as ' + BACKUP_SUFFIX + '.'
        )
    )
    parser.add_argument('roots', nargs='+', help='FITS files or directories to scan recursively')
    parser.add_argument('--cadence-seconds', type=float, default=10.0)
    actions = parser.add_mutually_exclusive_group()
    actions.add_argument('--write', action='store_true', help='Update mismatched FITS headers')
    actions.add_argument('--restore', action='store_true',
                         help='Restore original headers from the specified backup files')
    parser.add_argument('--verbose', action='store_true', help='Print every mismatched file')
    return parser


def main(argv=None):
    args = _parser().parse_args(argv)
    if args.restore:
        for backup_path in args.roots:
            print('Restored {}'.format(restore_backup(backup_path)))
        return 0
    try:
        if not math.isfinite(args.cadence_seconds) or args.cadence_seconds <= 0:
            raise ValueError('Cadence must be finite and greater than zero')
        paths = find_radiance_maps(args.roots)
    except (FileNotFoundError, ValueError) as error:
        print(str(error), file=sys.stderr)
        return 1

    mismatches = 0
    errors = 0
    for path in paths:
        try:
            result = repair_file(path, args.cadence_seconds, write=args.write)
        except (OSError, ValueError) as error:
            errors += 1
            print('Skipping {}: {}'.format(path, error), file=sys.stderr)
            continue

        if result.needs_update:
            mismatches += 1
            if args.verbose:
                action = 'updated' if args.write else 'would update'
                old_dates = sorted({str(value) for value in result.previous_dates})
                print('{}: {} {} -> {}'.format(
                    action, result.path, old_dates, result.expected_date_obs))

    action = 'Updated' if args.write else 'Would update'
    print('{} {} of {} files; {} errors.'.format(action, mismatches, len(paths), errors))
    if not args.write and mismatches:
        print('Dry run only. Re-run with --write to apply these changes.')
    return 1 if errors else 0


if __name__ == '__main__':
    raise SystemExit(main())
