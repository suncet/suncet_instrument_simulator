"""One-off audit/repair of output times and filenames, using cached radiance maps.

Run with one or more config paths. Defaults to a dry run. --write backs up every
changed file, preserves already-correct outputs on name collisions, and verifies
that image bytes and FITS checksums survive the repair.
"""

import argparse
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import os
from pathlib import Path
import re
import shutil
import tempfile

import astropy.units as u
from astropy.io import fits
from astropy.time import Time

from suncet_instrument_simulator.config_parser import Config
from suncet_instrument_simulator.make_radiance_maps import model_frame_observation_time
from suncet_instrument_simulator.stack_schedule import build_observation_sequence


@dataclass
class OutputRepair:
    source: Path
    destination: Path
    fields: dict
    action: str


def file_state(path):
    stat = path.stat()
    return stat.st_ino, stat.st_size, stat.st_mtime_ns


def pixel_hashes(hdul):
    return [None if hdu.data is None else hashlib.sha256(
        hdu.data.tobytes()).hexdigest() for hdu in hdul]


def validate_file(path, fields=None):
    with fits.open(path, memmap=False) as hdul:
        for hdu in hdul:
            if 'CHECKSUM' in hdu.header and hdu.verify_checksum() != 1:
                raise ValueError('Invalid CHECKSUM: {}'.format(path))
            if 'DATASUM' in hdu.header and hdu.verify_datasum() != 1:
                raise ValueError('Invalid DATASUM: {}'.format(path))
        matches = fields is None or all(
            hdul[0].header.get(key) == value for key, value in fields.items())
        return matches, pixel_hashes(hdul)


def plan_repairs(config_filename, data_root):
    config_filename = Path(config_filename)
    config = Config(str(config_filename))
    output_root = data_root / 'synthetic' / 'level0' / 'fits'
    radiance_root = (data_root / config.model_data_folder.strip('/')
                     / config.euv_radiance_map_directory_name.strip('/'))
    observations = {obs.output_index: obs for obs in build_observation_sequence(config)}
    pattern = re.compile(re.escape(config_filename.stem) + r'_OBS_.+_(\d+)\.fits')
    repairs = []
    for source in sorted(output_root.glob(config_filename.stem + '_OBS_*.fits')):
        match = pattern.fullmatch(source.name)
        if match is None or int(match[1]) not in observations:
            raise ValueError('Unrecognized output index: {}'.format(source))
        index = int(match[1])
        observation = observations[index]
        model_dt = config.model_timestep.to_value(u.s)
        model_index = int(observation.start_seconds // model_dt)
        radiance_file = radiance_root / 'radiance_maps_{:03d}.fits'.format(model_index)
        expected = model_frame_observation_time(
            config.model_directory_name, model_index, config.model_timestep)
        with fits.open(radiance_file) as hdul:
            for hdu in hdul:
                value = hdu.header.get('DATE-OBS')
                scale = str(hdu.header.get('TIMESYS', 'UTC')).strip().lower()
                if value is None or abs((Time(value, scale=scale).utc
                                         - Time(expected, scale='utc')).sec) > 0.0005:
                    raise ValueError('Repair cached radiance timestamps first: {}'.format(radiance_file))
            primary_scale = str(hdul[0].header.get('TIMESYS', 'UTC')).strip().lower()
            start = Time(hdul[0].header['DATE-OBS'], scale=primary_scale).utc
        start += (observation.start_seconds - model_index * model_dt) * u.s
        start.precision = 3
        destination = output_root / '{}_OBS_{}_{:03d}.fits'.format(
            config_filename.stem, start.isot, index)
        fields = {
            'DATE-OBS': start.isot, 'DATE-BEG': start.isot,
            'DATE-END': (start + config.observation_window).isot,
            'TIMESYS': 'UTC', 'FILENAME': destination.name,
        }
        matches, _ = validate_file(source, fields)
        action = 'correct' if matches and destination == source else 'repair'
        if destination != source and destination.exists():
            if not validate_file(destination, fields)[0]:
                raise ValueError('Destination exists but has incorrect metadata: {}'.format(destination))
            action = 'archive_duplicate'
        repairs.append(OutputRepair(source, destination, fields, action))
    return repairs


def apply_repair(repair, backup_dir):
    if repair.action == 'correct':
        return
    source = repair.source
    original_state = file_state(source)
    _, before = validate_file(source)
    backup_dir.mkdir(parents=True, exist_ok=True)
    backup = backup_dir / source.name
    if backup.exists():
        raise FileExistsError(backup)
    shutil.copy2(source, backup)
    if validate_file(backup)[1] != before or file_state(source) != original_state:
        raise RuntimeError('Source changed during backup: {}'.format(source))
    if repair.action == 'archive_duplicate':
        if not validate_file(repair.destination, repair.fields)[0]:
            raise ValueError('Destination changed: {}'.format(repair.destination))
        source.unlink()
        return

    descriptor, temporary_name = tempfile.mkstemp(prefix='.timestamp-repair-', suffix='.fits', dir=source.parent)
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        shutil.copy2(source, temporary)
        with fits.open(temporary, mode='update', memmap=False) as hdul:
            hdul[0].header.update(repair.fields)
            for hdu in hdul:
                hdu.add_checksum()
        matches, after = validate_file(temporary, repair.fields)
        if not matches or before != after or file_state(source) != original_state:
            raise RuntimeError('Repair verification failed: {}'.format(source))
        if repair.destination == source:
            os.replace(temporary, source)
        else:
            # Exclusive creation: a rerun must never overwrite a newer output.
            os.link(temporary, repair.destination)
            source.unlink()
    finally:
        temporary.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('configs', nargs='+', type=Path)
    parser.add_argument('--data-root', type=Path, default=os.environ.get('suncet_data'))
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    if args.data_root is None:
        parser.error('Set suncet_data or pass --data-root.')
    plans = []
    for config in args.configs:
        plan = plan_repairs(config, args.data_root)
        print(config.stem, dict(Counter(item.action for item in plan)), flush=True)
        plans.extend(plan)
    if args.write:
        backup_dir = (args.data_root / 'synthetic' / 'level0' / 'timestamp_repair_backups'
                      / datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
        for repair in plans:
            apply_repair(repair, backup_dir)
        print('Original files backed up in:', backup_dir)
        for config in args.configs:
            remaining = plan_repairs(config, args.data_root)
            if any(item.action != 'correct' for item in remaining):
                raise RuntimeError('Post-repair audit failed: {}'.format(config))
        print('Post-repair timestamps, filenames, checksums, and image data verified.')
    else:
        print('Dry run only. Use --write to apply repairs with backups.')


if __name__ == '__main__':
    main()
