"""Wrap synthetic SunCET images in pipeline-compatible Level 0.5 products."""

from __future__ import annotations

import argparse
import json
import os
import re
from importlib.resources import files
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np
from astropy.io import fits


TEMPLATE_RESOURCE = "data/csie_meta_template_v1_1_2.json"
ENCODING_VALUES = {"16": 0, "8": 1, "12": 2, "jpegls": 3}


def load_metadata_template(path: str | Path | None = None) -> dict[str, object]:
    """Load the bundled CSIE v1.1.2 schema or a pipeline-produced meta JSON."""
    if path is None:
        resource = files("suncet_instrument_simulator").joinpath(TEMPLATE_RESOURCE)
        with resource.open("r", encoding="utf-8") as stream:
            metadata = json.load(stream)
    else:
        with Path(path).open("r", encoding="utf-8") as stream:
            metadata = json.load(stream)
    if not isinstance(metadata, dict):
        raise ValueError("metadata template must contain a JSON object")
    metadata.pop("image_id", None)
    return metadata


def _read_image(path: Path) -> tuple[np.ndarray, fits.Header]:
    with fits.open(path, memmap=False) as hdul:
        for hdu in hdul:
            if hdu.data is not None and np.ndim(hdu.data) == 2:
                return np.array(hdu.data, copy=True), hdu.header.copy()
    raise ValueError(f"{path} does not contain a two-dimensional FITS image")


def _as_uint16(image: np.ndarray) -> np.ndarray:
    array = np.asanyarray(image)
    if not np.issubdtype(array.dtype, np.number):
        raise TypeError(f"image dtype must be numeric, got {array.dtype}")
    if not np.all(np.isfinite(array)):
        raise ValueError("image contains NaN or infinite values")
    rounded = np.rint(array)
    if rounded.size and (rounded.min() < 0 or rounded.max() > np.iinfo(np.uint16).max):
        raise ValueError("image values must fit in the Level 0.5 uint16 range 0..65535")
    return rounded.astype(np.uint16)


def build_metadata(
    image: np.ndarray,
    source_header: Mapping[str, object],
    *,
    image_id: int,
    template: Mapping[str, object],
    row_bin: int | None = None,
    col_bin: int | None = None,
    exposure_ms: int | None = None,
    encoding: str = "16",
) -> dict[str, object]:
    """Populate current CSIE metadata fields from an image and its simulator header."""
    rows, cols = image.shape
    if row_bin is None:
        row_bin = int(source_header.get("NBIN2", 1) or 1)
    if col_bin is None:
        col_bin = int(source_header.get("NBIN1", 1) or 1)
    if row_bin < 1 or col_bin < 1:
        raise ValueError("row and column bin factors must be positive")

    if exposure_ms is None:
        exposure_seconds = source_header.get("EXPTIME")
        if exposure_seconds is not None:
            exposure_ms = int(round(float(exposure_seconds) * 1000))
    if exposure_ms is not None and exposure_ms < 0:
        raise ValueError("exposure time cannot be negative")

    metadata = dict(template)
    metadata.update(
        {
            "PKT_APID": 538,
            "SEQ_CTR": int(image_id) & 0x3FFF,
            "csie_meta_img_id": int(image_id),
            "csie_meta_capture_id": int(image_id),
            "csie_meta_dfsu_start": 0,
            "csie_meta_dfsu_end": int(rows * row_bin - 1),
            "csie_meta_fpm_row_per_frame": int(rows * row_bin),
            "csie_meta_fpm_pix_per_row": int(cols * col_bin),
            "fpm_proc_cfg_row_bin_meta": int(row_bin - 1),
            "fpm_proc_cfg_col_bin_meta": int(col_bin - 1),
            "icm_proc_cfg_encoding_meta": ENCODING_VALUES[encoding],
            "csie_meta_fletcher32": 0,
        }
    )
    if exposure_ms is not None:
        metadata["csie_meta_intg_ms"] = int(exposure_ms)
    return metadata


def _fits_keyword(name: str, used_keys: set[str]) -> str | None:
    if name.startswith("csie_"):
        name = name[5:]
    if name.startswith("meta_"):
        name = name[5:]
    base_key = re.sub(r"[^A-Za-z0-9]", "", name).upper()[:8] or "CSIE"
    if base_key not in used_keys:
        return base_key
    for suffix in range(1, 10):
        candidate = f"{base_key[:7]}{suffix}"
        if candidate not in used_keys:
            return candidate
    return None


def write_level0_5_product(
    image: np.ndarray,
    output_fits: str | Path,
    *,
    image_id: int,
    metadata: Mapping[str, object],
) -> tuple[Path, Path]:
    """Write FITS and JSON with the current processing-pipeline Level 0.5 layout."""
    image = _as_uint16(image)
    output_fits = Path(output_fits)
    output_fits.parent.mkdir(parents=True, exist_ok=True)
    meta_path = output_fits.with_name(f"{output_fits.stem}_meta.json")

    hdu = fits.PrimaryHDU(image)
    header = hdu.header
    rows = int(image.shape[0])
    header["IMAGEID"] = int(image_id)
    header["ROWS"] = rows
    header["COLS"] = int(image.shape[1])
    header["ROWPKTS"] = rows
    header["MISROWS"] = 0
    header["DUPROWS"] = 0
    header["CHKFAIL"] = 0
    header["CHKMISS"] = 0
    header["SELVALID"] = rows
    header["SELFAIL"] = 0
    header["SELMISS"] = 0
    header["PARTIAL"] = False
    header["ZEROFILL"] = 0
    header["MISZERO"] = 0
    header["CHKZERO"] = 0
    header["OUTRANGE"] = 0

    used_keys = set(header)
    for attr_name, value in sorted(metadata.items()):
        if attr_name.startswith("_") or not isinstance(value, (str, int, float, bool)):
            continue
        key = _fits_keyword(attr_name, used_keys)
        if key is None:
            continue
        try:
            header[key] = value
        except (TypeError, ValueError):
            continue
        used_keys.add(key)

    hdu.writeto(output_fits, overwrite=True)
    payload = {"image_id": int(image_id), **dict(metadata)}
    with meta_path.open("w", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return output_fits, meta_path


def convert_images(
    inputs: Iterable[str | Path],
    output_dir: str | Path,
    *,
    image_id_start: int = 1,
    metadata_template: str | Path | None = None,
    row_bin: int | None = None,
    col_bin: int | None = None,
    exposure_ms: int | None = None,
    encoding: str = "16",
    suffix: str = "",
) -> list[tuple[Path, Path]]:
    """Convert one or more synthetic FITS images into Level 0.5 product pairs."""
    template = load_metadata_template(metadata_template)
    output_dir = Path(output_dir)
    products = []
    for offset, input_name in enumerate(inputs):
        input_path = Path(input_name)
        image, source_header = _read_image(input_path)
        image_id = image_id_start + offset
        metadata = build_metadata(
            image,
            source_header,
            image_id=image_id,
            template=template,
            row_bin=row_bin,
            col_bin=col_bin,
            exposure_ms=exposure_ms,
            encoding=encoding,
        )
        output_path = output_dir / f"image_{image_id}{suffix}.fits"
        products.append(
            write_level0_5_product(
                image, output_path, image_id=image_id, metadata=metadata
            )
        )
    return products


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Wrap synthetic FITS images in realistic SunCET Level 0.5 FITS/JSON pairs."
    )
    parser.add_argument("inputs", nargs="+", help="Input FITS image(s), in output order")
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Defaults to $suncet_data/synthetic/level0_5",
    )
    parser.add_argument("--image-id-start", type=int, default=1)
    parser.add_argument("--metadata-template", type=Path)
    parser.add_argument("--row-bin", type=int, help="Physical rows per output row; defaults to NBIN2")
    parser.add_argument("--col-bin", type=int, help="Physical columns per output column; defaults to NBIN1")
    parser.add_argument("--exposure-ms", type=int, help="Defaults to EXPTIME converted to milliseconds")
    parser.add_argument("--encoding", choices=sorted(ENCODING_VALUES), default="16")
    parser.add_argument("--suffix", default="", help="Optional filename suffix, e.g. -synthetic")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _parser()
    args = parser.parse_args(argv)
    output_dir = args.output_dir
    if output_dir is None:
        data_root = os.getenv("suncet_data")
        if not data_root:
            parser.error("--output-dir is required when suncet_data is not set")
        output_dir = Path(data_root) / "synthetic" / "level0_5"
    products = convert_images(
        args.inputs,
        output_dir,
        image_id_start=args.image_id_start,
        metadata_template=args.metadata_template,
        row_bin=args.row_bin,
        col_bin=args.col_bin,
        exposure_ms=args.exposure_ms,
        encoding=args.encoding,
        suffix=args.suffix,
    )
    for fits_path, json_path in products:
        print(f"Wrote {fits_path}")
        print(f"Wrote {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
