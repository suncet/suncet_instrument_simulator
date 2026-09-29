"""Local-only exposure-normalized frame stretch comparison; no publishing path."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from PIL import Image, ImageDraw

from exploration_notebooks.colortable_options import options, font, display_arrays
from suncet_instrument_simulator import make_movie


def variants(data):
    maximum = float(np.max(data))
    low, high, width = make_movie._asinh_limits([data])
    _, _, _, _, original = display_arrays(data, {"LEVEL": 1, "BUNIT": "DN/s"})
    result = []

    def add(slug, title, pixels, lo, hi, detail):
        result.append(dict(slug=slug, title=title, pixels=np.clip(pixels, 0, 1),
                           vmin=float(lo), vmax=float(hi), detail=detail))

    for slug, title, lo in [("log43", "Log10 / your setting", 43.2777),
                            ("log10", "Log10 / low black level", 10.),
                            ("log100", "Log10 / darker background", 100.),
                            ("log200", "Log10 / corona emphasis", 200.)]:
        add(slug, title, np.log10(np.maximum(data, lo) / lo) / np.log10(maximum / lo),
            lo, maximum, "log10(I / vmin) / log10(vmax / vmin)")
    lo = 43.2777
    for softening in [10., 50., 200., 1000.]:
        add(f"asinh{int(softening)}", f"Asinh / softening {softening:g} DN/s",
            np.arcsinh(np.maximum(data - lo, 0) / softening) /
            np.arcsinh((maximum - lo) / softening), lo, maximum,
            f"asinh((I - vmin) / {softening:g}); normalized at vmax")
    for power, name in [(.5, "Square root"), (.25, "Fourth root"),
                        (.125, "Eighth root"), (.1, "Tenth root")]:
        add(f"power{power}", name, np.clip((data-lo)/(maximum-lo), 0, 1)**power,
            lo, maximum, f"((I - vmin) / (vmax - vmin)) ^ {power:g}")
    log = np.log10(np.maximum(data, lo)/lo) / np.log10(maximum/lo)
    add("loggamma", "Log10 + gamma 0.8", log**.8, lo, maximum,
        "Your log10 setting, then display gamma 0.8")
    for key, title in [("current", "Gallery fourth-root formula"),
                       ("asinh", "Gallery asinh formula")]:
        display, vmin, vmax = original[key]
        add(f"reference-{key}", title, (display-vmin)/(vmax-vmin), low, high,
            f"Existing gallery formula; asinh softening {width:.6g} DN/s" if key == "asinh"
            else "Existing gallery formula: (I^0.25 - low^0.25) / (high^0.25 - low^0.25)")
    return result


def render(source, output):
    output.mkdir(parents=True, exist_ok=True)
    raw, header = fits.getdata(source, header=True)
    if header.get("LEVEL") != 1 or header.get("BUNIT") != "DN/s":
        raise ValueError("Supply the exposure-normalized Level 1 DN/s FITS")
    palettes = [p for slug in ["sunset-gold", "current-inferno", "aurora-mint"]
                for p in options() if p["slug"] == slug]
    manifest = {"input": str(source.resolve()), "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                "orientation": "Upper origin, same as voting gallery (old GUI used lower origin)",
                "palettes": [{"slug": p["slug"], "title": p["title"]} for p in palettes], "modes": {}}
    for mode in ["unfiltered", "filtered"]:
        data = raw if mode == "unfiltered" else make_movie.apply_radial_filter(raw, 300)
        entries = variants(data)
        manifest["modes"][mode] = [{k: v for k, v in e.items() if k != "pixels"} for e in entries]
        for p in palettes:
            folder = output / mode / p["slug"]
            folder.mkdir(parents=True, exist_ok=True)
            for e in entries:
                pixels = p["cmap"](e["pixels"], bytes=True)
                Image.fromarray(pixels).convert("RGB").save(folder / f'{e["slug"]}.png')
                thumb = Image.fromarray(pixels).convert("RGB")
                thumb.thumbnail((500, 375), Image.Resampling.LANCZOS)
                thumb.save(folder / f'{e["slug"]}.webp', quality=92)
        sheet = Image.new("RGB", (1500, 5*435), "#101113")
        draw = ImageDraw.Draw(sheet)
        for i, e in enumerate(entries):
            x, y = i % 3 * 500, i // 3 * 435
            with Image.open(output / mode / "sunset-gold" / f'{e["slug"]}.webp') as im:
                sheet.paste(im, (x, y))
            draw.text((x+12, y+380), f'{i+1:02}  {e["title"]}', font=font(17), fill="white")
            draw.text((x+12, y+407), f'{e["vmin"]:.4f} to {e["vmax"]:.2f} DN/s', font=font(14), fill="#aeb5bd")
        sheet.save(output / f"overview-{mode}.jpg", quality=93)
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2))
    template = Path(__file__).with_name("stretch_options.html").read_text()
    (output / "index.html").write_text(template.replace("__DATA__", json.dumps(manifest)))
    print(json.dumps({"gallery": str((output / "index.html").resolve()), "stretches": 15,
                      "pngs": 90, "maximum_dn_per_second": float(raw.max())}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fits", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("output/stretch-options/frame300"))
    args = parser.parse_args()
    render(args.fits, args.output)
