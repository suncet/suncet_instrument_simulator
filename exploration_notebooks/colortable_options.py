"""Render a reproducible color-table study without changing the movie defaults.

Run from the repository root:
    python -m exploration_notebooks.colortable_options --fits /path/to/frame.fits
"""

import argparse
import hashlib
import json
from pathlib import Path
import shutil

import matplotlib
import numpy as np
from astropy.io import fits
from matplotlib.colors import LinearSegmentedColormap
from PIL import Image, ImageDraw, ImageFont
import sunpy
from sunpy.visualization.colormaps import cm

from suncet_instrument_simulator import make_movie


SUNPY_SOURCE = "https://docs.sunpy.org/en/stable/_modules/sunpy/visualization/colormaps/cm.html"
BRAND = {
    "orange": "#ed8a00", "blue": "#0063ed", "yellow": "#daed00",
    "red": "#ed1300", "green": "#00ed8a", "violet": "#8a00ed",
}
# Exact, frequently occurring pixels measured from SunCET Poster v3.png.
POSTER = {
    "night": "#141234", "indigo": "#393b71", "violet": "#6a5391",
    "lavender": "#6969a2", "lilac": "#c276c3", "rose": "#ee97b7",
    "pink": "#f7b7c4", "highlight": "#fbdcde",
    "blue": "#339ce2",
}


def options():
    result = []

    def add(slug, title, family, description, source=None, colors=None):
        cmap = (LinearSegmentedColormap.from_list(slug, colors, N=256)
                if colors else matplotlib.colormaps[source])
        result.append(dict(id=len(result) + 1, slug=slug, title=title,
                           family=family, description=description,
                           source=source, anchors=colors, cmap=cmap))

    add("current-inferno", "Current / Inferno", "Reference",
        "The current make_movie.py color table.", "inferno")
    heritage = [
        ("aia171-eui174", "AIA 171 / EUI 174", "sdoaia171",
         "Standard SunPy table; also SUVI 171 and EUI FSI/HRI 174."),
        ("aia193-suvi195", "AIA 193 / SUVI 195", "goes-rsuvi195",
         "Standard SunPy table, shared by AIA 193 and SUVI 195."),
        ("aia211", "AIA 211", "sdoaia211", "Standard AIA violet table."),
        ("aia335-suvi284", "AIA 335 / SUVI 284", "goes-rsuvi284",
         "Standard SunPy table, shared by AIA 335 and SUVI 284."),
        ("aia304-eui304", "AIA 304 / EUI 304", "solar orbiterfsi304",
         "Standard SunPy table; also SUVI 304."),
        ("suvi131", "SUVI 131 / AIA 131", "goes-rsuvi131", "Standard teal mission table."),
        ("suvi94", "SUVI 94 / AIA 94", "goes-rsuvi94", "Standard green mission table."),
        ("eit171", "SOHO / EIT 171", "sohoeit171", "Standard SOHO/EIT table."),
        ("eit195", "SOHO / EIT 195", "sohoeit195", "Standard SOHO/EIT table."),
        ("eit284", "SOHO / EIT 284", "sohoeit284", "Standard SOHO/EIT table."),
    ]
    for slug, title, source, description in heritage:
        add(slug, title, "Missions", description, source)
    add("eui-amber", "EUI-inspired amber", "Missions",
        "Custom warm-gold variation; not an official EUI lookup table.",
        colors=["#000000", "#251003", "#63320b", "#af681d", "#e8ae44", "#ffdf87", "#fff7de"])
    custom = [
        ("sunset-gold", "Sunset / Golden hour", "Sunset",
         ["#05030a", "#26102c", "#65244e", "#b6444b", "#ee7950", "#ffbd69", "#fff2c9"]),
        ("sunset-rose", "Sunset / Rose", "Sunset",
         ["#03030b", "#24163f", "#653167", "#a94b83", "#e9799b", "#ffb3b2", "#fff0da"]),
        ("sunset-twilight", "Sunset / Twilight", "Sunset",
         ["#020309", "#151b43", "#493776", "#945784", "#e08483", "#ffc58a", "#fff4d5"]),
        ("sunset-ember", "Sunset / Ember", "Sunset",
         ["#000000", "#25070f", "#66121f", "#ad2b24", "#ed602b", "#ffa953", "#ffe8b0"]),
        ("sunset-coral", "Sunset / Coral", "Sunset",
         ["#030409", "#16263d", "#563f61", "#a6617e", "#e68e93", "#ffc4ad", "#fff4dc"]),
        ("euv-violet", "EUV / Electric violet", "Violet",
         ["#000005", "#16082e", "#3c126e", "#702cc0", "#a461ef", "#d5acff", "#faf4ff"]),
        ("euv-amethyst", "EUV / Amethyst", "Violet",
         ["#020106", "#1e102c", "#4c275f", "#7e438f", "#ad6db8", "#d4a8dc", "#f8edf9"]),
        ("euv-ice", "EUV / Violet ice", "Violet",
         ["#000005", "#180c38", "#43277d", "#6652bf", "#8596df", "#afd9f1", "#efffff"]),
        ("euv-lilac", "EUV / Lilac", "Violet",
         ["#030107", "#221331", "#543266", "#86588f", "#b987b6", "#dfbade", "#fff2fc"]),
        ("euv-magenta", "EUV / Violet to pink", "Violet",
         ["#000005", "#220833", "#551065", "#90228f", "#ca56ae", "#ed9bcf", "#fff0f9"]),
        ("suncet-brand-fire", "SunCET / Brand fire", "SunCET",
         ["#000000", "#210837", BRAND["violet"], BRAND["red"], BRAND["orange"], BRAND["yellow"], "#fff8dd"]),
        ("suncet-blue-gold", "SunCET / Blue + gold", "SunCET",
         [(0, "#000000"), (.16, "#081d45"), (.36, BRAND["blue"]),
          (.53, "#ad89a0"), (.70, BRAND["orange"]), (.88, "#ffcf73"), (1, "#fff5df")]),
        ("suncet-triad", "SunCET / Violet + green", "SunCET",
         [(0, "#000000"), (.16, "#1c103e"), (.36, BRAND["violet"]),
          (.55, "#6c91ca"), (.76, BRAND["green"]), (.91, "#bdedaf"), (1, "#f3ffe9")]),
        ("poster-rose", "Poster / Lilac + rose", "Poster",
         ["#010003", POSTER["night"], POSTER["violet"], POSTER["lilac"],
          POSTER["rose"], POSTER["pink"], POSTER["highlight"]]),
        ("poster-indigo", "Poster / Indigo glow", "Poster",
         ["#010003", POSTER["night"], POSTER["indigo"], POSTER["lavender"],
          POSTER["lilac"], POSTER["pink"], POSTER["highlight"]]),
    ]
    for slug, title, family, colors in custom:
        description = {
            "Sunset": "Custom sequential sunset colors.",
            "Violet": "Custom violet false-color table for EUV imagery.",
            "SunCET": "Uses printed hex anchors from the SunCET palette PDF, with dark/bright extensions.",
            "Poster": "Uses exact indigo, lilac, and pink colors sampled from SunCET Poster v3.png.",
        }[family]
        add(slug, title, family, description, colors=colors)
    for name in ["magma", "cividis", "viridis", "gray"]:
        add(name, "Grayscale" if name == "gray" else name.title(), "Alternatives",
            "Matplotlib reference table.", name)
    add("coronal-ice", "Coronal ice", "Alternatives", "Custom cool blue-to-cyan table.",
        colors=["#000005", "#0b1730", "#193c60", "#266e91", "#4ca3b3", "#9ed7d2", "#f0fff2"])
    # Append new candidates so existing favorite IDs and URLs stay stable.
    add("poster-blue-rose", "Poster / Blue to rose", "Poster",
        "Poster purple shadows, bright blue loops, and pink highlights, with a pale extension.",
        colors=[(0, "#010003"), (.16, POSTER["night"]), (.32, POSTER["violet"]),
                (.52, POSTER["blue"]), (.69, "#b69bdc"), (.84, POSTER["rose"]),
                (.94, POSTER["pink"]), (1, POSTER["highlight"])])
    add("poster-blue-lilac", "Poster / Blue + lilac", "Poster",
        "A cooler poster variation: indigo and blue through lilac to pale pink.",
        colors=[(0, "#010003"), (.13, POSTER["night"]), (.27, "#3d316b"),
                (.45, "#246bb0"), (.60, POSTER["blue"]), (.76, "#bfa0df"),
                (.9, POSTER["pink"]), (1, POSTER["highlight"])])
    add("aurora-mint", "Aurora / Mint", "Alternatives",
        "Violet-black through cool teal and luminous mint to pale gold.",
        colors=["#030208", "#201739", "#1c4560", "#217d83", "#65b798", "#bce3b7", "#fff5d1"])
    add("pearl-opal", "Pearl / Opal", "Alternatives",
        "Restrained cool gray, blue, and rose with pearlescent highlights.",
        colors=["#020306", "#18202c", "#3b5263", "#75879b", "#b8afc8", "#e8d2db", "#fff8f0"])
    add("copper-teal", "Copper / Teal", "Alternatives",
        "Deep teal shadows transitioning through muted copper to warm gold.",
        colors=["#000407", "#082c3c", "#26515a", "#776762", "#bd886d", "#e9be8a", "#fff1cd"])
    add("tequila-sunrise", "Tequila Sunrise", "Sunset",
        "Grenadine red through orange juice gold to a pale citrus highlight, inspired by a Tequila Sunrise drink.",
        colors=["#080103", "#490817", "#a41929", "#e94725", "#fa8b24", "#ffd354", "#fff4bb"])
    for wavelength in [171, 195, 284, 304]:
        add(f"euvi{wavelength}", f"STEREO / EUVI {wavelength}", "Missions",
            f"Standard SunPy STEREO/EUVI {wavelength} angstrom color table.", f"euvi{wavelength}")
    families = {"Missions": "existing EUV imager tables", "Sunset": "inspired by sunsets",
                "Violet": "(Almost) ultraviolet", "SunCET": "SunCET branding",
                "Poster": "SunCET NASA Poster"}
    for item in result:
        if item["family"] == "SunCET":
            item["title"] = item["title"].replace("SunCET /", "SunCET branding /")
        item["family"] = families.get(item["family"], item["family"])
    # Retain historical numbers and slugs so existing votes never change meaning.
    return [item for item in result if item["id"] not in {24, 29}]


def font(size):
    return ImageFont.truetype(str(Path(matplotlib.get_data_path()) / "fonts/ttf/DejaVuSans.ttf"), size)


def contact_sheet(items, output, stretch, *, columns=4, thumb_width=400, title=None):
    gap, margin = 16, 24
    with Image.open(output / items[0]["images"][stretch]) as first:
        thumb_height = round(thumb_width * first.height / first.width)
    tile_height = thumb_height + 64
    rows = (len(items) + columns - 1) // columns
    sheet = Image.new("RGB", (2 * margin + columns * thumb_width + (columns - 1) * gap,
                              108 + rows * tile_height + (rows - 1) * gap + 20), "#101113")
    draw = ImageDraw.Draw(sheet)
    draw.text((margin, 18), title or "SunCET / Frame 300 / Color table study", font=font(28), fill="white")
    label = "Log10" if stretch == "current" else "Asinh / 200 DN/s softening"
    draw.text((margin, 60), label + "  |  Shared limits / No radial filter", font=font(17), fill="#aeb5bd")
    for i, item in enumerate(items):
        x = margin + (i % columns) * (thumb_width + gap)
        y = 108 + (i // columns) * (tile_height + gap)
        with Image.open(output / item["images"][stretch]) as im:
            sheet.paste(im.convert("RGB").resize((thumb_width, thumb_height), Image.Resampling.LANCZOS), (x, y))
        draw.text((x, y + thumb_height + 10), f'{item["id"]:02d}  {item["title"]}', font=font(18), fill="white")
        ramp = np.asarray(item["cmap"](np.linspace(0, 1, thumb_width), bytes=True))[:, :3]
        sheet.paste(Image.fromarray(np.repeat(ramp[None], 10, axis=0)), (x, y + thumb_height + 40))
    return sheet


def display_arrays(data, header):
    """Keep legacy DN rendering, but derive shared limits for normalized DN/s."""
    normalized = float(header.get("LEVEL", 0)) == 1 and header.get("BUNIT") == "DN/s"
    low, high, width = make_movie._asinh_limits([data])
    root_limits = (max(low, 0.) ** .25, high ** .25) if normalized else (.08, 21.)
    return normalized, low, high, width, {
        "current": (data ** .25, *root_limits),
        "asinh": (np.arcsinh((data.astype(np.float64) - low) / width),
                  0., float(np.arcsinh((high - low) / width))),
    }


def versioned_asset(output, relative):
    digest = hashlib.sha256((output / relative).read_bytes()).hexdigest()[:12]
    return f"{relative}?v={digest}"


def voting_display_arrays(data):
    """Local experiment options 1 and 7, retaining the stable asset directory keys."""
    low, high, width = 43.2777, float(np.max(data)), 200.
    if not np.all(np.isfinite(data)) or high <= low:
        raise ValueError("Voting input must be finite with maximum above 43.2777 DN/s")
    return low, high, width, {
        "current": (np.log10(np.maximum(data, low) / low) / np.log10(high / low), 0., 1.),
        "asinh": (np.arcsinh(np.maximum(data - low, 0) / width) /
                  np.arcsinh((high - low) / width), 0., 1.),
    }


def render(fits_path, output, public_output=None, poster_path=None):
    output.mkdir(parents=True, exist_ok=True)
    for directory in ["current", "asinh", "luts", "ramps", "sheets", "thumbs/current", "thumbs/asinh"]:
        (output / directory).mkdir(parents=True, exist_ok=True)
    raw, header = fits.getdata(fits_path, header=True)
    if header.get("LEVEL") != 1 or header.get("BUNIT") != "DN/s":
        raise ValueError("Voting gallery requires an exposure-normalized Level 1 DN/s FITS")
    data = np.asarray(raw, dtype=np.float64)
    normalized = True
    low, high, width, displays = voting_display_arrays(data)
    items = options()
    if poster_path:
        (output / "references").mkdir(exist_ok=True)
        with Image.open(poster_path) as poster:
            poster = poster.convert("RGB")
            poster.thumbnail((1200, 1800), Image.Resampling.LANCZOS)
            poster.save(output / "references/suncet-nasa-poster.webp", quality=90)
    for item in items:
        name = f'{item["id"]:02d}-{item["slug"]}'
        item["images"] = {}
        item["thumbnails"] = {}
        for stretch, (display, vmin, vmax) in displays.items():
            relative = f"{stretch}/{name}.png"
            pixels = make_movie._render_frame(display, vmin=vmin, vmax=vmax, cmap=item["cmap"],
                                              output_filename=output / relative)
            assert pixels.shape == (*data.shape, 4)
            assert np.all(pixels[..., 3] == 255)
            if item["id"] == 1:
                reference = make_movie._render_frame(display, vmin=vmin, vmax=vmax, cmap="inferno")
                np.testing.assert_array_equal(pixels, reference)
            item["images"][stretch] = relative
            thumbnail = f"thumbs/{stretch}/{name}.webp"
            thumb = Image.fromarray(pixels).convert("RGB")
            thumb.thumbnail((440, 330), Image.Resampling.LANCZOS)
            thumb.save(output / thumbnail, quality=88)
            item["thumbnails"][stretch] = thumbnail
        rgb = item["cmap"](np.linspace(0, 1, 256), bytes=True)[:, :3]
        np.savetxt(output / "luts" / f"{name}.csv", rgb, delimiter=",", fmt="%d", header="red,green,blue", comments="")
        Image.fromarray(np.repeat(rgb[None], 16, axis=0)).save(output / "ramps" / f"{name}.png")
        item["ramp"] = f"ramps/{name}.png"
        item["lut"] = f"luts/{name}.csv"
    # Explicit Open Graph artwork prevents crawlers from selecting the poster.
    social = Image.new("RGB", (1200, 630), "#101113")
    draw = ImageDraw.Draw(social)
    draw.text((24, 20), "SunCET | Color Voting", font=font(36), fill="white")
    for x, slug in [(24, "current-inferno"), (612, "aurora-mint")]:
        item = next(item for item in items if item["slug"] == slug)
        with Image.open(output / item["images"]["asinh"]) as picture:
            social.paste(picture.convert("RGB").resize((564, 423), Image.Resampling.LANCZOS), (x, 96))
        draw.text((x, 545), item["title"], font=font(28), fill="white")
    social.save(output / "social-preview-v1.png")
    for stretch in displays:
        contact_sheet(items, output, stretch).save(output / f"overview-{stretch}.png")
        for family in dict.fromkeys(x["family"] for x in items if x["family"] != "Reference"):
            subset = [items[0]] + [x for x in items if x["family"] == family]
            contact_sheet(subset, output, stretch, columns=3,
                          title=f"SunCET / Frame 300 / {family}").save(output / "sheets" / f"{family.lower()}-{stretch}.png")
    shortlist = [item for item in items if item["id"] in [1, 4, 13, 15, 18, 23, 26, 32]]
    contact_sheet(shortlist, output, "asinh", title="SunCET / Eight directions / Frame 300").save(output / "highlights.png")
    for stretch in displays:
        contact_sheet([items[0]] + [item for item in items if item["id"] >= 33], output, stretch, columns=3,
                      title="SunCET / New directions / Frame 300").save(output / f"new-options-{stretch}.png")
    manifest = {
        "study_id": "suncet-frame300-v1",
        "input": str(fits_path.resolve()), "input_sha256": hashlib.sha256(fits_path.read_bytes()).hexdigest(),
        "frame": 300, "shape": list(raw.shape), "date_obs": header.get("DATE-OBS"),
        "rendering": {"radial_filter_sigma_pixels": None, "origin": "upper (current renderer)",
                      "current": {"stretch": "log10", "vmin": low, "vmax": high,
                                  "formula": "log10(max(I, vmin)/vmin) / log10(vmax/vmin)"},
                      "asinh": {"low": low, "high": high, "width": width,
                                "formula": "asinh(max(I-low, 0)/width) / asinh((high-low)/width)",
                                "limits": "43.2777 DN/s to image maximum; shared by every option"}},
        "versions": {"numpy": np.__version__, "matplotlib": matplotlib.__version__, "sunpy": sunpy.__version__},
        "brand_hex": BRAND, "poster_hex": POSTER, "mission_source": SUNPY_SOURCE,
        "options": [{k: v for k, v in item.items() if k != "cmap"} for item in items],
    }
    manifest["processing"] = {
        "exposure_normalized": normalized,
        "fits_metadata": {key: header.get(key) for key in [
            "LEVEL", "BUNIT", "PROCSTAT", "EXP_MASK", "EFFEXPI", "EFFEXPO",
            "INTTIMEI", "INTTIMEO", "NSTACKI", "NSTACKO", "STKNORMI", "STKNORMO"]},
        "note": ("Exposure normalization only; not a complete Level 1 calibration. "
                 "Saturated source pixels remain unrecoverable." if normalized else "Legacy stored-DN composite."),
    }
    # Version URLs, not palette identities: existing favorites still refer to the same slugs.
    for item in manifest["options"]:
        for key in ["images", "thumbnails"]:
            item[key] = {stretch: versioned_asset(output, relative)
                         for stretch, relative in item[key].items()}
    manifest["overviews"] = {stretch: versioned_asset(output, f"overview-{stretch}.png") for stretch in displays}
    manifest["social_image"] = versioned_asset(output, "social-preview-v1.png")
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    web_source = Path(__file__).resolve().parents[1] / "web/colortables"
    script_version = hashlib.sha256((web_source / "gallery.js").read_bytes()).hexdigest()[:12]
    template = Path(__file__).with_name("colortable_gallery.html").read_text().replace("__GALLERY_VERSION__", script_version)
    template = template.replace("__SOCIAL_IMAGE__", manifest["social_image"])
    (output / "index.html").write_text(template.replace("__STUDY_JSON__", json.dumps(manifest).replace("</", "<\\/")))
    shutil.copyfile(web_source / "gallery.js", output / "gallery.js")
    if not (output / "config.js").exists():
        shutil.copyfile(web_source / "config.js", output / "config.js")
    shutil.copytree(web_source / "vendor", output / "vendor", dirs_exist_ok=True)
    (output / "README.md").write_text(
        "# SunCET color table study\n\n"
        f"Open index.html locally. There are {len(items)} distinct tables, each rendered at {raw.shape[1]} x {raw.shape[0]} pixels "
        "with two shared stretches. The gallery starts with log10.\n\n"
        f"Input: `{fits_path.resolve()}`\n\n" +
        ("The input is the pipeline's provisional exposure-normalized Level 1 frame (DN/s). "
         "Both stretches use 43.2777 DN/s to the image maximum, with no radial filter. "
         "Effective exposures are recorded in manifest.json; this is not full Level 1 calibration. "
         if normalized else "The current reference is pixel-identical to make_movie.plot_scaled_image(scale='1/4'). ") +
        "Images retain the existing orientation and borderless renderer. "
        "Log10 and asinh (200 DN/s softening) match local experiment options 1 and 7. "
        "The directory key 'current' is retained for compatibility and now means log10.\n\n"
        "Standard mission tables come from SunPy. Shared mission tables are combined into one option: "
        "AIA171/SUVI171/EUI174, AIA193/SUVI195, AIA335/SUVI284, AIA304/SUVI304/EUI304, "
        "AIA131/SUVI131, and AIA94/SUVI94. EUI-inspired amber is custom. "
        "Instrument names describe color tables, not a change in the simulated passband.\n\n"
        f"Source: {SUNPY_SOURCE}\n\n"
        "SunCET brand anchors use the PDF's printed hex values (some printed RGB numbers differ). "
        "Poster anchors are exact frequently occurring RGB pixels from SunCET Poster v3.png. "
        "Custom tables are interpolated in sRGB and are aesthetic candidates, not asserted to be perceptually uniform.\n\n"
        "The manifest records color anchors and rendering settings; luts/ contains reusable 256-row RGB CSVs. "
        "Load a CSV with matplotlib.colors.ListedColormap(np.loadtxt(path, delimiter=',', skiprows=1)/255). "
        "The CSVs reproduce the exported 8-bit color mapping. No production defaults were changed.\n"
    )
    if public_output:
        public_output.mkdir(parents=True, exist_ok=True)
        # Remove only known retired generated assets, never unrelated output files.
        for root in [output, public_output]:
            for name in ["24-suncet-blue-gold", "29-cividis"]:
                for directory, extension in [("current", "png"), ("asinh", "png"),
                                             ("ramps", "png"), ("luts", "csv"),
                                             ("thumbs/current", "webp"), ("thumbs/asinh", "webp")]:
                    (root / directory / f"{name}.{extension}").unlink(missing_ok=True)
        for directory in ["current", "asinh", "luts", "ramps", "thumbs", "vendor"]:
            shutil.copytree(output / directory, public_output / directory, dirs_exist_ok=True)
        for filename in ["gallery.js", "overview-current.png", "overview-asinh.png", "social-preview-v1.png"]:
            shutil.copyfile(output / filename, public_output / filename)
        if (output / "references").exists():
            shutil.copytree(output / "references", public_output / "references", dirs_exist_ok=True)
        if not (public_output / "config.js").exists():
            shutil.copyfile(web_source / "config.js", public_output / "config.js")
        public_manifest = {**manifest, "input": fits_path.name}
        (public_output / "manifest.json").write_text(json.dumps(public_manifest, indent=2) + "\n")
        (public_output / "index.html").write_text(template.replace("__STUDY_JSON__", json.dumps(public_manifest).replace("</", "<\\/")))
        (public_output / "README.md").write_text((output / "README.md").read_text().replace(str(fits_path.resolve()), fits_path.name))
        (public_output / ".nojekyll").touch()
        seed = "-- Generated by exploration_notebooks.colortable_options. Safe to rerun.\n"
        seed += "insert into public.colortable_studies (id, title) values ('suncet-frame300-v1', 'SunCET frame 300') on conflict (id) do nothing;\n"
        seed += "update public.colortable_palettes set is_active=false where study_id='suncet-frame300-v1' and slug in ('suncet-blue-gold', 'cividis');\n"
        for item in items:
            title = item["title"].replace("'", "''")
            seed += ("insert into public.colortable_palettes (study_id, slug, number, title) "
                     f"values ('suncet-frame300-v1', '{item['slug']}', {item['id']}, '{title}') "
                     "on conflict (study_id, slug) do update set number=excluded.number, title=excluded.title, is_active=true;\n")
        (web_source / "catalog.sql").write_text(seed)
    print(json.dumps({"gallery": str((output / "index.html").resolve()), "tables": len(items),
                      "full_resolution_pngs": 2 * len(items), "reference_pixel_match": True}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fits", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("output/colortables/frame300"))
    parser.add_argument("--public-output", type=Path, help="Publishable static files with local paths removed")
    parser.add_argument("--poster", type=Path, help="SunCET NASA poster reference image")
    args = parser.parse_args()
    render(args.fits, args.output, args.public_output, args.poster)


if __name__ == "__main__":
    main()
