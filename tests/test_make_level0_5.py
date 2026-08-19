import json

import numpy as np
from astropy.io import fits

from suncet_instrument_simulator.make_level0_5 import (
    convert_images,
    load_metadata_template,
    main,
)


def test_convert_image_matches_pipeline_level0_5_layout(tmp_path):
    source = tmp_path / "synthetic.fits"
    source_header = fits.Header({"NBIN1": 2, "NBIN2": 3, "EXPTIME": 0.035})
    image = np.arange(20, dtype=np.uint16).reshape(4, 5)
    fits.PrimaryHDU(image, source_header).writeto(source)

    [(fits_path, json_path)] = convert_images(
        [source], tmp_path / "level0_5", image_id_start=4112, suffix="-synthetic"
    )

    assert fits_path.name == "image_4112-synthetic.fits"
    assert json_path.name == "image_4112-synthetic_meta.json"
    with fits.open(fits_path) as hdul:
        np.testing.assert_array_equal(hdul[0].data, image)
        header = hdul[0].header
        assert header["IMAGEID"] == 4112
        assert header["ROWS"] == 4
        assert header["COLS"] == 5
        assert header["ROWPKTS"] == 4
        assert header["MISROWS"] == 0
        assert header["PARTIAL"] is False
        assert header["PKTAPID"] == 538
        assert header["IMGID"] == 4112
        assert header["FPMROWPE"] == 12
        assert header["FPMPIXPE"] == 10
        assert header["INTGMS"] == 35

    payload = json.loads(json_path.read_text())
    assert set(payload) == {"image_id", *load_metadata_template()}
    assert payload["image_id"] == 4112
    assert payload["csie_meta_img_id"] == 4112
    assert payload["csie_meta_capture_id"] == 4112
    assert payload["csie_meta_fpm_row_per_frame"] == 12
    assert payload["csie_meta_fpm_pix_per_row"] == 10
    assert payload["fpm_proc_cfg_row_bin_meta"] == 2
    assert payload["fpm_proc_cfg_col_bin_meta"] == 1
    assert payload["csie_meta_intg_ms"] == 35
    assert payload["icm_proc_cfg_encoding_meta"] == 0


def test_real_meta_json_can_supply_housekeeping_defaults(tmp_path):
    source = tmp_path / "synthetic.fits"
    fits.PrimaryHDU(np.ones((2, 3), dtype=np.uint16)).writeto(source)
    template = load_metadata_template()
    template["csie_meta_bus_volt"] = 1234
    template["image_id"] = 999
    template_path = tmp_path / "real_meta.json"
    template_path.write_text(json.dumps(template))

    [(_, json_path)] = convert_images(
        [source],
        tmp_path / "out",
        image_id_start=7,
        metadata_template=template_path,
        exposure_ms=1000,
        encoding="jpegls",
    )

    payload = json.loads(json_path.read_text())
    assert payload["image_id"] == 7
    assert payload["csie_meta_img_id"] == 7
    assert payload["csie_meta_bus_volt"] == 1234
    assert payload["csie_meta_intg_ms"] == 1000
    assert payload["icm_proc_cfg_encoding_meta"] == 3


def test_cli_defaults_to_suncet_data(tmp_path, monkeypatch):
    source = tmp_path / "synthetic.fits"
    fits.PrimaryHDU(np.ones((2, 3), dtype=np.uint16)).writeto(source)
    monkeypatch.setenv("suncet_data", str(tmp_path))

    assert main([str(source), "--image-id-start", "42"]) == 0

    output_dir = tmp_path / "synthetic" / "level0_5"
    assert (output_dir / "image_42.fits").is_file()
    assert (output_dir / "image_42_meta.json").is_file()
