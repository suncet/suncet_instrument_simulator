import numpy as np

from exploration_notebooks.colortable_options import display_arrays, versioned_asset, voting_display_arrays


def test_voting_stretches_match_selected_local_options():
    from exploration_notebooks.stretch_options import variants
    data = np.geomspace(1, 936214.2857142857, 1000).reshape(20, 50)
    low, high, width, displays = voting_display_arrays(data)
    assert (low, high, width) == (43.2777, float(data.max()), 200.)
    experiment = variants(data)
    np.testing.assert_array_equal(displays["current"][0], experiment[0]["pixels"])
    np.testing.assert_array_equal(displays["asinh"][0], experiment[6]["pixels"])


def test_normalized_frame_uses_shared_percentile_limits():
    data = np.arange(1, 101, dtype=float).reshape(10, 10) * 10000
    normalized, low, high, _, displays = display_arrays(data, {"LEVEL": 1, "BUNIT": "DN/s"})
    assert normalized
    assert displays["current"][1:] == (low ** .25, high ** .25)
    assert displays["current"][2] > 21


def test_legacy_header_does_not_accidentally_enable_normalization():
    data = np.ones((10, 10))
    normalized, _, _, _, displays = display_arrays(data, {"LEVEL": .5, "BUNIT": "dn/sec"})
    assert not normalized
    assert displays["current"][1:] == (.08, 21.)


def test_asset_version_changes_with_content_not_filename(tmp_path):
    asset = tmp_path / "frame.png"
    asset.write_bytes(b"before")
    before = versioned_asset(tmp_path, asset.name)
    asset.write_bytes(b"after")
    after = versioned_asset(tmp_path, asset.name)
    assert before != after
    assert before.startswith("frame.png?v=")
    assert after.startswith("frame.png?v=")
