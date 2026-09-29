import numpy as np

from exploration_notebooks.stretch_options import variants


def test_stretches_are_finite_bounded_monotonic():
    data = np.concatenate(([0.], np.geomspace(1, 1e6, 1000))).reshape(1, -1)
    entries = variants(data)
    assert len(entries) == 15
    assert len({e["slug"] for e in entries}) == 15
    for entry in entries:
        pixels = entry["pixels"]
        assert np.all(np.isfinite(pixels))
        assert pixels.min() >= 0
        assert pixels.max() <= 1
        assert np.all(np.diff(pixels[0]) >= 0)


def test_log_matches_prior_gui_physical_limits():
    low, high = 43.2777, 936214.2857142857
    data = np.array([[0, low, np.sqrt(low * high), high]])
    entry = variants(data)[0]
    assert entry["vmin"] == low
    assert entry["vmax"] == high
    np.testing.assert_allclose(entry["pixels"], [[0, 0, .5, 1]])
