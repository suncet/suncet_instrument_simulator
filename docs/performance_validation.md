# Simulator performance validation — 2026-09-04

The performance work starts from commit `482c35a9d39ddd89a7f588bd6392a49f25589299`.
That commit was verified on `origin/main` with a clean working tree before edits;
this work is on `codex/simulator-performance`.

**Measured outcome:** a full limb observation fell from **459.1 to 96.4 seconds
(4.8× faster)** with the final image unchanged for the tested seed. Peak sampled
resident memory fell from **28.4 to 19.9 GiB (30% lower)**. Radiance generation
and PSF generation also passed separate comparisons against the old algorithms.

| Workload | Original | Optimized | Comparison |
|---|---:|---:|---|
| Full observation, default options | 459.1 s | 96.4 s | 13 noisy detector images and final image bitwise identical |
| Full observation, retaining the pure reference | 459.1 s | 108.3 s | Same noisy images; pure images differ only by small floating-point rounding |
| Real 1024² EM frame → 40 radiance planes | 5.48 s | 0.73 s | All output float32 values identical |
| 2000² PSF, 8,000 orders, 185 Å | 477.7 s | 0.657 s | Full composite, wings, encircled energy, and convolved scene identical |

These are individual measurements on a 64 GiB macOS machine, using the configured
Python 3.11.11 / NumPy 2.2.1 / SciPy 1.15.0 / Astropy 7.0.0 / SunPy 6.0.4
environment. Desktop activity and some independent validation overlapped the
runs. They are not isolated repeated throughput measurements; in particular,
the difference between the two optimized observation times should not all be
attributed to omitting the pure reference. Memory was sampled every 0.25 seconds;
RSS reflects macOS memory management rather than total allocated array bytes.

**Implementation.** Repeated integrations share deterministic optical results
and read-only expected-photon arrays. The standard nine short integrations and
four long integrations therefore require five optical scene calculations instead
of thirteen. Each member still receives independent metadata with its own
integration timestamp, and the random draw order and existing detector-noise
correlations are preserved. Distinct members are processed separately to bound
intermediate memory. Cache keys retain the exact ordered normalized model
contributions; weights are neither rounded nor algebraically rearranged.

FFT convolution retains the existing cropped kernels, centering, wrap boundary,
and normalization. Finite, normalizable inputs use SciPy FFTs with two workers
and omit Astropy's extra NaN-interpolation transforms. Nonfinite/masked inputs
and nonnormalizable kernels retain the original path. Cropped results own compact
real arrays rather than retaining the larger complex FFT workspace.

Electron conversion accumulates wavelength contributions into a 2-D buffer
using weights derived from SciPy's Simpson rule, including its even-sample and
nonuniform-grid policies. Compatible wavelength-map units are converted to the
first band's unit as before. Pure reference images are optional through
`run(retain_pure_reference=True)`. Particle rejection uses a running uint32 sum
and maximum, preserving digital overflow behavior. FITS headers are constructed
in memory.

Radiance generation processes 4096 spatial pixels per batch, retaining the
native float32 cast before the existing spectral slices and sums. It avoids the
1.56 GiB native-wavelength cube. Pre-binning emissivity was deliberately not used
because it changes the rounding stages. The production float64 emissivity table
matched exactly. Tests also cover float32 emissivity: all-float32 inputs can
change the last bits, bounded to two ULP in the tested nonnegative cases (the
specific 600-band regression differs by at most one ULP). The batched `np.dot`
path avoids spurious macOS BLAS status warnings; nonfinite products fall back to
`matmul` and genuine product/cast overflow diagnostics remain tested.

CPU PSF generation evaluates local Gaussian support extending beyond float64
underflow, with an extra pixel margin. Every omitted term already evaluates to
zero in the dense original. Order iteration, coordinate arithmetic,
accumulation, normalization, and final FFT composition remain the same. This is
not the earlier approximate tail-cutoff prototype. Nonfinite intensities and
unusual widths retain dense evaluation; the GPU path is unchanged.

Movie generation streams the same rendered RGBA pixels into the encoder, with
optional PNG export and reuse of the previous input frame. The optional SNR RMS
calculation uses compiled local correlation instead of per-pixel Python. It
retains the historical footprint, edge counts, integer overflow, and local
nonfinite behavior; float32 rounding is tested within four ULP. SNR remains
disabled in the ordinary pipeline. Earlier isolated measurements were about
8.8× faster for movie rendering plus PNG roundtrip, and about 497× for the RMS
kernel; neither is included in the full-observation timing above.

**Full-observation numerical comparison.** Both versions used limb observation
0, 2000 × 1500 detector pixels, 40 wavelengths, nine 0.035-second integrations,
four 15-second integrations, both saved PSFs, and seed `20260904`. The baseline
ran directly from an archived copy of the pushed commit. Outputs and full
intermediate arrays were written to temporary directories; no production
science products were overwritten. Array snapshots added roughly five seconds
per run, measured separately and subtracted from the compute times above.

- Compared all 520 exposed wavelength maps, all 13 noisy detector maps, all 13
  pure electron maps when retained, and the final 750 × 1000 uint16 image.
- All noisy detector maps and the final image were bitwise identical in both the
  retained-reference and default optimized runs.
- All recorded map metadata matched, including per-member timestamps. Random
  generator states matched after every recorded pipeline stage.
- The largest per-map relative L2 optical difference was `1.71e-15`. The largest
  absolute optical difference was `7.46e-8` in the expected-photon map units.
- Pure electron maps differed by at most `3.35e-7` electron-equivalent counts per
  pixel; the largest per-map relative L2 difference was `4.54e-13`.
- FITS data and scientific header values matched; CHECKSUM values differ because
  Astropy includes the actual write time in checksum-related card comments.
  Every output passed CHECKSUM and DATASUM verification.
- The default comparison intentionally omits the 13 pure-reference maps. They
  are the only missing entries in its evidence manifest.

**PSF comparison.** Loaded the old PSF module directly from `482c35a`, with source
hash recorded in the evidence. The complete 2000 × 2000, 8,000-order, 185 Å,
20-lines-per-inch, cardinal-arm composite matched bitwise: zero differing pixels.
Both normalized sums were `1.0000000000002116`. Encircled energy at radii 1–1400
pixels and wing energy beyond 10–1000 pixels matched exactly. Convolution with a
real 1024² temperature-summed EM scene also matched bitwise. Additional bounded
tests cover other wavelengths, oblique arms, support-edge/distant orders,
representable subnormal tails, and exceptional inputs. This validates one full
production-sized wavelength case, not every possible user configuration.

**Regression checks.** The complete suite passed: **192 tests, no warnings**.
The real-data instrument integration test ran against a temporary data root
containing read-only links to calibrations and one EM input, with generated
radiance files confined to that temporary root. The CI blocking flake8 check
(`E9,F63,F7,F82`) and `git diff --check` passed. The CI advisory style check was
also run; the repository still has nonblocking style warnings.

Dedicated tests cover cache reuse and metadata independence, fractional starts,
terminal padding, single-integration stacks, independent Fano draws, unit
conversion, finite/nonfinite FFT paths, PSF storage ownership, Simpson grids,
digital overflow, radiance dtype/rounding, PSF support, movie pixels/transparency,
and optional SNR edge/dtype behavior. Independent code reviews checked the cache
and the radiance/PSF implementations.

**Limits.** The end-to-end comparison covers one representative full observation
with one fixed seed, not a complete CME sequence or every seed. Small floating
changes could cross a stochastic or rounding boundary in other runs, so
bit-identical output is not guaranteed universally. The evidence supports
preserved numerical accuracy and noise behavior in the tested cases; it is not
a new validation of the underlying physical model. Existing detector binning,
SNR reference scaling/window conventions, radiance bin endpoints, and movie
difference arithmetic were kept separate from this performance work.

Raw measurements are in [performance_evidence](performance_evidence/):
[default comparison](performance_evidence/default_comparison.json),
[retained-reference comparison](performance_evidence/retained_reference_comparison.json),
[radiance](performance_evidence/radiance.json),
[PSF](performance_evidence/psf.json), and
[runtime source hashes](performance_evidence/runtime_source_hashes.json).

Reproducible full-array runners, the baseline/optimized source snapshots, and
large comparison arrays were kept under `/tmp/suncet-validation` rather than in
the repository. The original assessment and bounded benchmark scripts are also
available in the earlier performance-assessment artifact.
