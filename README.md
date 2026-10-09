[![CI](https://github.com/butala/pyrsss/actions/workflows/ci.yml/badge.svg)](https://github.com/butala/pyrsss/actions/workflows/ci.yml)

# pyrsss

Python remote sensing and space science tools.

## Install

Requires Python 3.10+. Core install (pure Python, with [uv](https://docs.astral.sh/uv/)):

```bash
uv sync
```

Domain-specific extras (see `pyproject.toml` for what each pulls in):

```bash
uv sync --extra full          # or --extra solar, --extra mag, ...
```

### External tools (RINEX front end)

The RINEX front end (`pyrsss-gnss-preprocess`, `pyrsss-gnss-phase-edit`,
`pyrsss-gnss-rinex`, `pyrsss-gnss-process`) shells out to external
executables:

- [teqc](https://www.unavco.org/software/data-processing/teqc) (RINEX
  normalization and summaries)
- GNSSTk/GPSTk `DiscFix` (cycle-slip detection) and `RinDump`
  (observable dump) --- set the `GNSSTK_BUILD` environment variable to
  the GNSSTk build directory containing these tools

Everything else (leveling, bias estimation and calibration, IONEX
handling, ...) is pure Python and needs no external tools.

RINEX can also be read without the GNSSTk tools via
`pyrsss.gnss.rinex.RinexDump.from_rinex(obs_fname, nav_fname)`
(georinex + hatanaka; RINEX 2/3 including Hatanaka-compressed input,
with optional broadcast-orbit geometry for the az/el and satellite ECEF
columns).

## Solar data (`pyrsss.solar`)

Archive acquisition for solar tomography, on the division of labour that
keeps it small: **sunpy FIDO is the engine** for the missions it serves
(SOHO, STEREO, SDO, PSP/WISPR, Hinode), and thin mission shims cover what
FIDO cannot reach -- `soho` (IDOC LASCO C2 *polarized* brightness and the
SOHO orbit `.DAT` files), `mlso` (K-Cor L2 `pbavg`), `punch` (full-Sun
polarized mosaics, via the mission's `punchbowl`). Every cache gets a
`SHA256SUMS` manifest (`pyrsss.solar.manifest`) so a campaign's data can be
re-verified.

```bash
uv sync --extra solar
pyrsss-solar-fetch-orbit 2008-02-01 --out data/orbit
pyrsss-solar-fetch-kcor 2019-02-28 --out data/kcor
```

The geometry half of tomography lives beside it: `registry` is the
catalog of every calibrated instrument with a knowable position (FOV
model, archive, product level), `fov` draws them as rectangular pyramids
and coronagraph annuli in direction space, `ephemeris` places the
observer (SOHO orbit files, ground sites, or JPL HORIZONS), `overlap`
finds the shared sky and `overlap.compare` returns the intercalibration
constant of two images over it. `viz4d` puts all of it in one 4-D scene
with pyviz4d -- source locations and FOV pyramids animated over time
(`uv sync --extra solar-viz`, then
`pyrsss-solar-fov-scene --ids lasco_c2,kcor --live`).

On the shelf: ASO-S/LST, Aditya-L1 VELC/SUIT, PROBA-3/ASPIICS,
Solar Orbiter/Metis.

## Notebooks

The `notebooks/` directory contains exploratory research notebooks and
reference documents predating the current package layout; they are kept
for historical context and are not maintained against the library API.
