[![CI](https://github.com/butala/pyrsss/actions/workflows/ci.yml/badge.svg)](https://github.com/butala/pyrsss/actions/workflows/ci.yml)

# pyrsss

Python remote sensing and space science tools.

## Install

Requires Python 3.10+. Core install (pure Python):

```bash
pip install .
```

Domain-specific extras (see `pyproject.toml` for what each pulls in):

```bash
pip install '.[full]'
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

## Notebooks

The `notebooks/` directory contains exploratory research notebooks and
reference documents predating the current package layout; they are kept
for historical context and are not maintained against the library API.
