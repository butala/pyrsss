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

### Optional `pyrsss.gnsstk` extension

A few modules (`pyrsss.gnss.bias`, `pyrsss.gnss.ipp`, `pyrsss.igrf.line`,
`pyrsss.iri.iri_stec`, `pyrsss.gnss.iri_stec`, `pyrsss.util.los_integrator`)
need the Cython `gnsstk` extension built against a local GNSSTk/GPSTk tree.
Build with the environment variables set:

```bash
export GNSSTK_SRC=/path/to/gnsstk/src   # source tree
export GNSSTK_BUILD=/path/to/gnsstk/build  # built libgnsstk
pip install .
```

Without these variables a pure-Python build is performed and the modules
above are unavailable.
