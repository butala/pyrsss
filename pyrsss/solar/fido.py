"""sunpy FIDO as the engine, for the missions it already serves.

Everything sunpy's Federated Internet Data Objects client can find should be
found through it -- it already speaks VSO and JSOC for SOHO/LASCO-EIT,
STEREO/SECCHI-EUVI, SDO/AIA, PSP/WISPR and Hinode/XRT, and re-implementing
search semantics would be the one unforgivable thing a fetch layer could do.
This module is the standardisation around it: one call shape, one manifest,
one cache directory.

The mission shims in this package exist only for the archives FIDO does not
reach -- IDOC pB and the SOHO orbits (``soho``), MLSO (``mlso``), PUNCH
(``punch``).

Import of ``sunpy`` is deferred so ``pyrsss.solar`` works without it.
"""

import logging
from pathlib import Path

logger = logging.getLogger('pyrsss.solar.fido')

# The instruments FIDO covers that tomography cares about, as
# (fido instrument, what it is).
MISSIONS = {
    'lasco': ('LASCO', 'SOHO/LASCO C2,C3 -- total brightness (pB: see soho)'),
    'eit': ('EIT', 'SOHO/EIT EUV'),
    'secchi': ('SECCHI', 'STEREO/COR1,COR2,EUVI -- L2 available'),
    'aia': ('AIA', 'SDO/AIA -- L1.5 science-ready'),
    'wispr': ('WISPR', 'Parker Solar Probe white light -- L2'),
    'xrt': ('XRT', 'Hinode/XRT soft X-ray'),
}


def _fido():
    try:
        from sunpy.net import fido
    except ImportError as e:
        raise ImportError(
            'FIDO support needs sunpy: uv sync --extra solar') from e
    return fido


def search(instrument, start, end, **attrs):
    """
    Query FIDO for *instrument* (a ``MISSIONS`` key or a sunpy instrument
    string) between *start* and *end*. Returns the sunpy query result;
    callers wanting just files can iterate its ``file`` column.
    """
    fido = _fido()
    name = MISSIONS.get(instrument, (instrument,))[0]
    logger.info('FIDO search %s in [%s, %s]', name, start, end)
    return fido.search(fido.attrs.Time(start, end),
                       fido.attrs.Instrument(name),
                       **attrs)


def fetch(query, out_dir=Path('.')):
    """
    Download a FIDO *query* into *out_dir*. Returns sunpy's fetch result,
    whose paths are what the manifest should cover.
    """
    fido = _fido()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    logger.info('fetching %d file(s) into %s', len(query), out_dir)
    return fido.fetch(query, path=str(out_dir / '{file}'))
