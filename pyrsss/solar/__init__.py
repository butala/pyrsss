"""Solar data acquisition for pyrsss: the archive layer tomography needs.

The division of labour is the whole design:

* **sunpy FIDO is the engine** for every archive it serves (SOHO/LASCO-EIT,
  STEREO/SECCHI-EUVI, SDO/AIA, PSP/WISPR, Hinode/XRT) -- see :mod:`fido`.
  Re-implementing search semantics is the one unforgivable thing a fetch
  layer could do.
* **Mission shims** cover what FIDO cannot reach, each small and honest
  about what it knows: :mod:`soho` (IDOC LASCO C2 *polarized* brightness
  and the SOHO orbit .DAT files -- ported from the ZJUI/zheyuan fetcher),
  :mod:`mlso` (K-Cor L2 ``pbavg``, the only ground-based pB), :mod:`punch`
  (full-Sun polarized mosaics, via the mission's own ``punchbowl``).
* **Provenance** is a first-class output: :mod:`manifest` writes and
  checks ``SHA256SUMS`` beside every cache, so a campaign's data can be
  re-verified and a re-fetch told from a silently different product.

Tomography's quantity of interest is polarized brightness (Thomson
scattering), which is why the shims lead with pB products: IDOC kfcorona,
K-Cor pbavg, PUNCH L3_PAM. The calibrated FIDO missions cover total
brightness and EUV.

On the shelf for the next wave (archives known, shims not yet written):
ASO-S/LST (China, Ly-alpha + white-light coronagraph,
``aso-s.pmo.ac.cn/sodc``), Aditya-L1 VELC/SUIT (India,
``pradan.issdc.gov.in/al1``), PROBA-3/ASPIICS (ESA inner corona) and
Solar Orbiter/Metis (SOAR). Install the extras group to pull the
dependencies this family wants::

    pip install 'pyrsss[solar]'

PUNCH support additionally wants ``pip install punchbowl`` (its dependency
footprint is heavier than the rest and stays its own choice).
"""

from .manifest import sha256, write_manifest, read_manifest, verify
from . import fido, mlso, punch, soho

__all__ = ['sha256', 'write_manifest', 'read_manifest', 'verify',
           'fido', 'mlso', 'punch', 'soho']
