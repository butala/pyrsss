"""The instrument registry: every calibrated source with a knowable position.

Two requirements put an instrument in this catalog:

1. **Calibrated data exist** -- a documented level-1.5/2-equivalent product in
   physical units, served by an archive that will still be there next year.
2. **The source position is knowable** -- observer ephemeris and pointing
   sufficient to place the instrument's field of view on the sky: a
   spacecraft state (orbit file, HORIZONS, SPICE) or a ground site, plus a
   plate scale / angular extent.

That is the tomography contract: a measurement is only a measurement if you
know where it looked from. ``FOV`` below is the geometric half of it (see
:mod:`pyrsss.solar.fov` for the models); ``Ephemeris`` names where the
position comes from (see :mod:`pyrsss.solar.ephemeris`).

Values are *nominal* instrument design numbers -- the ones a survey draws
with -- and every entry says its source archive and product so a fetcher
(:mod:`pyrsss.solar.fido`, :mod:`pyrsss.solar.soho`, ...) can be written
against it. Entries marked ``status='historical'`` fly no longer but their
calibrated archives are open (long-baseline tomography wants them).
"""

from dataclasses import dataclass, field


@dataclass(frozen=True)
class FOV:
    """An instrument's field of view, as geometry.

    ``kind='rect'`` is a rectangular pyramid: a tangent-plane imager with
    ``nx`` x ``ny`` pixels of ``plate_scale`` arcsec each, looking along its
    boresight (towards the Sun's centre for every instrument in the table).

    ``kind='annulus'`` is a coronagraph's circular annulus on the sky:
    radii in R_sun as seen from the observer, with an optional azimuth
    sector (``sector_deg`` = None means full 360). A coronagraph *is* a
    rectangular pyramid of rays with the annulus mask laid over it; this
    entry records the mask, and the plate scale records the pixels.
    """

    kind: str                     # 'rect' | 'annulus'
    plate_scale_arcsec: float = None
    nx: int = None
    ny: int = None
    r_inner_rsun: float = None    # annulus only
    r_outer_rsun: float = None    # annulus only
    sector_deg: float = None      # annulus only; None = full circle

    def __post_init__(self):
        if self.kind not in ('rect', 'annulus'):
            raise ValueError(f'unknown FOV kind {self.kind!r}')
        if self.kind == 'rect' and (self.plate_scale_arcsec is None
                                    or self.nx is None):
            raise ValueError('rect FOV needs plate_scale_arcsec and nx')
        if self.kind == 'annulus' and (self.r_inner_rsun is None
                                       or self.r_outer_rsun is None):
            raise ValueError('annulus FOV needs r_inner_rsun and r_outer_rsun')


@dataclass(frozen=True)
class Instrument:
    """One entry of the registry."""

    id: str                       # stable short name, e.g. 'lasco_c2'
    mission: str
    instrument: str
    kind: str                     # 'coronagraph-pb' | 'coronagraph' | 'euv' | ...
    archive: str                  # where the calibrated product lives
    level: str                    # the calibrated product's level name
    product: str                  # the product to fetch for tomography
    observer: str                 # spacecraft name or 'ground:MLSO' etc.
    ephemeris: str                # 'horizons' | 'soho-orbit-file' | 'earth'
    fov: FOV
    active: str                   # years flown, or 'ongoing'
    status: str = 'operational'   # 'operational' | 'historical' | 'commissioning'
    notes: str = ''


def _rect(scale, nx, ny=None):
    return FOV('rect', plate_scale_arcsec=scale, nx=nx, ny=ny or nx)


def _annulus(r_in, r_out, scale, nx):
    return FOV('annulus', plate_scale_arcsec=scale, nx=nx, ny=nx,
               r_inner_rsun=r_in, r_outer_rsun=r_out)


# ---------------------------------------------------------------- the table
#
# Coronagraphs first (polarized brightness where it exists -- Thomson
# scattering tomography's actual measurement), then the EUV/SXR imagers
# that constrain the inner corona, then the historical archives.

INSTRUMENTS = (
    # ---- coronagraphs ----
    Instrument('lasco_c2', 'SOHO', 'LASCO C2', 'coronagraph-pb',
               'IDOC / CDAW / SOHO archive', 'L0.5-L1 (pB: L2 synoptic, IDOC)',
               'kfcorona pB (IDOC) or L2 pB synoptic', 'spacecraft:SOHO',
               'soho-orbit-file', _annulus(2.0, 6.0, 23.8, 1024),
               '1996-ongoing',
               notes='the historical pB workhorse; the IDOC products are the '
                     'calibrated pB used for tomography'),
    Instrument('lasco_c3', 'SOHO', 'LASCO C3', 'coronagraph-pb',
               'IDOC / CDAW / SOHO archive', 'L0.5-L1', 'L1 radiance',
               'spacecraft:SOHO', 'soho-orbit-file',
               _annulus(3.7, 30.0, 56.0, 1024), '1996-ongoing',
               notes='wider annulus, coarser plate scale'),
    Instrument('secchi_cor1', 'STEREO', 'COR1', 'coronagraph-pb',
               'STEREO Science Center', 'L2', 'L2 polarized', 'spacecraft:STEREO-A',
               'horizons', _annulus(1.4, 4.0, 75.0, 512),
               '2006-ongoing', notes='STEREO-A still flying; B lost 2014'),
    Instrument('secchi_cor2', 'STEREO', 'COR2', 'coronagraph-pb',
               'STEREO Science Center', 'L2', 'L2 polarized', 'spacecraft:STEREO-A',
               'horizons', _annulus(2.5, 15.0, 14.7, 512),
               '2006-ongoing'),
    Instrument('punch_wfi', 'PUNCH', 'WFI x4', 'coronagraph-pb',
               'punchbowl / data.nasa.gov', 'L2 (PTM), L3_PAM',
               'L3_PAM polarized mosaic', 'spacecraft:PUNCH',
               'earth', _annulus(3.5, 30.0, 288.0, 1024),
               '2025-ongoing', notes='four wide-field imagers, full-Sun pB '
                                     'mosaics over 3.5-30 R_sun; Earth-orbit '
                                     'constellation, modelled at Earth'),
    Instrument('punch_nfi', 'PUNCH', 'NFI', 'coronagraph-pb',
               'punchbowl / data.nasa.gov', 'L2 (PTM)', 'L2 rectified NFI',
               'spacecraft:PUNCH', 'earth', _annulus(5.0, 15.0, 94.0, 1024),
               '2025-ongoing', notes='narrow-field inner corona 5-15 R_sun'),
    Instrument('wispr', 'Parker Solar Probe', 'WISPR', 'coronagraph',
               'wispr.nrl.navy.mil / SDAC', 'L2', 'L2 FITS (MSB)',
               'spacecraft:PSP', 'horizons', _rect(380.0, 1024),
               '2018-ongoing', notes='deep-space viewpoints; combined inner/'
                                     'outer telescopes span ~13-108 R_sun'),
    Instrument('metis', 'Solar Orbiter', 'Metis', 'coronagraph-pb',
               'SOAR (esa)', 'L2', 'L2 VL+UV, polarimetric mode',
               'spacecraft:SolarOrbiter', 'horizons',
               _annulus(1.6, 9.0, 85.0, 1024), '2020-ongoing',
               notes='VL + Ly-alpha; polarimetric sequences give pB'),
    Instrument('aspiics', 'PROBA-3', 'ASPIICS', 'coronagraph',
               'ESA SOC', 'L2', 'L2 VL', 'spacecraft:PROBA-3', 'horizons',
               _annulus(1.05, 3.0, 4.5, 1024), '2024-ongoing',
               status='commissioning',
               notes='formation-flying occulter: the first 1.05 R_sun '
                     'inner edge from space'),
    Instrument('kcor', 'MLSO', 'K-Cor (COSMO)', 'coronagraph-pb',
               'mlso.hao.ucar.edu', 'L1/L2', 'L2 pbavg', 'ground:MLSO',
               'earth', _annulus(1.05, 3.0, 5.6, 1024), '2013-ongoing',
               notes='the only ground-based pB; cloud/quality gated'),
    Instrument('ucomp', 'MLSO', 'UCoMP', 'euv',
               'mlso.hao.ucar.edu', 'L1', 'L1 intensity', 'ground:MLSO',
               'earth', _rect(156.0, 1024), '2021-ongoing',
               notes='coronal magnetography (Fe XIII); tomography '
                     'constraints rather than density'),
    Instrument('velc', 'Aditya-L1', 'VELC', 'coronagraph-pb',
               'PRADAN/ISSDC', 'L2', 'L2 imaging 500 nm + polarimetry',
               'spacecraft:Aditya-L1', 'horizons',
               _annulus(1.05, 3.0, 5.0, 1024), '2023-ongoing',
               notes='internally occulted; imaging + spectro-polarimetry'),
    Instrument('lst', 'ASO-S', 'LST', 'coronagraph-pb',
               'aso-s.pmo.ac.cn SODC', 'L1/L2', 'L2 Ly-alpha + WL',
               'spacecraft:ASO-S', 'horizons',
               _annulus(1.05, 2.5, 6.0, 1024), '2022-ongoing',
               notes='Ly-alpha and white-light inner corona together; '
                     'perihelion quarters carry stray-light caveats'),
    # ---- EUV / SXR imagers (inner-corona constraints, limb views) ----
    Instrument('aia', 'SDO', 'AIA', 'euv',
               'JSOC', 'L1.5 (aia_prep)', 'L1.5 FITS',
               'spacecraft:SDO', 'horizons', _rect(60.0, 4096),
               '2010-ongoing', notes='9 EUV channels; science-ready is L1.5'),
    Instrument('secchi_euvi', 'STEREO', 'EUVI', 'euv',
               'STEREO Science Center', 'L2', 'L2 FITS', 'spacecraft:STEREO-A',
               'horizons', _rect(56.0, 2048), '2006-ongoing'),
    Instrument('eit', 'SOHO', 'EIT', 'euv',
               'SOHO archive', 'L1', 'L1 calibrated', 'spacecraft:SOHO',
               'soho-orbit-file', _rect(102.0, 1024), '1996-ongoing'),
    Instrument('swap', 'PROBA-2', 'SWAP', 'euv',
               'PROBA-2 archive (ROB)', 'L1/L2', 'L2 FITS',
               'spacecraft:PROBA-2', 'horizons', _rect(192.0, 1024),
               '2009-ongoing'),
    Instrument('suvi', 'GOES-16..18', 'SUVI', 'euv',
               'NOAA NCEI / NGDC', 'L2', 'L2 netCDF', 'spacecraft:GOES',
               'horizons', _rect(192.0, 1024), '2016-ongoing',
               notes='geostationary: fixed Earth view, 6 EUV channels'),
    Instrument('xrt', 'Hinode', 'XRT', 'sxr',
               'Hinode Science Centre', 'L1', 'L1 calibrated',
               'spacecraft:Hinode', 'horizons', _rect(512.0, 512),
               '2006-ongoing'),
    Instrument('suit', 'Aditya-L1', 'SUIT', 'euv',
               'PRADAN/ISSDC', 'L1/L2', 'L2 UV images', 'spacecraft:Aditya-L1',
               'horizons', _rect(40.0, 2048), '2023-ongoing',
               notes='200-400 nm UV disk imager beside VELC'),
    Instrument('sxt_asos', 'ASO-S', 'SXT', 'sxr',
               'aso-s.pmo.ac.cn SODC', 'L1/L2', 'L2 SXR',
               'spacecraft:ASO-S', 'horizons', _rect(100.0, 1024),
               '2022-ongoing'),
    # ---- historical (calibrated archives, no longer flying) ----
    Instrument('sxt_yohkoh', 'Yohkoh', 'SXT', 'sxr',
               'ISAS/NAOJ archive', 'L1', 'L1 calibrated',
               'spacecraft:Yohkoh', 'horizons', _rect(204.0, 1024),
               '1991-2001', status='historical',
               notes='solar-maximum SXR tomography baseline'),
    Instrument('tesis', 'CORONAS-F/Photon', 'TESIS EVT', 'euv',
               'IKI archive', 'L1', 'L1 EUV', 'spacecraft:CORONASF',
               'horizons', _rect(180.0, 1024), '2009-2010',
               status='historical'),
)


_BY_ID = {inst.id: inst for inst in INSTRUMENTS}


def get(instrument_id):
    """The ``Instrument`` for *instrument_id*; raises ``KeyError``."""
    return _BY_ID[instrument_id]


def by_mission(mission):
    """Every instrument of a mission (exact match on ``mission``)."""
    return tuple(i for i in INSTRUMENTS if i.mission == mission)


def by_kind(kind):
    """Every instrument of one kind ('coronagraph-pb', 'euv', ...)."""
    return tuple(i for i in INSTRUMENTS if i.kind == kind)


def coronagraphs(polarized_only=False):
    """The coronagraphs -- where tomography's pB measurements come from."""
    return tuple(i for i in INSTRUMENTS
                 if i.kind.startswith('coronagraph')
                 and (not polarized_only or i.kind == 'coronagraph-pb'))


def ids():
    """Every registry id, in table order."""
    return tuple(i.id for i in INSTRUMENTS)
