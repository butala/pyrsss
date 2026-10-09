"""SOHO data: LASCO C2 polarized brightness (IDOC) and the SOHO orbit files.

Two products, one mission, and both are outside what sunpy's FIDO serves --
which is why this module exists at all:

* **pB**: the IDOC ``kfcorona_sph`` LASCO C2 orange polarized-brightness
  products (IAS/U-PSud, ``idoc-lasco.ias.u-psud.fr``). These are the
  calibrated pB synoptic products -- the quantity Thomson-scattering
  tomography actually measures -- rather than the total-brightness L0.5/L1
  LASCO frames FIDO hands back.
* **Orbit**: the NRL predictive orbit ``.DAT`` files under
  ``soho.nascom.nasa.gov/data/ancillary/orbit/predictive/<year>/``, whose 36
  fields give the spacecraft state in GCI/GSE/GSM/HEC plus Carrington
  coordinates -- the ephemeris SphericalCT's ``get_builda_coordinates``
  consumes.

Ported from ``ZJUI/zheyuan/graduate_thesis_project/python/fetch_lasco_c2.py``
(2026-10), which carried the two pieces of knowledge written down here: the
orbit record's exact field list, and the WCS recipe that turns an IDOC pB
image into an angular map (``create_pB_map``'s constants, kept as
``PB_WCS`` below). The fetchers cache by name -- an existing local file is
never re-fetched -- and everything is plain ``requests`` so the tests can run
against fixtures with no network.
"""

import logging
from dataclasses import dataclass, fields
from datetime import datetime
from pathlib import Path

try:
    import requests
except ImportError:  # core install: the parsers still work, fetching does not
    requests = None

logger = logging.getLogger('pyrsss.solar.soho')


def _requests():
    """``requests``, or the error that says how to get it."""
    if requests is None:
        raise ImportError("fetching needs requests: pip install 'pyrsss[solar]'")
    return requests


# IDOC product tree (LASCO C2, orange, spherical pB). Verified 2026-10 by the
# zheyuan fetcher that used it daily; the directory is <year>/<fname>.
IDOC_HOST = 'http://idoc-lasco.ias.u-psud.fr'
PB_TEMPLATE = (IDOC_HOST + '/sitools/datastorage/user/results/'
               'kfcorona_sph_new_optimized/C2/Orange/{dt:%Y}/{pB_fname}')

# SOHO NRL predictive orbit listings; one .DAT per day, revised in place
# (the revision is the two digits before ".DAT", highest wins).
ORBIT_DIR_TEMPLATE = ('https://soho.nascom.nasa.gov/data/ancillary/'
                      'orbit/predictive/{dt:%Y}')

# The angular model for an IDOC pB frame. IDOC ships the image with XSUN/
# YSUN (the Sun's centre in pixels) but not a projection: these turn it into
# a tangent-plane map at LASCO C2's plate scale. RSUN's units are the one
# trap -- see the IDOC "Important Warning Concerning Units" note; both
# numbers are the ones the zheyuan code shipped with.
PB_WCS = dict(CRVAL1=1.00000, CRVAL2=1.00000,
              CDELT1=23.799999, CDELT2=23.799999,
              CTYPE1='HPLN-TAN', CTYPE2='HPLT-TAN',
              CUNIT1='arcsec', CUNIT2='arcsec',
              RSUN_KM=695990, RSUN_ARCSEC=959.63,
              SCALE=1e-10)


@dataclass(frozen=True)
class OrbitRecord:
    """One line of a SOHO NRL orbit ``.DAT`` file: 38 tokens, in order.

    The file is the BINTABLE its own header documents (36 numeric columns:
    the state in GCI/GSE/GSM/HEC, the Sun vector, and the Carrington pairs)
    with ``date`` and ``time`` spelled at the front, so a line splits into
    38 whitespace tokens. Positions are km, velocities km/s; ``hec`` is the
    Heliocentric Ecliptic frame -- the one with the fewest steps to
    HeliographicCarrington, which is why the pB map builder starts there.
    ``date``/``time`` are strings exactly as the file spells them
    (``15-Sep-2008  00:00:00.000``).
    """

    date: str
    time: str
    year: int
    doy: int
    ellapsed_ms: int
    gci_x_km: float
    gci_y_km: float
    gci_z_km: float
    gci_vx_kms: float
    gci_vy_kms: float
    gci_vz_kms: float
    gse_x_km: float
    gse_y_km: float
    gse_z_km: float
    gse_vx_kms: float
    gse_vy_kms: float
    gse_vz_kms: float
    gsm_x_km: float
    gsm_y_km: float
    gsm_z_km: float
    gsm_vx_kms: float
    gsm_vy_kms: float
    gsm_vz_kms: float
    gci_sun_x_km: float
    gci_sun_y_km: float
    gci_sun_z_km: float
    hec_x_km: float
    hec_y_km: float
    hec_z_km: float
    hec_vx_kms: float
    hec_vy_kms: float
    hec_vz_kms: float
    cr_earth: int
    heliographic_lon_earth: float
    heliographic_lat_earth: float
    cr_soho: int
    heliographic_lon_soho: float
    heliographic_lat_soho: float

    @classmethod
    def parse_line(cls, line):
        """
        Parse one whitespace-separated record. Field order is the class's,
        and the declared types are applied -- the file is text, so without
        the casts every field arrives as a string (the trap the original
        zheyuan fetcher papered over with a dataclass ``__post_init__``).
        """
        spec = fields(cls)
        values = line.split()
        if len(values) != len(spec):
            raise ValueError(
                f'expected {len(spec)} fields, got {len(values)}')
        return cls(*(f.type(v) for f, v in zip(spec, values)))

    @property
    def obstime(self):
        """The record's epoch as a ``datetime`` (the file's ``date time``)."""
        return datetime.strptime(self.date + ' ' + self.time + '000',
                                 '%d-%b-%Y %H:%M:%S.%f')


def fetch_url(url, dest, force=False):
    """
    Download *url* to *dest* unless it is already there (or *force*).
    Returns the local path; a non-200 is a ``RuntimeError``.
    """
    dest = Path(dest)
    if dest.exists() and not force:
        logger.info('using cached %s', dest)
        return dest
    logger.info('fetching %s -> %s', url, dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    r = _requests().get(url, timeout=60)
    if r.status_code != 200:
        raise RuntimeError(f'{url}: status code = {r.status_code}')
    with open(dest, 'wb') as fid:
        fid.write(r.content)
    return dest


def parse_orbit_file(path):
    """
    Parse a SOHO orbit ``.DAT`` into a list of ``OrbitRecord``, one per line.
    """
    records = []
    with open(path) as fid:
        for i, line in enumerate(fid, 1):
            line = line.strip()
            if not line:
                continue
            try:
                records.append(OrbitRecord.parse_line(line))
            except ValueError as e:
                raise ValueError(f'{path}:{i}: {e}')
    return records


def closest_orbit(records, target):
    """
    The record nearest *target* in time, and the offset (``record - target``).
    """
    if not records:
        raise ValueError('no orbit records')
    best = min(records,
               key=lambda x: (x.obstime > target,
                              abs((x.obstime - target).total_seconds())))
    return best, best.obstime - target


def find_orbit_dat(directory, date):
    """
    The newest ``*.DAT`` in *directory* covering *date* (the revision is the
    two digits before the extension; highest wins). ``directory`` may be a
    local directory or a ``requests``-fetchable listing saved by
    ``fetch_orbit_listing``. Raises ``ValueError`` if the day is absent.
    """
    directory = Path(directory)
    wanted = date.strftime('%Y%m%d')
    hits = []
    for path in directory.iterdir():
        name = path.name
        if wanted in name and name.upper().endswith('.DAT'):
            hits.append(name)
    if not hits:
        raise ValueError(f'{date:%Y-%m-%d} not found in {directory}')
    return directory / max(hits, key=lambda name: int(name[-6:-4]))


def orbit_url(date):
    """The predictive-orbit directory for *date*'s year."""
    return ORBIT_DIR_TEMPLATE.format(dt=date)


def pb_url(dt, pB_fname):
    """The IDOC pB product URL for the frame named *pB_fname* on *dt*."""
    return PB_TEMPLATE.format(dt=dt, pB_fname=pB_fname)


def main(argv=None):
    """``pyrsss-solar-fetch-orbit``: fetch one day's SOHO orbit file."""
    import argparse

    parser = argparse.ArgumentParser(
        description='Fetch one day of SOHO NRL predictive orbit data.')
    parser.add_argument('date', help='YYYY-MM-DD')
    parser.add_argument('--out', type=Path, default=Path('.'),
                        help='cache directory (default: .)')
    args = parser.parse_args(argv)

    date = datetime.strptime(args.date, '%Y-%m-%d')
    listing_name = f'NRL_orbits_{date:%Y}.txt'
    listing = fetch_url(orbit_url(date), args.out / listing_name)
    # The listing is an HTML directory page: the .DAT names are its hrefs.
    dat_name = None
    wanted = date.strftime('%Y%m%d')
    for line in open(listing):
        if wanted in line and '.DAT' in line and 'href="' in line:
            name = line.split('href="')[1].split('"')[0]
            dat_name = max(filter(None, [dat_name, name]),
                           key=lambda n: int(n[-6:-4]))
    if dat_name is None:
        raise SystemExit(f'no .DAT for {wanted} in {listing}')
    path = fetch_url(f'{orbit_url(date)}/{dat_name}', args.out / dat_name)
    records = parse_orbit_file(path)
    print(f'{path}: {len(records)} records, '
          f'{records[0].obstime} .. {records[-1].obstime}')
    return 0


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    raise SystemExit(main())
