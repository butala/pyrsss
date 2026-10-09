"""MLSO (Mauna Loa) data: K-Cor polarized brightness and UCoMP.

Ground-based, the only pB source that is not a spacecraft -- and the only one
that watches the inner corona (1.05-3 R_sun) continuously in good weather.
Products are served by HAO's ``mlso_data_get.php`` front end (and, since
2025, a documented web-service API); the tomography-relevant product is
K-Cor **L2 ``pbavg``** -- the polarized-brightness average, already
photometrically calibrated and cloud-quality gated.

Quality is the catch with any ground-based series: the ``qual`` parameter
exists because clouds and sky brightness gate whole frames, so a fetcher
must record the quality it took, not just the files.

Verified 2026-10 against the MLSO data-access pages; the PHP endpoint is
the stable public interface and is what this module drives.
"""

import logging
from datetime import datetime
from pathlib import Path
from urllib.parse import urlencode

try:
    import requests
except ImportError:  # core install: the parsers still work, fetching does not
    requests = None

logger = logging.getLogger('pyrsss.solar.mlso')


def _requests():
    """``requests``, or the error that says how to get it."""
    if requests is None:
        raise ImportError("fetching needs requests: pip install 'pyrsss[solar]'")
    return requests

# The public data front end. One call lists every file for an instrument /
# date / level / product; the response is an HTML page of download links.
MLSO_GET = 'https://mlso.hao.ucar.edu/mlso_data_get.php'

# The products tomography wants. ``proc`` is the processing level's product
# name; K-Cor's L2 polarized-brightness average is ``pbavg``.
PRODUCTS = {
    'kcor': {'inst': 'kcor', 'level': 'l2', 'proc': 'pbavg'},
    'ucomp': {'inst': 'ucomp', 'level': 'l1', 'proc': 'intensity'},
}


def query(instrument, date, level=None, proc=None, qual='all'):
    """
    The ``mlso_data_get.php`` query string for one instrument and date.

    *instrument* is a ``PRODUCTS`` key ('kcor', 'ucomp') or a raw
    ``inst=`` value; *level*/*proc* default to the product's calibrated
    choice. Returns the fully qualified URL.
    """
    spec = PRODUCTS.get(instrument, {'inst': instrument})
    params = {
        'date1': date.strftime('%Y-%m-%d'),
        'inst': spec['inst'],
        'level': level or spec.get('level', 'l2'),
        'qual': qual,
        'proc': proc or spec.get('proc', 'pbavg'),
    }
    return MLSO_GET + '?' + urlencode(params)


def parse_listing(html):
    """
    The file names a ``mlso_data_get`` page links to, in page order.

    The page is a plain list of download links; anything ending in FITS or
    fits.gz is a product (the page also links its own header rows).
    """
    import re
    names = re.findall(r'href="([^"]+\.(?:fts|fits)(?:\.gz)?)"', html,
                       flags=re.IGNORECASE)
    return [Path(name).name for name in names]


def fetch_day(instrument, date, out_dir, qual='all', force=False):
    """
    Fetch every listed product for one instrument-day into *out_dir*.

    Returns the local paths; the listing page itself is cached beside them
    as ``mlso_<inst>_<date>.html`` so the manifest can cover what was taken.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    url = query(instrument, date, qual=qual)
    listing_path = out_dir / (f'mlso_{instrument}_{date:%Y%m%d}.html')
    if not listing_path.exists() or force:
        logger.info('fetching listing %s', url)
        r = _requests().get(url, timeout=60)
        r.raise_for_status()
        listing_path.write_text(r.text)
    html = listing_path.read_text()
    paths = []
    for name in parse_listing(html):
        dest = out_dir / name
        if dest.exists() and not force:
            logger.info('using cached %s', dest)
        else:
            file_url = url.split('?')[0].replace('mlso_data_get.php',
                                                 f'data/{name}')
            logger.info('fetching %s -> %s', file_url, dest)
            r = _requests().get(file_url, timeout=60)
            r.raise_for_status()
            dest.write_bytes(r.content)
        paths.append(dest)
    return listing_path, paths


def main(argv=None):
    """``pyrsss-solar-fetch-kcor``: one day of K-Cor L2 pbavg."""
    import argparse

    parser = argparse.ArgumentParser(
        description='Fetch MLSO K-Cor L2 polarized-brightness data.')
    parser.add_argument('date', help='YYYY-MM-DD')
    parser.add_argument('--instrument', default='kcor',
                        choices=sorted(PRODUCTS), help='(default: kcor)')
    parser.add_argument('--qual', default='all',
                        help="quality gate (default: all)")
    parser.add_argument('--out', type=Path, default=Path('.'),
                        help='cache directory (default: .)')
    args = parser.parse_args(argv)

    date = datetime.strptime(args.date, '%Y-%m-%d')
    listing, paths = fetch_day(args.instrument, date, args.out, qual=args.qual)
    print(f'{listing}: {len(paths)} products')
    for path in paths:
        print(f'  {path}')
    return 0


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    raise SystemExit(main())
