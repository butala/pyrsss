"""
Tests for pyrsss.solar: the fetch layer's parse/manifest logic, with every
network call mocked out. The suite must run offline -- these modules are
only trustworthy if their *handling* is pinned independent of the archives
being up.
"""

from datetime import datetime
from types import SimpleNamespace
from pathlib import Path

import pytest

from pyrsss.solar import manifest, mlso, soho

FIXTURES = Path(__file__).parent / 'fixtures'


# ------------------------------------------------------------------ manifest

def test_manifest_round_trip(tmp_path):
    for name in ('a.fits', 'b.fits'):
        (tmp_path / name).write_bytes(name.encode())
    m = manifest.write_manifest(tmp_path)
    ok, missing, changed = manifest.verify(tmp_path)
    assert ok == ['a.fits', 'b.fits']
    assert missing == [] and changed == []
    assert manifest.read_manifest(m)['a.fits'] == manifest.sha256(
        tmp_path / 'a.fits')


def test_manifest_reports_missing_and_changed(tmp_path):
    (tmp_path / 'a.fits').write_bytes(b'a')
    (tmp_path / 'b.fits').write_bytes(b'b')
    manifest.write_manifest(tmp_path)
    (tmp_path / 'b.fits').write_bytes(b'tampered')
    (tmp_path / 'a.fits').unlink()
    ok, missing, changed = manifest.verify(tmp_path)
    assert ok == [] and missing == ['a.fits'] and changed == ['b.fits']


def test_manifest_covers_one_subdirectory_tree(tmp_path):
    (tmp_path / 'x.bin').write_bytes(b'x')
    manifest.write_manifest(tmp_path)
    names = manifest.read_manifest(tmp_path / manifest.MANIFEST_NAME)
    assert list(names) == ['x.bin']


# ---------------------------------------------------------------------- soho

def test_orbit_record_parses_a_dat_line():
    line = (FIXTURES / 'soho_orbit_sample.DAT').read_text().splitlines()[0]
    rec = soho.OrbitRecord.parse_line(line)
    assert rec.date == '01-Feb-2008'
    assert rec.year == 2008 and rec.doy == 32 and rec.ellapsed_ms == 0
    assert rec.cr_soho == 20
    assert rec.hec_x_km == pytest.approx(-1.33e8)
    assert rec.obstime == datetime(2008, 2, 1, 0, 0, 0)


def test_orbit_record_rejects_a_short_line():
    with pytest.raises(ValueError, match=r'expected \d+ fields'):
        soho.OrbitRecord.parse_line('01-Feb-2008 00:00:00.000 2008 32 0')


def test_parse_orbit_file_and_closest_record():
    records = soho.parse_orbit_file(FIXTURES / 'soho_orbit_sample.DAT')
    assert len(records) == 2
    rec, offset = soho.closest_orbit(
        records, datetime(2008, 2, 1, 0, 0, 30))
    assert rec.ellapsed_ms == 15000
    assert abs(offset.total_seconds()) == 15


def test_find_orbit_dat_takes_the_newest_revision(tmp_path):
    for name in ('SOHO_20080201.DAT', 'SOHO_20080201_02.DAT',
                 'SOHO_20080202.DAT'):
        (tmp_path / name).write_text('')
    assert soho.find_orbit_dat(
        tmp_path, datetime(2008, 2, 1)).name == 'SOHO_20080201_02.DAT'
    with pytest.raises(ValueError, match='not found'):
        soho.find_orbit_dat(tmp_path, datetime(2008, 3, 1))


def test_urls_are_the_documented_templates():
    dt = datetime(2008, 2, 1)
    assert soho.pb_url(dt, 'foo.fts') == (
        'http://idoc-lasco.ias.u-psud.fr/sitools/datastorage/user/results/'
        'kfcorona_sph_new_optimized/C2/Orange/2008/foo.fts')
    assert soho.orbit_url(dt) == (
        'https://soho.nascom.nasa.gov/data/ancillary/orbit/predictive/2008')


def test_fetch_url_caches_and_does_not_refetch(tmp_path, monkeypatch):
    calls = []

    class FakeResponse:
        status_code = 200
        content = b'payload'

    def fake_get(url, timeout=None):
        calls.append(url)
        return FakeResponse()

    monkeypatch.setattr(soho, 'requests', SimpleNamespace(get=fake_get))
    dest = tmp_path / 'sub' / 'file.bin'
    assert soho.fetch_url('http://example/file.bin', dest) == dest
    assert dest.read_bytes() == b'payload'
    assert soho.fetch_url('http://example/file.bin', dest) == dest
    assert len(calls) == 1, 'a cached file must not be re-fetched'


def test_fetch_url_raises_on_http_error(tmp_path, monkeypatch):
    class FakeResponse:
        status_code = 404
        content = b''

    monkeypatch.setattr(
        soho, 'requests',
        SimpleNamespace(get=lambda url, timeout=None: FakeResponse()))
    with pytest.raises(RuntimeError, match='404'):
        soho.fetch_url('http://example/gone', tmp_path / 'gone')


# ---------------------------------------------------------------------- mlso

def test_mlso_query_is_the_pbavg_product():
    url = mlso.query('kcor', datetime(2019, 2, 28))
    assert url.startswith(mlso.MLSO_GET)
    assert 'date1=2019-02-28' in url
    assert 'inst=kcor' in url and 'level=l2' in url
    assert 'proc=pbavg' in url and 'qual=all' in url


def test_mlso_listing_extracts_fits_names():
    html = (FIXTURES / 'mlso_listing.html').read_text()
    assert mlso.parse_listing(html) == [
        '20190228_174307_kcor_l2_pbavg.fts.gz',
        '20190228_174807_kcor_l2_pbavg.fts.gz']


def test_mlso_listing_ignores_non_products():
    assert mlso.parse_listing('<a href="index.html">x</a>') == []
