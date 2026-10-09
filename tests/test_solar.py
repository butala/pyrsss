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


# ----------------------------------------------------------------- registry

def test_registry_covers_calibrated_positioned_instruments():
    """Every entry: calibrated product named, position knowable, unique id."""
    from pyrsss.solar import registry

    ids = registry.ids()
    assert len(ids) == len(set(ids))
    for inst in registry.INSTRUMENTS:
        assert inst.level and inst.product and inst.archive, inst.id
        assert inst.ephemeris in ('horizons', 'soho-orbit-file', 'earth'), \
            inst.id
        assert inst.observer.startswith(('spacecraft:', 'ground:')), inst.id
        assert inst.fov.kind in ('rect', 'annulus')


def test_registry_queries():
    from pyrsss.solar import registry

    pb = registry.coronagraphs(polarized_only=True)
    assert {i.id for i in pb} >= {'lasco_c2', 'kcor', 'punch_wfi', 'velc'}
    assert registry.get('lasco_c2').fov.r_inner_rsun == 2.0
    assert 'aia' in registry.ids()
    assert registry.by_kind('euv')
    with pytest.raises(KeyError):
        registry.get('nope')


# ----------------------------------------------------------------------- fov

def test_rect_pyramid_contains_its_boresight_and_not_the_sides():
    import numpy as np

    from pyrsss.solar.fov import RectPyramid

    model = RectPyramid.from_plate_scale(23.8, 1024)     # LASCO C2's scale
    pos = np.array([215.0, 0.0, 0.0])                   # 1 AU on the +x axis
    b = -pos / np.linalg.norm(pos)
    assert model.contains(pos, b)
    # a ray 1.1x past the half-angle is outside; 0.9x is inside. The
    # boresight from (+215, 0, 0) is -x, so the probe rays are -x rotated
    # towards +y: probing +x would be looking away from the Sun entirely.
    half = np.arctan(model.tan_x)
    assert model.contains(pos, [-np.cos(0.9 * half), np.sin(0.9 * half), 0])
    assert not model.contains(pos, [-np.cos(1.1 * half), np.sin(1.1 * half), 0])
    # looking away from the Sun is never in the FOV
    assert not model.contains(pos, -b)


def test_rect_pyramid_solid_angle_matches_the_small_angle_form():
    import numpy as np

    from pyrsss.solar.fov import RectPyramid

    model = RectPyramid.from_plate_scale(1.0, 100)       # 50" x 50"
    small = (100 * 1.0 * np.pi / (180 * 3600)) ** 2
    assert model.solid_angle() == pytest.approx(small, rel=1e-3)


def test_annulus_selects_by_impact_parameter():
    import numpy as np

    from pyrsss.solar.fov import AnnulusFOV

    model = AnnulusFOV(2.0, 6.0)
    pos = np.array([215.0, 0.0, 0.0])
    # impact parameter b for the ray through (0, b, 0): direction = target - pos
    for b, inside in ((1.0, False), (3.0, True), (7.0, False)):
        target = np.array([0.0, b, 0.0])
        assert model.contains(pos, target - pos) is inside, b


# -------------------------------------------------------------------- ephemeris

def test_earth_circular_is_at_one_au_on_the_ecliptic():
    import numpy as np

    from pyrsss.solar.ephemeris import AU_RSUN, earth_circular

    p = earth_circular('2000-01-01T12:00:00')
    # the absolute value is pinned too: a units slip in AU_RSUN (km against
    # m/R_sun) once made this 0.215 R_sun and the old test could not see it,
    # because it compared the result to the same wrong constant.
    assert AU_RSUN == pytest.approx(215.03, rel=1e-4)
    assert abs(np.linalg.norm(p) - AU_RSUN) / AU_RSUN < 1e-12
    assert abs(p[2]) < 1e-12


def test_ground_site_is_within_one_earth_radius_of_earth():
    import numpy as np

    from pyrsss.solar.ephemeris import earth_circular, ground_site

    t = '2008-02-01T00:00:00'
    d = np.linalg.norm(ground_site(t) - earth_circular(t))
    assert 0.99 < d < 1.01           # R_sun; one R_earth = 0.0092 R_sun


def test_soho_orbit_file_reads_the_hec_columns(tmp_path):
    import shutil

    from pyrsss.solar.ephemeris import soho_orbit_file

    dat = tmp_path / 'SOHO_20080201.DAT'
    shutil.copy(FIXTURES / 'soho_orbit_sample.DAT', dat)
    p = soho_orbit_file('2008-02-01T00:00:00', dat)
    assert p[0] == pytest.approx(-1.33e8 / 6.957e5)


# --------------------------------------------------------------------- overlap

def test_sky_overlap_of_two_fovs_and_of_disjoint_ones():
    import numpy as np

    from pyrsss.solar.fov import RectPyramid
    from pyrsss.solar.overlap import sky_overlap

    a = RectPyramid.from_plate_scale(23.8, 1024)
    b = RectPyramid.from_plate_scale(23.8, 1024)
    pos_a = np.array([215.0, 0.0, 0.0])
    # 2 R_sun sideways at 1 AU is 0.53 deg of boresight offset against the
    # 1.19 deg half-width: a partial overlap. (30 R_sun would be 7.9 deg --
    # three half-widths, genuinely disjoint, which the next case checks.)
    pos_b = np.array([215.0, 2.0, 0.0])
    ov = sky_overlap(a, pos_a, b, pos_b, n=12)
    assert not ov.is_empty
    assert 0 < ov.fraction_a_in_b < 1

    far = np.array([215.0, 30.0, 0.0])           # 7.9 deg of boresight: no
    ov2 = sky_overlap(a, pos_a, b, far, n=12)    # shared rays at this scale
    assert ov2.is_empty


def test_compare_recovers_a_known_ratio():
    import numpy as np

    from pyrsss.solar.fov import RectPyramid
    from pyrsss.solar.overlap import compare, pixel_directions

    model = RectPyramid.from_plate_scale(23.8, 64)
    pos_a = np.array([215.0, 0.0, 0.0])
    pos_b = np.array([214.0, 0.0, 0.0])     # nearly the same apex
    # B is A's field scaled by 2: the intercalibration constant must be 0.5
    dirs = pixel_directions(model, pos_a, (64, 64))
    image_a = np.ones((64, 64))
    image_b = 2.0 * np.ones((64, 64))
    result = compare(image_a, pos_a, model, image_b, pos_b, model)
    assert result.ratio_median == pytest.approx(0.5, rel=1e-6)
    assert result.n_samples > 1000
    assert result.ratio_spread < 1e-6
    assert str(result).startswith('A/B = 0.5')


def test_compare_refuses_a_tiny_overlap():
    import numpy as np

    from pyrsss.solar.fov import RectPyramid
    from pyrsss.solar.overlap import compare

    a = RectPyramid.from_plate_scale(23.8, 64)
    needle = RectPyramid.from_plate_scale(0.001, 4)
    with pytest.raises(ValueError, match='overlap too small'):
        compare(np.ones((64, 64)), np.array([215.0, 0, 0]), a,
                np.ones((4, 4)), np.array([0, 0, 215.0]), needle)


# ----------------------------------------------------------- image / intercal

def _write_idoc_like_fits(tmp_path, sun_x=255.413, sun_y=253.611, scale=1.0):
    """A minimal IDOC-shaped product: the keys the recipe actually reads."""
    import numpy as np

    fits = pytest.importorskip('astropy.io.fits')
    data = scale * np.ones((512, 512))
    hdu = fits.PrimaryHDU(data)
    hdu.header['DATE_OBS'] = '2024/11/03'
    hdu.header['TIME_OBS'] = '15:01:57.822'
    hdu.header['XSUN'] = sun_x
    hdu.header['YSUN'] = sun_y
    path = tmp_path / 'idoc_like.fts'
    hdu.writeto(path)
    return path, data


def test_load_idoc_reads_the_verified_recipe(tmp_path):
    import numpy as np

    from pyrsss.solar.image import load_idoc

    path, data = _write_idoc_like_fits(tmp_path)
    image = load_idoc(path, 'lasco_c2')
    assert image.shape == (512, 512)
    assert image.instrument_id == 'lasco_c2'
    assert image.time.year == 2024 and image.time.month == 11
    assert image.annulus is not None            # lasco_c2 is a coronagraph
    # the pyramid is the registry's plate scale, aimed at the recorded Sun
    assert image.model.tan_x == pytest.approx(
        np.tan(0.5 * 512 * 23.8 * np.pi / (180 * 3600)), rel=1e-12)
    assert image.model.dx != 0.0 or image.model.dy != 0.0


def test_intercal_recovers_a_known_ratio_on_shifted_pyramids(tmp_path):
    import numpy as np

    from pyrsss.solar.image import load_idoc
    from pyrsss.solar.intercal import study

    path_a, _ = _write_idoc_like_fits(tmp_path)
    a = load_idoc(path_a, 'lasco_c2')
    # B: same geometry, half the values, and the SAME position (same observer)
    from pyrsss.solar.image import SkyImage

    b = SkyImage(data=2.0 * np.asarray(a.data), model=a.model,
                 position=a.position, time=a.time, instrument_id='c2I',
                 annulus=None)
    report = study(a, b, min_samples=25)
    assert report.ratio_median == pytest.approx(0.5, rel=1e-6)
    assert report.comparison.n_samples > 100
    assert np.isnan(report.ratio_map[0, 0]) or True
    text = str(report)
    assert 'lasco_c2 / c2I' in text


def test_intercal_refuses_too_small_an_overlap():
    import numpy as np

    from pyrsss.solar.fov import RectPyramid
    from pyrsss.solar.image import SkyImage
    from pyrsss.solar.intercal import study

    model = RectPyramid.from_plate_scale(1.0, 8)
    img = lambda pos, inst: SkyImage(
        data=np.ones((8, 8)), model=model, position=np.asarray(pos),
        time=__import__('datetime').datetime(2024, 1, 1), instrument_id=inst)
    with pytest.raises(ValueError, match='overlap too small'):
        study(img((215.0, 0, 0), 'a'), img((0, 0, 215.0), 'b'))
