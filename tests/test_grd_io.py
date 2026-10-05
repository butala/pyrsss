"""Round-trip tests for the .grd reader/writer (pyrsss.emtf.grdio)."""
import numpy as np

from pyrsss.emtf.grdio.grd_io import grd_read, grd_write


def test_roundtrip(tmp_path):
    n_lons, n_lats, n_times = 6, 5, 3
    lon_res, lon_west = 0.5, -120.0
    lat_res, lat_south = 0.25, 40.0
    lon_grid = lon_west + lon_res * np.arange(n_lons)
    lat_grid = lat_south + lat_res * np.arange(n_lats)
    time_grid = [100.0, 200.0, 300.0]
    rng = np.random.default_rng(42)
    DATA = rng.standard_normal((n_lons, n_lats, n_times))

    fname = str(tmp_path / 'test.grd')
    grd_write(fname, lon_grid, lat_grid, time_grid, DATA)
    lon_grid2, lat_grid2, times2, DATA2 = grd_read(fname)

    np.testing.assert_allclose(lon_grid2, lon_grid)
    np.testing.assert_allclose(lat_grid2, lat_grid)
    np.testing.assert_allclose(times2, time_grid)
    # .grd stores 11 significant digits
    np.testing.assert_allclose(DATA2, DATA, rtol=1e-9)


def test_roundtrip_17_lat_block(tmp_path):
    # the original reader hard-coded 17 latitude points per block
    n_lons, n_lats, n_times = 2, 17, 1
    lon_grid = -100.0 + 0.5 * np.arange(n_lons)
    lat_grid = 30.0 + 0.1 * np.arange(n_lats)
    time_grid = [1.0]
    rng = np.random.default_rng(7)
    DATA = rng.standard_normal((n_lons, n_lats, n_times))

    fname = str(tmp_path / 'test17.grd')
    grd_write(fname, lon_grid, lat_grid, time_grid, DATA)
    _, _, times2, DATA2 = grd_read(fname)

    np.testing.assert_allclose(times2, time_grid)
    np.testing.assert_allclose(DATA2, DATA, rtol=1e-9)
