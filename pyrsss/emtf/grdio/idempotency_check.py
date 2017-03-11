import grd_io as gi
import numpy as np

TEST_INPUT  = 'TEST_INPUT.grd'
TEST_OUTPUT = 'TEST_OUTPUT.grd'

lon_grid, lat_grid, time_grid, DATA = gi.grd_read(TEST_INPUT)
gi.grd_write(TEST_OUTPUT, lon_grid, lat_grid, time_grid, DATA)
lon_grid2, lat_grid2, time_grid2, DATA2 = gi.grd_read(TEST_OUTPUT)

print 'longitude grid idempotency:', np.all(lon_grid2 == lon_grid)
print 'latitude grid idempotency:',  np.all(lat_grid2 == lat_grid)
print 'data idempotency:',           np.all(DATA      == DATA2)
