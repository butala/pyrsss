import IPython
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import read_ionex

def plot_grid2D(lons, lats, tec_grid2D, datetime, title_label = ''):

	LATS, LONS = np.meshgrid(lats, lons)

	m = Basemap(llcrnrlon=-180,
	            llcrnrlat=-55,
	            urcrnrlon=180,
	            urcrnrlat=75,
	            projection='merc',
	            area_thresh=1000,
	            resolution='i')

	m.drawstates()
	m.drawcountries()
	m.drawcoastlines()

	parallels = np.arange(-90,90,20)
	m.drawparallels(parallels,labels=[True,False,False,True])
	meridians = np.arange(0,360,40)
	m.drawmeridians(meridians,labels=[True,False,False,True])

	m.scatter(LONS, LATS, c=tec_grid2D, latlon = True, linewidths=0, s=5)
	m.colorbar()

	plt.title('%s\n%s' % (title_label, datetime.isoformat(' ')))

FILE = 'igrg0010.16i'
print FILE

matplotlib.rcParams['savefig.dpi'] = 500
matplotlib.rcParams['figure.figsize'] = (5,4)
matplotlib.rcParams['font.size'] = 4

## RAW ##
print "Generating plot of raw data...",
raw_lons, raw_lats, time_grid, raw_tec_grids, raw_rms_grids, satellite_biases, station_biases = read_ionex.raw(FILE)
plt.subplot(211)
plot_grid2D(raw_lons, raw_lats, raw_tec_grids[:,:,0], time_grid[0], title_label = "Raw IONEX Data")
plt.subplot(212)
plot_grid2D(raw_lons, raw_lats, raw_rms_grids[:,:,0], time_grid[0])
plt.savefig('raw.png')
plt.close()
print "done."

# ## SPATIAL INTERPOLATION ##
print "Generating plot of spatially interpolated (nearest neighbor) data...",
dense_lons  = np.arange(-180, 180, 1)
dense_lats  = np.arange(-90, 90, 1)

time_grid, si_grid_tec, si_grid_rms, satellite_biases, station_biases = read_ionex.interpolate2D_spatial(FILE, spatial_grid = (dense_lons, dense_lats), method = 'nearest')
plt.subplot(211)
plot_grid2D(dense_lons, dense_lats, si_grid_tec[:,:,0], time_grid[0], title_label = "Spatial Interpolation (Nearest Neighbor); dLat, dLon = 1")
plt.subplot(212)
plot_grid2D(dense_lons, dense_lats, si_grid_rms[:,:,0], time_grid[0], title_label = "Spatial Interpolation (Nearest Neighbor); dLat, dLon = 1")
plt.savefig('nearest.png')
plt.close()
print "done."

print "Generating plot of spatially interpolated (bilinear) data...",
time_grid, si_grid_tec, si_grid_rms, satellite_biases, station_biases = read_ionex.interpolate2D_spatial(FILE, spatial_grid = (dense_lons, dense_lats), method = 'linear')
plt.subplot(211)
plot_grid2D(dense_lons, dense_lats, si_grid_tec[:,:,0], time_grid[0], title_label = "Spatial Interpolation (Bilinear); dLat, dLon = 1")
plt.subplot(212)
plot_grid2D(dense_lons, dense_lats, si_grid_rms[:,:,0], time_grid[0])
plt.savefig('bilinear.png')
plt.close()
print "done."

# ## TEMPORAL INTERPOLATION ##
print "Generating plot of temporally interpolated (linear) data...",
import datetime
 # 24 hours, 15 minute increments starting from the first epoch in the IONEX file
 # future versions will allow you to specify a 'dt' value instead of manually creating the temporal_grid
 # if this is a critically needed feature, let me know and I will prioritize it - Matt
temporal_grid = [time_grid[0] + datetime.timedelta(minutes=minute) for minute in range(0, 60*24+15, 15)]
ti_lons, ti_lats, ti_grid_tec, ti_grid_rms, satellite_biases, station_biases = read_ionex.interpolate2D_temporal(FILE, temporal_grid = temporal_grid, method = 'linear')
# plot the last entry
plot_grid2D(ti_lons, ti_lats, ti_grid_tec[:,:,-1], temporal_grid[-1], title_label = "Temporal Interpolation (Linear); dt = 15 Minutes")
plt.savefig('temporal.png')
plt.close()
print "done."

# ## SPATIOTEMPORAL INTERPOLATION ##
print "Generating plot of spatiotemporally interpolated (bilinear, linear) data...",
sti_grid_tec, sti_grid_rms, satellite_biases, station_biases = read_ionex.interpolate2D_spatiotemporal(FILE, temporal_grid = temporal_grid, spatial_grid = (dense_lons, dense_lats), \
													spatial_method = 'linear', temporal_method = 'linear')
plot_grid2D(dense_lons, dense_lats, sti_grid_tec[:,:,-1], temporal_grid[-1], title_label = "Spatiotemporal Interpolation (Bilinear, Linear); dLat, dLon = 1; dt = 15 Minutes")
plt.savefig('spatiotemporal.png')
plt.close()
print "done."
