"""
Read/write tools for nonuniform electric field .grd format.
Matthew Grawe, grawe2 (at) illinois.edu
January 2017
"""

import numpy as np

def next_line(grd_file):
	"""
	next_line
	Function returns the next line in the file
	that is not a blank line, unless the line is
	'', which is a typical EOF marker.
	"""

	done   = False
	while not done:
		line = grd_file.readline()
		if line == '':
			return line, False
		elif line.strip():
			return line, True


def read_block(grd_file, n_lats):

	lats = []

	# read+store until we have collected n_lats
	go = True
	while go:

		fline, status = next_line(grd_file)
		line          = fline.split()

		# the line hats lats in it
		lats.extend(np.array(line).astype('float'))

		if len(lats) == 17:
			go = False
			return np.array(lats)

def grd_read(grd_filename):
	"""
	Opens the .grd file grd_file and returns the following:
	lon_grid :  1D numpy array of lons
	lat_grid :  1D numpy array of lats
	time_grid:  1D numpy array of times
	DATA     :  3D numpy array of the electric field data, such that
				the electric field at (lon, lat) for time t
				is accessed via DATA[lon, lat, t].
	"""
	with open(grd_filename, 'rb') as grd_file:

		# read the header line
		fline, status = next_line(grd_file)
		line          = fline.split()

		lon_res      = float(line[0])
		lon_west     = float(line[1])
		n_lons       = int(line[2])
		lat_res      = float(line[3])
		lat_south    = float(line[4])
		n_lats       = int(line[5])

		DATA  = []
		times = []

		go = True
		while go:
			# get the time index line
			ftline, status = next_line(grd_file)
			tline          = ftline.split()
			t              = float(tline[0])

			times.append(t)

			SLICE = np.zeros([n_lons, n_lats])
			for lon_index in range(0, n_lons):

				data_slice = read_block(grd_file, n_lats)

				SLICE[lon_index, :] = data_slice

			DATA.append(SLICE.T)

			# current line should have length one to indicate next time index
			# make sure, then back up

			before        = grd_file.tell()
			fline, status = next_line(grd_file)
			line          = fline.split()

			if len(line) != 1:
				if status == False:
					# EOF, leave
					break
				else:
					raise Exception('Unexpected number of lat entries.')

			grd_file.seek(before)

	DATA      = np.array(DATA).T
	lon_grid  = np.arange(lon_west,  lon_west  + lon_res*n_lons, lon_res)
	lat_grid  = np.arange(lat_south, lat_south + lat_res*n_lats, lat_res)
	time_grid = np.array(times)

	return lon_grid, lat_grid, times, DATA

def write_lon_block(grd_file, n_lats, data):
	"""
	len(data) == n_lats should be True
	"""

	current_index =  0
	go1 = True
	while go1:
		line = ['']*81
		go2  = True
		internal_index = 0
		while go2:
			datum = data[current_index]

			line[16*internal_index:16*internal_index+16] = str(datum).rjust(16)
			current_index += 1
			internal_index += 1			

			if(current_index >= len(data)):
				line[80] = '\n'
				grd_file.write("".join(line))
				go2 = False
				go1 = False
			elif(internal_index >= 5):
				line[80] = '\n'
				grd_file.write("".join(line))
				go2 = False

def grd_write(grd_filename, lon_grid, lat_grid, time_grid, DATA):
	"""
	Writes out DATA corresponding to the locations
	specified by lon_grid, lat_grid in the .grd format.

	lon_grid must have the westmost point as lon_grid[0].
	lat_grid must have the southmost point as lat_grid[0].

	Assumptions made:
	latitude/longitude resolutions are positive
	number of latitude/longitude points in header is positive
	at least one space must be between each number
	data lines have no more than 5 entries

	Assumed structure of header line:
	# first 16 spaces allocated as longitude resolution
	# next 16 spaces allocated as westmost longitude
	# next 5 spaces allocated as number of longitude points
	# next 11 spaces allocated as latitude resolution
	# next 16 spaces allocated as southmost latitude
	# next 16 spaces allocated for number of latitude points
	# TOTAL: 80 characters

	Assumed stucture of time line:
	# 5 blank spaces
	# rest of line allocated for time

	Assumed structure a data line:
	# 16 spaces allocated for data entry
	# .. .. ..

	"""

	with open(grd_filename, 'wb') as grd_file:

		lon_res      = np.abs(lon_grid[1] - lon_grid[0])
		lon_west     = lon_grid[0]
		n_lons       = len(lon_grid)
		lat_res      = np.abs(lat_grid[1] - lat_grid[0])
		lat_south    = lat_grid[0]
		n_lats       = len(lat_grid)
		n_times      = len(time_grid)

		# write the header: 80 characters
		header = ['']*81
		header[0:16]  = str(lon_res).rjust(16)
		header[16:32] = str(lon_west).rjust(16)
		header[32:37] = str(n_lons).rjust(5)
		header[37:48] = str(lat_res).rjust(11)
		header[48:64] = str(lat_south).rjust(16)
		header[64:80] = str(n_lats).rjust(16)
		header[80]    = '\n'

		header_str = "".join(header)
		grd_file.write(header_str)

		for i, t in enumerate(time_grid):
			# write the time line
			timeline = ['']*14
			timeline[5:] = str(t).rjust(9)
			timeline[13] = '\n'
			timeline_str = "".join(timeline)
			grd_file.write(timeline_str)

			for j, lon in enumerate(lon_grid):

				# write the lon blocks
				write_lon_block(grd_file, n_lats, DATA[j, :, i])

		grd_file.close()
