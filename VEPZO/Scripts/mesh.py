import math
import datetime

import numpy as np
import pandas as pd

import shapely
from shapely.geometry import Polygon, Point, LineString
from shapely.geometry import box
from matplotlib import pyplot as plt
import matplotlib.animation as animation

# solar calculator
def getSEA(struct_time, latitude, longitude, utc):
    date = struct_time
    # print('date', date)
    hour = date[3]
    minute = date[4]
    # Check your timezone to add the offset
    hour_minute = (hour + minute / 60) - utc
    day_of_year = date[7]

    g = (360 / 365.25) * (day_of_year + hour_minute / 24)

    g_radians = math.radians(g)

    declination = 0.396372 - 22.91327 * math.cos(g_radians) + 4.02543 * math.sin(g_radians) - 0.387205 * math.cos(
        2 * g_radians) + 0.051967 * math.sin(2 * g_radians) - 0.154527 * math.cos(3 * g_radians) + 0.084798 * math.sin(
        3 * g_radians)

    time_correction = 0.004297 + 0.107029 * math.cos(g_radians) - 1.837877 * math.sin(g_radians) - 0.837378 * math.cos(
        2 * g_radians) - 2.340475 * math.sin(2 * g_radians)

    SHA = (hour_minute - 12) * 15 + longitude + time_correction

    if SHA > 180:
        SHA_corrected = SHA - 360
    elif SHA < -180:
        SHA_corrected = SHA + 360
    else:
        SHA_corrected = SHA

    lat_radians = math.radians(latitude)
    d_radians = math.radians(declination)
    SHA_radians = math.radians(SHA)

    SZA_radians = math.acos(
        math.sin(lat_radians) * math.sin(d_radians) + math.cos(lat_radians) * math.cos(d_radians) * math.cos(
            SHA_radians))

    SZA = math.degrees(SZA_radians)

    SEA = 90 - SZA

    return SEA

def getAZ(struct_time, latitude, longitude, utc):
    date = struct_time
    hour = date[3]
    minute = date[4]
    # Check your timezone to add the offset
    hour_minute = (hour + minute / 60) - utc
    day_of_year = date[7]

    g = (360 / 365.25) * (day_of_year + hour_minute / 24)

    g_radians = math.radians(g)

    declination = 0.396372 - 22.91327 * math.cos(g_radians) + 4.02543 * math.sin(g_radians) - 0.387205 * math.cos(
        2 * g_radians) + 0.051967 * math.sin(2 * g_radians) - 0.154527 * math.cos(3 * g_radians) + 0.084798 * math.sin(
        3 * g_radians)

    time_correction = 0.004297 + 0.107029 * math.cos(g_radians) - 1.837877 * math.sin(g_radians) - 0.837378 * math.cos(
        2 * g_radians) - 2.340475 * math.sin(2 * g_radians)

    SHA = (hour_minute - 12) * 15 + longitude + time_correction

    if SHA > 180:
        SHA_corrected = SHA - 360
    elif SHA < -180:
        SHA_corrected = SHA + 360
    else:
        SHA_corrected = SHA

    lat_radians = math.radians(latitude)
    d_radians = math.radians(declination)
    SHA_radians = math.radians(SHA)

    SZA_radians = math.acos(
        math.sin(lat_radians) * math.sin(d_radians) + math.cos(lat_radians) * math.cos(d_radians) * math.cos(
            SHA_radians))

    SZA = math.degrees(SZA_radians)

    cos_AZ = (math.sin(d_radians) - math.sin(lat_radians) * math.cos(SZA_radians)) / (
    math.cos(lat_radians) * math.sin(SZA_radians))

    AZ_rad = math.acos(cos_AZ)
    AZ = math.degrees(AZ_rad)

    if hour >= 12:
        AZ = 360 - AZ

    return AZ

def FitDimension(edge, dim):
	floor_value = math.floor(edge / dim)
	ceil_value = math.ceil(edge / dim)
	# if edge / dim is exact division
	if floor_value == ceil_value:
		return floor_value
	floor_check = abs(edge / floor_value - dim)
	ceil_check = abs(edge / ceil_value - dim)
	if floor_check < ceil_check:
		return floor_value
	else:
		return ceil_value

def PileUpList(flatlist, dimx, dimy):
	reclist = []
	for i in range(dimx * dimy):
		if i < len(flatlist):
			reclist.append(flatlist[i])
		else:
			reclist.append(None)
	nests = []
	sub = []
	for i in range(len(reclist)):
		sub.append(reclist[i])
		if len(sub) == dimx:
			nests.insert(0, sub)
			sub = []
	return nests

# wall dimension:				______
#  ________ 1111111111111111111  dimz
#	win_h	1111111110000001111
#			1111111110000001111
#  ________	1111111110000001111
#	sill	ptA1111111111111ptB
#			|- gap -|
def ProjectWall(ptA, ptB, zelta, theta, gap, sill, win_h)->" \
	Returns the shadow (shapely.geometry.Polygon) of the wall under the sun at \
	certain azimuth (zelta) and altitude (theta, elevation angle from horizon). \
	Returns tao, the incident beam angle to the wall surface":
	length = math.sqrt(
		math.pow(ptB[0] - ptA[0], 2) + 
		math.pow(ptB[1] - ptA[1], 2))
	dirx = (ptB[0] - ptA[0]) / length
	diry = (ptB[1] - ptA[1]) / length
	# dir0 = math.acos(diry)
	dir0 = math.atan2(1, 0) - math.atan2(diry, dirx)
	if dir0 < 0:
		dir0 += math.pi * 2
	tao = abs(zelta - dir0)
	if tao > math.pi:
		tao -= math.pi
	# print(zeltas[i]/math.pi * 180)
	# print(dir0/math.pi * 180)
	# print(tao/math.pi * 180)

	# calculate unit projector with unit length
	projx = math.sin(zelta - math.pi) / math.tan(theta)
	projy = math.cos(zelta - math.pi) / math.tan(theta)
	# consider the sill, gap and wwr to construct the projection
	pt1 = Point(ptA[0] + dirx * gap + projx * sill, 
		ptA[1] + diry * gap + projy * sill)
	pt2 = Point(ptB[0] - dirx * gap + projx * sill, 
		ptB[1] - diry * gap + projy * sill)
	pt3 = Point(ptB[0] - dirx * gap + projx * (sill + win_h), 
		ptB[1] - diry * gap + projy * (sill + win_h))
	pt4 = Point(ptA[0] + dirx * gap + projx * (sill + win_h), 
		ptA[1] + diry * gap + projy * (sill + win_h))
	shadow = Polygon([pt1, pt2, pt3, pt4, pt1])

	return shadow, tao


def GetMesh(
	polyloops:"nested lists of vertice coords representing exterior/interior boundary", 
	mask_wwr:"list of window to wall ratio of each wall. duplicate the last if not enough", 
	mask_adia:"list of boundary condition of each wall. duplicate the last if not enough", 
	scale_factor:"the default mesh cell will be 1m * 1m times this value", 
	location:"tuple for your location: (latitude, longitude, utc)", 
	stopwatch:"tuple for simulation time: (time_start, time_end, time_delta)", 
	solarload:"solar radiation, apply a discount on 1380 W/m2 according to lattitude and clearness")->" \
	Returns snapshots, a nested lists representing solar load of each time step. \
	Returns mesh_mask, a 2d array representing floorplan mesh and the boundary condition. \
	Returns cell_dim, a tuple of cell dimension: (x, y)":

	# unzip polyloops to Shapely Polygon
	# the polyloops must be a nested loops (list)
	floorplan = Polygon(polyloops[0])
	# if more than 1 loop, update the floorplan with holes in it
	if len(polyloops) > 1:
		innerpolys = []
		for i in range(1, len(polyloops)):
			innerpolys.append(Polygon(polyloops[i]))
		floorplan = Polygon(floorplan.exterior.coords, \
			[poly.exterior.coords for poly in innerpolys])

	while len(mask_wwr) < len(polyloops[0]):
		mask_wwr.append(mask_wwr[len(mask_wwr) - 1])
	while len(mask_adia) < len(polyloops[0]):
		mask_adia.append(0)

	dimz = 3
	gap = 0.2
	sill = 1

	# x-y bounding box is a (minx, miny, maxx, maxy) tuple
	dimx = floorplan.bounds[2] - floorplan.bounds[0]
	dimy = floorplan.bounds[3] - floorplan.bounds[1]
	tickx = FitDimension(dimx, scale_factor)
	ticky = FitDimension(dimy, scale_factor)
	# update floorplan by moving it to the origin
	floorplan = shapely.affinity.translate(
		floorplan, -floorplan.bounds[0], -floorplan.bounds[1])
	# move out a little bit for the boundary mesh
	floorplan = shapely.affinity.translate(floorplan, 
		dimx / tickx, dimy / ticky)
	# cache shell/holes list for coords tuples
	# convert shapely.coords.CoordinateSequence obj to list
	outerloop = list(floorplan.exterior.coords)
	print(outerloop)
	outerpatch = Polygon(floorplan.exterior)
	outeredge_adia = []
	for i in range(len(outerloop) - 1):
		if mask_adia[i] == 1:
			outeredge_adia.append(LineString([outerloop[i], outerloop[i + 1]]))
		else:
			outeredge_adia.append(None)
	innerloops = []
	innerpatchs = []
	for hole in floorplan.interiors:
		innerloops.append(list(hole.coords))
		innerpatchs.append(Polygon(hole))
	

	# generate the mask
	mesh = []
	cell_area = dimx / tickx * dimy / ticky

	# note that a perimeter mesh has grown based on the actual floorplan
	# to represent the boundary condition that can be customized by user
	# so the tickx and ticky plus 2. for example:
	#							0 0 0 0 0
	#   1 0 0					0 1 0 0 0
	#   1 1 0   grows into ->	0 1 1 0 0
	#   1 1 1 					0 1 1 1 0
	#							0 0 0 0 0
	# OUTPUT VARIABLE
	mask_mesh = np.zeros((ticky + 2, tickx + 2), dtype=int)

	# grow a matrix from bottom to top (positive y axis)
	# note that the mesh follows the same order
	for i in range((tickx + 2) * (ticky + 2)):
	    col = i % (tickx + 2)
	    row = (ticky + 2) - i // (tickx + 2) - 1
	    cell = box(col * dimx / tickx, row * dimy / ticky, \
	        (col + 1) * dimx / tickx, (row + 1) * dimy / ticky)
	    mesh.append(cell)
	    # note that this returns the intersection of cell and the mcr
	    sect_positive = outerpatch.intersection(cell)
	    if sect_positive.area / cell_area > 0.5:
	    	mask_mesh[(ticky + 2) - row - 1][col] += 1
	    	for innerpatch in innerpatchs:
	    		sect_negative = innerpatch.intersection(cell)
	    		if sect_negative.area / cell_area > 0.5:
	    			mask_mesh[(ticky + 2) - row - 1][col] += 1
	    else:
	    	for edge in outeredge_adia:
	    		if edge is not None:
	    			sect_line = cell.intersection(edge)
		    		if sect_line.length > 0.000001:
		    			mask_mesh[(ticky + 2) - row - 1][col] = 3

	# OUTPUT VARIABLE
	cell_dimension = (dimx / tickx, dimy / ticky)

	print(mask_mesh)

	# sum(sum(mask_mesh))
	mask_mesh_1d = mask_mesh.flatten()

	# projection vector
	zeltas = []         # in degree North-0 East-90 South-180 West-270
	thetas = []         # in degree Horizontal-0 Perpendicular-90
	time_ellapse = 0    # accumulated time in second
	time_current = stopwatch[0]
	while time_current < stopwatch[1]:
	    time_current += stopwatch[2]
	    thetas.append(math.radians(getSEA(
	    	time_current.timetuple(), location[0], location[1], location[2])))
	    zeltas.append(math.radians(getAZ(
	    	time_current.timetuple(), location[0], location[1], location[2])))
	# print(thetas)
	# print(zeltas)

	# a 1d array that recreates the matrix from top to bottom (negative y axis)
	mask_mesh_compressed = mask_mesh.reshape(1, (tickx + 2) * (ticky + 2))
	# initialize zeros for solar mask at each snapshot
	frame_solar = []
	checker = ""
	for i in range(len(zeltas)):
		# print("timestep: " + str(i))

		mask_solar = [ 0 for j in range(len(mesh))]

		# skip if the sun still below the horizon
		if thetas[i] < 0:
			frame_solar.append(PileUpList(mask_solar, tickx + 2, ticky + 2))
			continue

		# iterate each edge of the polygon
		for v in range(len(outerloop) - 1): # len(outerloop) - 1
			if mask_wwr[v] == 0:
				continue

			length = math.sqrt(
				math.pow(outerloop[v + 1][0] - outerloop[v][0], 2) + 
				math.pow(outerloop[v + 1][1] - outerloop[v][1], 2))
			win_height = length * dimz * mask_wwr[v] / (length - 2 * gap) - sill
			(shadow, tao) = ProjectWall(outerloop[v], outerloop[v + 1], 
				zeltas[i], thetas[i], gap, sill, win_height)
			shadow_area = shadow.area

			# checker += "{{{0:2f}, {1:2f}, 0}}\n".format(ptA.x, ptA.y) \
			# 	+ "{{{0:2f}, {1:2f}, 0}}\n".format(ptB.x, ptB.y) \
			# 	+ "{{{0:2f}, {1:2f}, 0}}\n".format(ptC.x, ptC.y) \
			# 	+ "{{{0:2f}, {1:2f}, 0}}\n".format(ptD.x, ptD.y)

			# filter out invalid/self-intersected polygon
			if shadow_area <= 0.000001:
				continue

			# calculate the equivalent area of solar beam
			area = (length - 2 * gap) * math.sin(tao) * win_height * math.cos(thetas[i])
			# print(area)
			# print("-----")

			if area < 0.000001:
				continue

			# remove blocking from self shadowing by outer loop
			# pending work

			# remove blocking from self shadowing by inner loop
			for w in range(len(outerloop) - 1):
				if w != v:
					(_shadow, _tao) = ProjectWall(outerloop[w], outerloop[w + 1], 
						zeltas[i], thetas[i], 0, 0, dimz - sill)
					shadow = shadow.difference(_shadow)

			# remove blocking from inner loops shadowing
			for innerloop in innerloops:
				for w in range(len(innerloop) - 1):
					(_shadow, _tao) = ProjectWall(innerloop[w], innerloop[w + 1], 
						zeltas[i], thetas[i], 0, 0, dimz - sill)
					shadow = shadow.difference(_shadow)

			# iterate each cell, following the mask_mesh_1d
			for j in range(len(mesh)):
				if mask_mesh_1d[j] != 1:
					mask_solar[j] += 0
					continue

				sect = shadow.intersection(mesh[j])
				if sect.area < 0.000001:
					continue

				# for pt in list(sect.exterior.coords):
				# 	checker += "{{{0:2f}, {1:2f}, 0}}# ".format(pt[0], pt[1])
				# checker += "\n"

				mask_solar[j] += (sect.area / shadow_area) * area * solarload

		frame_solar.append(PileUpList(mask_solar, tickx + 2, ticky + 2))


	# print(checker)

	return frame_solar, mask_mesh, cell_dimension


if __name__ =='__main__':

	# DIAMOND
	vertexloops = [[[3, 0], [10, 0], [10, 7], [7, 10], [0, 10], [0, 3], [3, 0]], 
		 [[4, 4], [6, 4], [6, 6], [4, 6], [4, 4]]]

	# RAND
	# vertexloops = [[[4, 0], [10, 0], [10, 5], [6, 5], [6, 10], [0, 10], [0, 5], [4, 5], [4, 0]]]

	# SNAKE
	# vertexloops = [[[0, 0], [11, 0], [11, 3], [2, 3], [2, 4], [11, 4], [11, 11], [0, 11], [0, 8], \
	#    [9, 8], [9, 7], [0, 7], [0, 0]]]

	# FRAME
	# vertexloops = [[[0, 0], [10, 0], [10, 10], [0, 10], [0, 0]], \
	# 	[[3, 3], [7, 3], [7, 7], [3, 7], [3, 3]]]

	# vertexloops = [[[0, 0], [10, 0], [10, 10], [0, 10], [0, 0]]]

	# mask_wwr = [0.8]
	mask_wwr = [0.8, 0]
	mask_adia = [0, 0, 0, 1, 1, 0]

	scale_factor = 1

	latitude = 33
	longitude = 122
	utc = 8      # time zone

	time_step = 600					# in seconds
	time_start = datetime.datetime.strptime("2022-01-01 00:00:00", "%Y-%m-%d %H:%M:%S")
	time_end = datetime.datetime.strptime("2022-01-01 23:59:59", "%Y-%m-%d %H:%M:%S")
	time_delta = datetime.timedelta(seconds=time_step)

	solarload = 500

	(snapshots, mask_mesh, mesh_dim) = GetMesh(
		vertexloops, mask_wwr, mask_adia, scale_factor, 
		(latitude, longitude, utc), 
		(time_start, time_end, time_delta), 
		solarload)

	# visualize radiation intensity on the floorplan
	maxRadiation = 0
	for snapshot in snapshots:
		if np.max(snapshot) > maxRadiation:
			maxRadiation = np.max(snapshot)

	fig = plt.figure(figsize=(5,5), dpi=72)
	xlattice = np.arange(0, len(snapshots[0][0]) + 1)
	ylattice = np.arange(0, len(snapshots[0]) + 1)
	cmap = plt.cm.get_cmap('inferno').copy()
	cmap.set_bad(color = 'w', alpha = 1.)
	plt.axis('equal')
	plt.axis('off') # remove all axes
	plt.xticks([]) # remove ticks on x-axis
	plt.yticks([]) # remove ticks on y-axis
	plt.subplots_adjust(left=0, bottom=0, right=1, top=1)
	imgs = []
	for snapshot in snapshots:
		for i in range(len(snapshot)):
			for j in range(len(snapshot[0])):
				if mask_mesh[len(snapshot) - i - 1][j] != 1:
					snapshot[i][j] = np.nan
		snapshot = np.ma.masked_invalid(snapshot)
		# ax = plt.pcolor(xlattice, ylattice, snapshot, norm=plt.Normalize(0, maxRadiation))
		ax = plt.pcolormesh(xlattice, ylattice, snapshot, cmap=cmap, edgecolors='None', 
			norm=plt.Normalize(0, maxRadiation))
		imgs.append((ax,))
	img_ani = animation.ArtistAnimation(fig, imgs, interval=50, repeat_delay=0, blit=True)
	plt.show()

	# pip install pillow
	# https://stackoverflow.com/questions/51512141/how-to-make-matplotlib-saved-gif-looping
	class LoopingPillowWriter(animation.PillowWriter):
		def finish(self):
			self._frames[0].save(
				self.outfile, save_all=True, append_images=self._frames[1:],
					duration=int(1000 / self.fps), loop=0)
	img_ani.save(r"animation.gif", writer=LoopingPillowWriter(fps=20))
	# img_ani.save(r"animation.gif", writer=animation.PillowWriter(fps=30))
	# img_ani.save(r"animation.gif", writer='imagemagick')