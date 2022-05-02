import math
import datetime

import numpy as np
import pandas as pd

import shapely
from shapely.geometry import Polygon, Point
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

    if (SHA > 180):
        SHA_corrected = SHA - 360
    elif (SHA < -180):
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

    if (SHA > 180):
        SHA_corrected = SHA - 360
    elif (SHA < -180):
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

def GetMesh(polyloop, wwr, cell_dim, location, stopwatch, solarload):

	floorplan = Polygon(polyloop)

	while len(wwr) < len(polyloop):
		wwr.append(wwr[len(wwr) - 1])

	dimz = 3
	gap = 0.2
	sill = 0.5

	# x-y bounding box is a (minx, miny, maxx, maxy) tuple
	dimx = floorplan.bounds[2] - floorplan.bounds[0]
	dimy = floorplan.bounds[3] - floorplan.bounds[1]
	tickx = FitDimension(dimx, cell_dim)
	ticky = FitDimension(dimy, cell_dim)
	floorplan_offset = shapely.affinity.translate(
		floorplan, -floorplan.bounds[0], -floorplan.bounds[1])
	shell = list(floorplan_offset.exterior.coords)

	# generate the mask
	mesh = []
	cell_area = dimx / tickx * dimy / ticky
	# OUTPUT VARIABLE
	mask_mesh = np.zeros((ticky, tickx), dtype=int)
	# grow a matrix from bottom to top (positive y axis)
	# note that the mesh follows the same order
	for i in range(tickx * ticky):
	    col = i % tickx
	    row = ticky - i // tickx - 1
	    cell = box(col * dimx / tickx, row * dimy / ticky, \
	        (col + 1) * dimx / tickx, (row + 1) * dimy / ticky)
	    mesh.append(cell)
	    sect = floorplan_offset.intersection(cell)
	    if sect.area / cell_area > 0.5:
	    	mask_mesh[ticky - row - 1][col] = 1
	# OUTPUT VARIABLE
	mesh_dimension = (dimx, dimy)

	# print(mask_mesh)

	# trim the mask
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
	mask_mesh_compressed = mask_mesh.reshape(1, tickx * ticky)
	# initialize zeros for solar mask at each snapshot
	frame_solar = []
	checker = ""
	for i in range(len(zeltas)):
		# print("timestep: " + str(i))

		mask_solar = [ 0 for j in range(len(mesh))]

		# skip if the sun still below the horizon
		if thetas[i] < 0:
			frame_solar.append(PileUpList(mask_solar, tickx, ticky))
			continue

		# iterate each edge of the polygon
		for v in range(len(shell) - 1): # len(shell) - 1
			if wwr[v] == 0:
				continue

			length = math.sqrt(
				math.pow(shell[v + 1][0] - shell[v][0], 2) + 
				math.pow(shell[v + 1][1] - shell[v][1], 2))
			height = length * dimz * wwr[v] / (length - 2 * gap) - sill
			dirx = (shell[v + 1][0] - shell[v][0]) / length
			diry = (shell[v + 1][1] - shell[v][1]) / length
			# dir0 = math.acos(diry)
			dir0 = math.atan2(1, 0) - math.atan2(diry, dirx)
			if dir0 < 0:
				dir0 += math.pi * 2
			tao = abs(zeltas[i] - dir0)
			if tao > math.pi:
				tao -= math.pi
			# print(zeltas[i]/math.pi * 180)
			# print(dir0/math.pi * 180)
			# print(tao/math.pi * 180)

			# calculate unit projector with unit length
			projx = math.sin(zeltas[i] - math.pi) / math.tan(thetas[i])
			projy = math.cos(zeltas[i] - math.pi) / math.tan(thetas[i])
			# consider the sill, gap and wwr to construct the projection
			ptA = Point(shell[v][0] + dirx * gap + projx * sill, 
				shell[v][1] + diry * gap + projy * sill)
			ptB = Point(shell[v + 1][0] - dirx * gap + projx * sill, 
				shell[v + 1][1] - diry * gap + projy * sill)
			ptC = Point(shell[v + 1][0] - dirx * gap + projx * (sill + height), 
				shell[v + 1][1] - diry * gap + projy * (sill + height))
			ptD = Point(shell[v][0] + dirx * gap + projx * (sill + height), 
				shell[v][1] + diry * gap + projy * (sill + height))
			shadow = Polygon([ptA, ptB, ptC, ptD, ptA])

			# checker += "{{{0:2f}, {1:2f}, 0}}\n".format(ptA.x, ptA.y) \
			# 	+ "{{{0:2f}, {1:2f}, 0}}\n".format(ptB.x, ptB.y) \
			# 	+ "{{{0:2f}, {1:2f}, 0}}\n".format(ptC.x, ptC.y) \
			# 	+ "{{{0:2f}, {1:2f}, 0}}\n".format(ptD.x, ptD.y)

			# filter out invalid/self-intersected polygon
			if shadow.area <= 0.000001:
				continue
			area = (length - 2 * gap) * math.sin(tao) * height * math.cos(thetas[i])
			# print(area)
			# print("-----")

			if area < 0.000001:
				continue

			# iterate each cell, following the mask_mesh_1d
			for j in range(len(mesh)):
				if mask_mesh_1d[j] == 0:
					mask_solar[j] += 0
					continue

				sect = shadow.intersection(mesh[j])
				if sect.area < 0.000001:
					continue

				# for pt in list(sect.exterior.coords):
				# 	checker += "{{{0:2f}, {1:2f}, 0}}# ".format(pt[0], pt[1])
				# checker += "\n"

				mask_solar[j] += sect.area / area * solarload

		frame_solar.append(PileUpList(mask_solar, tickx, ticky))


	# print(checker)

	return frame_solar, mask_mesh, mesh_dimension


if __name__ =='__main__':

	# vertices = [[2, 1], [8, 1], [8, 4], [4, 4], [4, 7], [2, 7], [2, 1]]
	vertices = [[0, 0], [10, 0], [10, 10], [0, 10], [0, 0]]
	
	mask_wwr = [0.8, 0, 0, 0]

	cell_dim = 2

	latitude = 33
	longitude = 122
	utc = 8      # time zone

	time_start = datetime.datetime.strptime("2022-01-01 00:00:00", "%Y-%m-%d %H:%M:%S")
	time_end = datetime.datetime.strptime("2022-01-01 23:59:59", "%Y-%m-%d %H:%M:%S")
	time_delta = datetime.timedelta(minutes=10)

	solarload = 1000

	(snapshots, mask_mesh, mesh_dim) = GetMesh(
		vertices, mask_wwr, cell_dim, 
		(latitude, longitude, utc), 
		(time_start, time_end, time_delta), 
		solarload)

	# visualize radiation intensity on the floorplan
	maxRadiation = 0
	for snapshot in snapshots:
		if np.max(snapshot) > maxRadiation:
			maxRadiation = np.max(snapshot)

	fig = plt.figure()
	xlattice = np.arange(0, len(snapshots[0][0]) + 1)
	ylattice = np.arange(0, len(snapshots[0]) + 1)
	imgs = []
	for snapshot in snapshots:
		imgs.append((plt.pcolor(xlattice, ylattice, snapshot, 
			norm=plt.Normalize(0, maxRadiation)),))
	img_ani = animation.ArtistAnimation(fig, imgs, interval = 50, repeat_delay=0, blit=True)
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