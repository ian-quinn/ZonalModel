import math
from datetime import datetime, timedelta

import numpy as np
from shapely.geometry import Polygon, Point, LineString
from shapely.geometry import box
from shapely import affinity
from matplotlib import pyplot as plt
import matplotlib.animation as animation

__tol__ = 0.000001


# solar calculator
def get_sun_position(time: datetime, location: tuple[float, float, int]) -> tuple[float, float]:
    latitude = location[0]
    longitude = location[1]
    utc = location[2]

    date = time.timetuple()
    # in case you need to parse the timestamp first to create a datetime
    hour = date[3]
    minute = date[4]
    # Check your timezone to add the offset
    hour_minute = (hour + minute / 60) - utc
    day_of_year = date[7]

    g = (360 / 365.25) * (day_of_year + hour_minute / 24)

    g_rad = math.radians(g)

    declination = 0.396372 - 22.91327 * math.cos(g_rad) + 4.02543 * math.sin(g_rad) - 0.387205 * math.cos(
        2 * g_rad) + 0.051967 * math.sin(2 * g_rad) - 0.154527 * math.cos(3 * g_rad) + 0.084798 * math.sin(
        3 * g_rad)

    time_correction = 0.004297 + 0.107029 * math.cos(g_rad) - 1.837877 * math.sin(g_rad) - 0.837378 * math.cos(
        2 * g_rad) - 2.340475 * math.sin(2 * g_rad)

    # solar hour angle measures the sun position relative to the observer on earth
    # (-180, 180) -180: sunrise, 0: noon, 180: sunset
    solar_hour_angle = (hour_minute - 12) * 15 + longitude + time_correction
    solar_hour_angle_corrected = solar_hour_angle
    if solar_hour_angle > 180:
        solar_hour_angle_corrected = solar_hour_angle - 360
    elif solar_hour_angle < -180:
        solar_hour_angle_corrected = solar_hour_angle + 360

    lat_rad = math.radians(latitude)
    dec_rad = math.radians(declination)
    sha_rad = math.radians(solar_hour_angle_corrected)

    # zenith measures the angle from sun position to the vertical top
    # (0, 90) 0: right above head, 90: at horizon
    zenith_rad = math.acos(math.sin(lat_rad) * math.sin(dec_rad) +
                           math.cos(lat_rad) * math.cos(dec_rad) * math.cos(sha_rad))
    zenith = math.degrees(zenith_rad)
    elevation = 90 - zenith

    # azimuth angle runs clockwise from NORTH on horizontal plane
    # (0, 360) 0: north, 90: east, 180: south, 270: west
    azimuth_cos = ((math.sin(dec_rad) - math.sin(lat_rad) * math.cos(zenith_rad)) /
                   (math.cos(lat_rad) * math.sin(zenith_rad)))
    azimuth_rad = math.acos(azimuth_cos)
    azimuth = math.degrees(azimuth_rad)
    if hour >= 12:
        azimuth = 360 - azimuth

    return elevation, azimuth


def fit_dimension(edge, dim):
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


# use np.reshape() for a workflow with np.array
def pileup_list(flatlist, dim_x, dim_y):
    reclist = []
    for i in range(dim_x * dim_y):
        if i < len(flatlist):
            reclist.append(flatlist[i])
        else:
            reclist.append(None)
    nests = []
    sub = []
    for i in range(len(reclist)):
        sub.append(reclist[i])
        if len(sub) == dim_x:
            nests.insert(0, sub)
            sub = []
    return nests


"""
 wall dimension:				______
  ________ 1111111111111111111  dimz
    win_h  1111111110000001111
           1111111110000001111
  ________ 1111111110000001111
    sill   ptA1111111111111ptB
           |- gap -|
"""


def project_wall(pt_a, pt_b, zeta, theta, gap, sill, win_h):
    """
    :param pt_a: start point of wall location line
    :param pt_b: end point of wall location line
    :param zeta: solar azimuth angle
    :param theta: solar altitude angle
    :param gap: distance from window to ptA
    :param sill: height of window sill
    :param win_h: height of window from wall bottom
    :return:
    shadow <shapely.geometry.Polygon> the shadow of the wall under the sun at certain azimuth (zeta)
    and altitude (theta, elevation angle from horizon).
    tao <number>, the incident beam angle to the wall surface
    """
    length = math.sqrt(
        math.pow(pt_b[0] - pt_a[0], 2) +
        math.pow(pt_b[1] - pt_a[1], 2))
    dir_x = (pt_b[0] - pt_a[0]) / length
    dir_y = (pt_b[1] - pt_a[1]) / length
    # dir0 = math.acos(dir_y)
    dir0 = math.atan2(1, 0) - math.atan2(dir_y, dir_x)
    if dir0 < 0:
        dir0 += math.pi * 2
    tao = abs(zeta - dir0)
    if tao > math.pi:
        tao -= math.pi
    # print(zetas[i]/math.pi * 180)
    # print(dir0/math.pi * 180)
    # print(tao/math.pi * 180)

    # calculate unit projector with unit length
    projx = math.sin(zeta - math.pi) / math.tan(theta)
    projy = math.cos(zeta - math.pi) / math.tan(theta)
    # consider the sill, gap and wwr to construct the projection
    pt1 = Point(pt_a[0] + dir_x * gap + projx * sill,
                pt_a[1] + dir_y * gap + projy * sill)
    pt2 = Point(pt_b[0] - dir_x * gap + projx * sill,
                pt_b[1] - dir_y * gap + projy * sill)
    pt3 = Point(pt_b[0] - dir_x * gap + projx * (sill + win_h),
                pt_b[1] - dir_y * gap + projy * (sill + win_h))
    pt4 = Point(pt_a[0] + dir_x * gap + projx * (sill + win_h),
                pt_a[1] + dir_y * gap + projy * (sill + win_h))
    shadow = Polygon([pt1, pt2, pt3, pt4, pt1])

    return shadow, tao


def get_mesh(polyloops, mask_wwr, mask_adiabatic, mesh_scale,
             location: tuple[float, float, int],
             stopwatch: tuple[datetime, datetime, timedelta],
             solarload, SHGC, epw_path: str, fakeslice: int) \
        -> tuple[np.array, np.array, tuple[int, int]]:
    """
    :param polyloops: nested lists of vertices coords representing exterior/interior boundary
    :param mask_wwr: list of window to wall ratio of each wall. duplicate the last if not enough
    :param mask_adiabatic: list of boundary condition of each wall. duplicate the last if not enough
    :param mesh_scale: the default mesh cell will be 1m * 1m times this value
    :param location: tuple for your location: (latitude, longitude, utc)
    :param stopwatch: tuple for simulation time: (time_start, time_end, time_delta)
    :param solarload: solar radiation, 1380 W/m2 * latitude * clearness * SHGC
    :param SHGC: solar heat gain coefficient
    :param epw_path: load solar beam & diffusion from TMY if provided
    :param fakeslice: > 0 to fake some typical days representing the whole year run
    :return:
    snapshots, a nested lists representing solar load of each time step.
    mesh_mask, a 2d array representing floorplan mesh and the boundary condition.
    cell_dim, a tuple of cell dimension: (x, y)
    """
    time_flag = datetime.strptime("2022-01-01 00:00:00", "%Y-%m-%d %H:%M:%S")

    # unzip polyloops to Shapely Polygon
    # the polyloops must be a nested loops (list)
    floorplan = Polygon(polyloops[0])
    # if more than 1 loop, update the floorplan with holes in it
    if len(polyloops) > 1:
        innerpolys = []
        for i in range(1, len(polyloops)):
            innerpolys.append(Polygon(polyloops[i]))
        floorplan = Polygon(floorplan.exterior.coords,
                            [poly.exterior.coords for poly in innerpolys])

    while len(mask_wwr) < len(polyloops[0]):
        mask_wwr.append(mask_wwr[len(mask_wwr) - 1])
    while len(mask_adiabatic) < len(polyloops[0]):
        mask_adiabatic.append(0)

    dim_z = 3
    gap = 0
    sill = 0.5

    # x-y bounding box is a (min_x, min_y, max_x, max_y) tuple
    dim_x = floorplan.bounds[2] - floorplan.bounds[0]
    dim_y = floorplan.bounds[3] - floorplan.bounds[1]
    tick_x = fit_dimension(dim_x, mesh_scale)
    tick_y = fit_dimension(dim_y, mesh_scale)
    # update floorplan by moving it to the origin
    floorplan = affinity.translate(
        floorplan, -floorplan.bounds[0], -floorplan.bounds[1])
    # move out a little bit for the boundary mesh
    floorplan = affinity.translate(floorplan,
                                   dim_x / tick_x, dim_y / tick_y)
    # cache shell/holes list for coords tuples
    # convert shapely.coords.CoordinateSequence obj to list
    outerloop = list(floorplan.exterior.coords)
    # print(outerloop)
    outerpatch = Polygon(floorplan.exterior)
    outeredge_adia = []
    for i in range(len(outerloop) - 1):
        if mask_adiabatic[i] == 1:
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
    cell_area = dim_x / tick_x * dim_y / tick_y

    """
    note that a perimeter mesh has grown based on the actual floorplan
    to represent the boundary condition that can be customized by user
    so the tick_x and tick_y plus 2. for example:
                            0 0 0 0 0
      1 0 0					0 1 0 0 0
      1 1 0   grows into ->	0 1 1 0 0
      1 1 1                 0 1 1 1 0
                            0 0 0 0 0
    """
    # OUTPUT VARIABLE
    mask_mesh = np.zeros((tick_y + 2, tick_x + 2), dtype=int)

    # grow a matrix from bottom to top (positive y-axis)
    # note that the mesh follows the same order
    for i in range((tick_x + 2) * (tick_y + 2)):
        col = i % (tick_x + 2)
        row = (tick_y + 2) - i // (tick_x + 2) - 1
        cell = box(col * dim_x / tick_x, row * dim_y / tick_y,
                   (col + 1) * dim_x / tick_x, (row + 1) * dim_y / tick_y)
        mesh.append(cell)
        # note that this returns the intersection of cell and the mcr
        sect_positive = outerpatch.intersection(cell)
        if sect_positive.area / cell_area > 0.5:
            mask_mesh[(tick_y + 2) - row - 1][col] += 1
            for innerpatch in innerpatchs:
                sect_negative = innerpatch.intersection(cell)
                if sect_negative.area / cell_area > 0.5:
                    mask_mesh[(tick_y + 2) - row - 1][col] += 1
        else:
            for edge in outeredge_adia:
                if edge is not None:
                    sect_line = cell.intersection(edge)
                    if sect_line.length > __tol__:
                        mask_mesh[(tick_y + 2) - row - 1][col] = 3

    # OUTPUT VARIABLE
    cell_dimension = (dim_x / tick_x, dim_y / tick_y)

    # print(mask_mesh)

    # sum(sum(mask_mesh))
    mask_mesh_1d = mask_mesh.flatten()

    # projection vector
    thetas = []  # in radians Horizontal-0 Perpendicular-90
    zetas = []  # in radians North-0 East-90 South-180 West-270
    # time_elapse = 0    # accumulated time in second

    # under fake mode, the sun position will be sampled evenly according to slice number
    # note that the June/December 22 will always be sampled (the highest/lowest) position
    if fakeslice:
        # convert 1-366 julian days to standard datetime
        sampledays = [datetime.strptime(str(int(366 / fakeslice) * (n + 1)), '%j') for n in range(fakeslice)]
        for date in sampledays:
            # offset a little bit to make June 22 at the center
            time_current = date - timedelta(days=7)
            time_stop = time_current + timedelta(hours=24)
            while time_current < time_stop:
                time_current += stopwatch[2]
                thetas.append(math.radians(get_sun_position(time_current, location)[0]))
                zetas.append(math.radians(get_sun_position(time_current, location)[1]))
    else:
        time_current = stopwatch[0]
        while time_current < stopwatch[1]:
            time_current += stopwatch[2]
            thetas.append(math.radians(get_sun_position(time_current, location)[0]))
            zetas.append(math.radians(get_sun_position(time_current, location)[1]))
    # print(thetas)
    # print(zetas)

    # if beam/diffuse radiation from TMY is loaded
    # overwrite the solarload variable
    series_beam = []
    series_time = []
    if epw_path:
        hour_tick = 0
        f = open(epw_path, mode="r")
        for line in f.readlines():
            datalist = line.split(',')
            if len(datalist) == 35:
                hour_tick += 1
                if hour_tick * 3600 > (stopwatch[0] - time_flag).total_seconds() and \
                        (hour_tick - 1) * 3600 < (stopwatch[1] - time_flag).total_seconds():
                    series_beam.append(int(datalist[13]))
                    if len(series_beam) == 0:
                        series_time.append(0.0)
                    else:
                        series_time.append((hour_tick - 1) * 3600 - (stopwatch[0] - time_flag).total_seconds())

    # a 1d array that recreates the matrix from top to bottom (negative y-axis)
    mask_mesh_compressed = mask_mesh.reshape(1, (tick_x + 2) * (tick_y + 2))
    # initialize zeros for solar mask at each snapshot
    frame_solar = []
    for i in range(len(zetas)):
        # print("timestep: " + str(i))

        mask_solar = [0 for j in range(len(mesh))]

        # skip if the sun still below the horizon
        if thetas[i] < 0:
            frame_solar.append(pileup_list(mask_solar, tick_x + 2, tick_y + 2))
            continue

        # iterate each edge of the polygon
        for v in range(len(outerloop) - 1):  # len(outerloop) - 1
            if mask_wwr[v] == 0:
                continue

            length = math.sqrt(
                math.pow(outerloop[v + 1][0] - outerloop[v][0], 2) +
                math.pow(outerloop[v + 1][1] - outerloop[v][1], 2))
            win_height = length * dim_z * mask_wwr[v] / (length - 2 * gap)
            (shadow, tao) = project_wall(outerloop[v], outerloop[v + 1],
                                         zetas[i], thetas[i], gap, sill, win_height)
            shadow_area = shadow.area

            # checker += "{{{0:2f}, {1:2f}, 0}}\n".format(ptA.x, ptA.y) \
            # 	+ "{{{0:2f}, {1:2f}, 0}}\n".format(ptB.x, ptB.y) \
            # 	+ "{{{0:2f}, {1:2f}, 0}}\n".format(ptC.x, ptC.y) \
            # 	+ "{{{0:2f}, {1:2f}, 0}}\n".format(ptD.x, ptD.y)

            # filter out invalid/self-intersected polygon
            if shadow_area <= __tol__:
                continue

            # calculate the equivalent area of solar beam
            area = (length - 2 * gap) * math.sin(tao) * win_height * math.cos(thetas[i])
            # print(area)
            # print("-----")

            if area < __tol__:
                continue

            # remove blocking from self shadowing by outer loop
            # pending work

            # remove blocking from self shadowing by inner loop
            for w in range(len(outerloop) - 1):
                if w != v:
                    (_shadow, _tao) = project_wall(outerloop[w], outerloop[w + 1],
                                                   zetas[i], thetas[i], 0, 0, dim_z - sill)
                    shadow = shadow.difference(_shadow)

            # remove blocking from inner loops shadowing
            for innerloop in innerloops:
                for w in range(len(innerloop) - 1):
                    (_shadow, _tao) = project_wall(innerloop[w], innerloop[w + 1],
                                                   zetas[i], thetas[i], 0, 0, dim_z - sill)
                    shadow = shadow.difference(_shadow)

            # iterate each cell, following the mask_mesh_1d
            for j in range(len(mesh)):
                if mask_mesh_1d[j] != 1:
                    mask_solar[j] += 0
                    continue

                sect = shadow.intersection(mesh[j])
                if sect.area < __tol__:
                    continue

                # for pt in list(sect.exterior.coords):
                # 	checker += "{{{0:2f}, {1:2f}, 0}}# ".format(pt[0], pt[1])
                # checker += "\n"

                # if use TMY and not in fake mode, overwrite the solarload
                if epw_path and not fakeslice:
                    for k in range(len(series_time)):
                        if i * stopwatch[2].total_seconds() < series_time[k]:
                            solarload = series_beam[k - 1]
                            # print("update solar as " + str(solarload))
                            break

                mask_solar[j] += round((sect.area / shadow_area) * area * solarload * SHGC, 2)

        frame_solar.append(pileup_list(mask_solar, tick_x + 2, tick_y + 2))

    # print(checker)
    # for snapshot in frame_solar:
    # print(np.array(snapshot))

    return np.array(frame_solar), mask_mesh, cell_dimension


# test this module
if __name__ == '__main__':
    '''
    DIAMOND
    POLYGON = [[[3, 0], [10, 0], [10, 7], [7, 10], [0, 10], [0, 3], [3, 0]], 
                [[4, 4], [6, 4], [6, 6], [4, 6], [4, 4]]]

    RAND
    POLYGON = [[[12, 0], [30, 0], [30, 15], [18, 15], [18, 30], 
                [0, 30], [0, 15], [12, 15], [12, 0]]]

    SNAKE
    POLYGON = [[[0, 0], [11, 0], [11, 3], [2, 3], [2, 4], [11, 4], 
                [11, 11], [0, 11], [0, 8], [9, 8], [9, 7], [0, 7], [0, 0]]]

    FRAME
    POLYGON = [[[0, 0], [21, 0], [21, 21], [0, 21], [0, 0]], 
                [[6, 6], [15, 6], [15, 15], [6, 15], [6, 6]]]
    '''

    POLYGON = [[[0, 0], [10, 0], [10, 10], [0, 10], [0, 0]]]

    # order follows the edge of the polygon
    MASK_WWR = [0.5]
    MASK_ADIABATIC = [0, 0, 0, 0, 0]

    MESH_SCALE = 1

    LATITUDE = 31.17
    LONGITUDE = 121.43
    UTC = 8

    TIME_STEP = 3600  # in seconds
    TIME_START = datetime.strptime("2022-01-22 00:00:00", "%Y-%m-%d %H:%M:%S")
    TIME_END = datetime.strptime("2022-01-22 23:59:59", "%Y-%m-%d %H:%M:%S")
    TIME_DELTA = timedelta(seconds=TIME_STEP)

    BEAM_RAD = 500

    (snapshots, mask_mesh, mesh_dim) = get_mesh(
        POLYGON, MASK_WWR, MASK_ADIABATIC, MESH_SCALE,
        (LATITUDE, LONGITUDE, UTC),
        (TIME_START, TIME_END, TIME_DELTA),
        BEAM_RAD, 0.6, "", 0)

    # visualize radiation intensity on the floorplan
    maxRadiation = np.max(snapshots)

    fig = plt.figure(figsize=(5, 5), dpi=72)
    lattice_x = np.arange(0, len(snapshots[0][0]) + 1)
    lattice_y = np.arange(0, len(snapshots[0]) + 1)
    cmap = plt.cm.get_cmap('inferno').copy()
    cmap.set_bad(color='w', alpha=1.)
    plt.axis('equal')
    plt.axis('off')  # remove all axes
    plt.xticks([])  # remove ticks on x-axis
    plt.yticks([])  # remove ticks on y-axis
    plt.subplots_adjust(left=0, bottom=0, right=1, top=1)
    images = []
    for snapshot in snapshots:
        for i in range(len(snapshot)):
            for j in range(len(snapshot[0])):
                if mask_mesh[len(snapshot) - i - 1][j] != 1:
                    snapshot[i][j] = np.nan
        snapshot = np.ma.masked_invalid(snapshot)
        # ax = plt.pcolor(lattice_x, lattice_y, snapshot, norm=plt.Normalize(0, maxRadiation))
        ax = plt.pcolormesh(lattice_x, lattice_y, snapshot, cmap=cmap, edgecolors='None',
                            norm=plt.Normalize(0, maxRadiation))
        images.append((ax,))
    img_ani = animation.ArtistAnimation(fig, images, interval=100, repeat_delay=0, blit=True)
    plt.show()

    # pip install pillow
    # https://stackoverflow.com/questions/51512141/how-to-make-matplotlib-saved-gif-looping
    # class LoopingPillowWriter(animation.PillowWriter):
    #     def finish(self):
    #         self._frames[0].save(
    #             self.outfile, save_all=True, append_images=self._frames[1:],
    #             duration=int(500 / self.fps), loop=0)
    # img_ani.save(r"animation.gif", writer=LoopingPillowWriter(fps)

    img_ani.save(r"animation.gif", writer='imagemagick')
