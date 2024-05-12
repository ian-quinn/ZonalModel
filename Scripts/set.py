import math
import datetime

import numpy as np
from mesh import get_mesh
from fake import fake_input
from rad import get_view_factor, get_mrt_factor


# assistants
def retrieve_id(coord, x, y):
    return x * y * coord[2] + y * coord[1] + coord[0]


def retrieve_coords(index, x, y):
    return [index % (x * y) % x, index % (x * y) // x, index // (x * y)]


# Retrieve the adjacent cell in the mesh mask
# -1 for indexation failure
# 0 for void region / 1 for floor plan occupied region
# 2 for void region enclosed / 3 for outer adiabatic boundary
def retrieve_adjacency_mask(index, x, y, mask, port_id):
    pt = retrieve_coords(index, x, y)
    if port_id == 0:
        if pt[0] - 1 >= 0:
            adj = [pt[0] - 1, pt[1]]
            return mask[y - adj[1] - 1, adj[0]]
    if port_id == 1:
        if pt[0] + 1 < x:
            adj = [pt[0] + 1, pt[1]]
            return mask[y - adj[1] - 1, adj[0]]
    if port_id == 2:
        if pt[1] - 1 >= 0:
            adj = [pt[0], pt[1] - 1]
            return mask[y - adj[1] - 1, adj[0]]
    if port_id == 3:
        if pt[1] + 1 < y:
            adj = [pt[0], pt[1] + 1]
            return mask[y - adj[1] - 1, adj[0]]
    # print("WARNING retrieving outside the mesh boundary")
    return -1


def zip_time_table(time, value):
    content = ""
    for i in range(len(time)):
        content += "{}, {:.2f}".format(time[i], value[i])
        if i + 1 != len(time):
            content += "; "
    return content


class Zone(object):
    load = 0.0
    t_0 = 293.15

    def __init__(self, idx, coords, dims):
        self._idx = idx

        self.coords = coords
        self.connect = [0, 0, 0, 0, 0, 0]
        self.port = ["", "", "", "", "", ""]
        self.name = "zone" + str(idx)
        self.dx = dims[0]
        self.dy = dims[1]
        self.dz = dims[2]

    def __repr__(self):
        return "zone{0} ({1}, {2}, {3}) load({4})".format(
            self._idx, self.coords[0], self.coords[1], self.coords[2], self.load)

    def serialize(self):
        return "  Zone {}(IsSource = {}, dx = {}, dy = {}, dz = {}, T_0 = {:.2f});\n" \
            .format(self.name, str(bool(self.load)).lower(), self.dx, self.dy, self.dz, self.t_0)


class Flow(object):
    def __init__(self, idx, direction, idx_a, idx_b, port_a, port_b):
        self._idx = idx

        self.zone_i = idx_a
        self.zone_o = idx_b
        self.port_i = port_a
        self.port_o = port_b
        self.name = "flow" + str(idx)
        self.direction = direction

    def __repr__(self):
        return "flow{0} direction{1} {2}/{3}".format(
            self.zone_i, self.direction, self.zone_i, self.zone_o)

    def serialize(self):
        return "  Flow {}(Direction = {});\n".format(self.name, self.direction)


class Wall(object):
    resistance = 1
    capacity = 10000
    t_0 = 293.15
    is_sunlit = False
    is_adiabatic = False
    is_radiated = False
    is_transparent = False

    def __init__(self, idx, direction, idx_a, port_a, coords, dims):
        area = 1
        centroid = [0, 0, 0]
        if direction == 0:
            area = dims[1] * dims[2]
            centroid = [coords[0], coords[1] + dims[1] / 2, coords[2] + dims[2] / 2]
        if direction == 1:
            area = dims[0] * dims[2]
            centroid = [coords[0] + dims[0] / 2, coords[1], coords[2] + dims[2] / 2]
        if direction == 2:
            area = dims[0] * dims[1]
            centroid = [coords[0] + dims[0] / 2, coords[1] + dims[1] / 2, coords[2]]
        normal = [0, 0, 0]
        normal[port_a // 2] = 1 - port_a % 2 * 2

        self._idx = idx

        self.zone_i = idx_a  # index of the connected zone
        self.port_i = port_a  # one of the zone port [0, 0, 0, 0, 0, 0]
        self.coords = coords  # np.array coordinates of the wall location point (min xyz)
        self.centroid = centroid  # np.array coordinates of the wall centroid
        self.name = "wall" + str(idx)
        self.direction = direction
        self.normal = np.array(normal)
        self.area = area

    def __repr__(self):
        return "wall{0} direction{1} {2}".format(self._idx, self.direction, self.zone_i)

    def serialize(self):
        return (("  Wall {}(IsSource = {}, IsAdiabatic = {}, IsRadiated = {}, Direction = {}, "
                 "X = {:.2f}, Y = {:.2f}, Z = {:.2f}, R = {:.2f}, C = {:.0f}, T_0 = {:.2f});\n")
                .format(self.name, str(self.is_sunlit).lower(),
                        str(self.is_adiabatic).lower(), str(self.is_radiated).lower(),
                        self.direction, self.centroid[0], self.centroid[1], self.centroid[2],
                        self.resistance / self.area, self.capacity * self.area, self.t_0))


# #################################### USER SETUP ######################################

"""
DIAMOND
[[[3, 0], [10, 0], [10, 7], [7, 10], [0, 10], [0, 3], [3, 0]], 
   [[4, 4], [6, 4], [6, 6], [4, 6], [4, 4]]]

RAND
[[[12, 0], [30, 0], [30, 15], [18, 15], [18, 30], [0, 30], [0, 15], [12, 15], [12, 0]]]

SNAKE
[[[0, 0], [33, 0], [33, 9], [6, 9], [6, 12], [33, 12], [33, 33], [0, 33], [0, 24], 
    [27, 24], [27, 21], [0, 21], [0, 0]]]

FRAME
POLYGON = [[[0, 0], [21, 0], [21, 21], [0, 21], [0, 0]], 
    [[6, 6], [15, 6], [15, 15], [6, 15], [6, 6]]]
"""

POLYGON = [[[0, 0], [9, 0], [9, 9], [0, 9], [0, 0]]]
MESH_SCALE = 3
# the wall boundary mask follows each edge in POLYGON[0]
MASK_WWR = [0.8, 0.8, 0.8, 0.8]
MASK_ADIABATIC = [0, 0, 0, 0]
EQUIP_LOC = [(3.5, 4.7)]

# if you need to fake some day slices do set this as 6, 12, 18 or 24
# in this mode, beam radiation is manually set and the diffuse radiation is calculated from it
FAKE_SLICE = 0
BEAM_RAD = 800
# turn on the radiation mode will calculate the view factors of internal surfaces
# then output Radbox.mo for radiation heat transfer. this is time-consuming
RADIATED_MODE = False
PATH_VIEW_FACTOR = ""

MODEL_NAME = "HybriZo"
PATH_EPW = "../Weather/Shanghai.epw"
LATITUDE = 31.17
LONGITUDE = 121.43
UTC = 8

TIME_STEP = 3600  # in seconds, must <= 3600
TIME_INIT = datetime.datetime.strptime("2022-01-01 00:00:00", "%Y-%m-%d %H:%M:%S")
TIME_START = datetime.datetime.strptime("2022-07-22 00:00:00", "%Y-%m-%d %H:%M:%S")
TIME_END = datetime.datetime.strptime("2022-07-23 00:00:00", "%Y-%m-%d %H:%M:%S")
TIME_DELTA = datetime.timedelta(seconds=TIME_STEP)

# presume all equipment has the same load
EQUIP_LOAD = 200
# a 24h schedule for overall internal load
# EQUIP_SCHE = [0, 0, 0, 0, 0, 0, 0, 0.2, 0.5, 1, 1, 1, 0.5, 0.5, 0.8, 1, 1, 1, 1, 0.5, 0.5, 0.5, 0, 0]
EQUIP_SCHE = [1, 1, 1, 1, 1, 1, 1, 0.8, 0.5, 0, 0, 0, 0.5, 0.5, 0.2, 0, 0, 0, 0, 0.5, 0.5, 0.5, 1, 1]

# general heat capacity of wall material (J/K*m2) (cp * rho * thickness)
SHGC = 0.6
C_WIN = 40000
C_WALL = 80000
# general thermal resistance of wall material (m2*K/W) (thickness / conductance + 1 / convection rate)
R_WIN = 0.5
R_WALL = 1

T_DELTA = 0

# ------------------------------ SOME PREPROCESSING ---------------------------
# just for now, include some geometry stuff

pts_outer = []
offset = np.array([MESH_SCALE, MESH_SCALE])
for coords in POLYGON[0]:
    pts_outer.append(np.array(coords) + offset)

# ################################### MESHING #######################################

# snapshots / 3d array in nested lists form representing solar load matrix of each step
# mesh_mask / 2d array representing solid cell in bounding box of target polygon
(snapshots, mesh_mask, cell_dim) = get_mesh(
    POLYGON, MASK_WWR, MASK_ADIABATIC, MESH_SCALE,
    (LATITUDE, LONGITUDE, UTC),
    (TIME_START, TIME_END, TIME_DELTA),
    BEAM_RAD, SHGC, PATH_EPW, FAKE_SLICE)

# passed down geometry information
dim_x = cell_dim[0]
dim_y = cell_dim[1]
dim_z = 1.0  # default
dims = np.array([dim_x, dim_y, dim_z])

tick_x = np.size(mesh_mask, 1)
tick_y = np.size(mesh_mask, 0)
tick_z = 3  # default

# simulation time for Modelica settings
time_simulation = (TIME_END - TIME_START).total_seconds()
if FAKE_SLICE:
    time_simulation = 3600 * 24 * FAKE_SLICE

# flip the snapshots to get the beam load
# axis-0: time dimension, the number of snapshots
# axis-1: space dimension, each row numbered from south to north ↑
# axis-2: space dimension, each column from west to east →
# if you flip the snapshots[i] upsidedown, you get the solar distribution on the floorplan at time i
solar_series = (snapshots
                .reshape(snapshots.shape[0], snapshots.shape[1] * snapshots.shape[2])
                .transpose())
solar_lumped_series = np.sum(snapshots, axis=(1, 2))

# internal heat source, directly to zone air volume, universal time schedule
# zone.load * series_equip_sche is the series of zone load
num_increments = int(math.ceil(time_simulation / TIME_STEP))
series_equip_sche = np.zeros(num_increments)
series_equip_zero = np.zeros(num_increments)
for i in range(num_increments):
    sche_span = int(TIME_STEP * i / 3600 % len(EQUIP_SCHE))
    series_equip_sche[i] = EQUIP_SCHE[sche_span]
# for i in range(np.size(equip_series, 0)):
#     for j in range(np.size(equip_series, 1)):
#         sche_span = int(TIME_STEP * j / 3600 % len(EQUIP_SCHE))
#         equip_series[i, j] = EQUIP_LOAD * EQUIP_SCHE[sche_span]

# grab timestamps
# the timestamps returned by get_mesh() should correspond to the following series_time
# make sure all functions follow the TIME_STEP setting
# the question is, time step smaller than 3600s worth consideration?
timestamps = [TIME_STEP * tick for tick in range(len(snapshots))]

# series of time marks the timestamps of each simulation step
# series of temp marks the temperature increment of each simulation step
# series of solar diffusion marks the solar evenly loaded on each tile
series_time = []
series_temp = []
series_solar_diff = []

if FAKE_SLICE:
    series_temp, _ = fake_input(PATH_EPW, FAKE_SLICE, TIME_STEP)
    series_time = [TIME_STEP * i for i in range(len(series_temp))]
    series_temp = series_temp + T_DELTA
    room_area = tick_x * tick_y * MESH_SCALE * MESH_SCALE
    for solar_lumped in solar_lumped_series:
        solar_per_area = solar_lumped / room_area
        # not valid. just one plausible relation between beam and diffusion
        if solar_per_area > 0:
            series_solar_diff.append(math.sqrt(solar_per_area) * 10 * SHGC)
        else:
            series_solar_diff.append(0.0)

# grab temperature series from epw file
# each row for series data in epw has 35 columns
# epw has 8760 rows for each hour in total
# No.6 is the outdoor dry-bulb temperature
# No.14 is the solar diffusion radiation
else:
    hour_tick = 0
    f = open(PATH_EPW, mode="r")
    for line in f.readlines():
        data_list = line.split(',')
        if len(data_list) == 35:
            hour_tick += 1
            if hour_tick * 3600 > (TIME_START - TIME_INIT).total_seconds() and \
                    (hour_tick - 1) * 3600 < (TIME_END - TIME_INIT).total_seconds():
                series_temp.append(round(float(data_list[6]) + 273.15 + T_DELTA, 2))
                series_solar_diff.append(int(data_list[14]) * SHGC)
                if len(series_temp) == 0:
                    series_time.append(0.0)
                else:
                    series_time.append((hour_tick - 1) * 3600 - (TIME_START - TIME_INIT).total_seconds())

temp_init = series_temp[0]

if len(timestamps) != len(series_time):
    print("WARNING! Timestamps not match {} - {}".format(len(timestamps), len(series_time)))

# ################################# MODULE GENERATION #################################

"""
   x→               None?            index                 coords
 y 0 0 0 0 0        N N N N N        20 21 22 23 24        * *
 ↓ 0 1 0 0 0        N Z N N N        15 16 17 18 19        * (1,3) *
   0 1 1 0 0   ->   N Z Z N N   ->   10 11 12 13 14   ->   * (1,2) (2,2) *
   0 1 1 1 0        N Z Z Z N        05 06 07 08 09        * (1,1) (2,1) (3,1) *
   0 0 0 0 0        N N N N N        00 01 02 03 04        * *     *     *     * *
"""
# zone will cover all over the mesh regardless of inside/outside the floor plan
# the index of each zone must be projected to a coordinate one and only
zones = []
for k in range(tick_z):
    for j in range(tick_y):
        for i in range(tick_x):
            if mesh_mask[tick_y - j - 1][i] == 1:
                new_zone = Zone(len(zones), [i, j, k], [dim_x, dim_y, dim_z])
                new_zone.t_0 = temp_init

                # check if this zone contains equipment
                equip_load = 0
                for loc in EQUIP_LOC:
                    if (i - 1) * dim_x < loc[0] < i * dim_x and \
                            (j - 1) * dim_y < loc[1] < j * dim_y and \
                            k == 0:
                        equip_load += EQUIP_LOAD
                if equip_load > 0:
                    new_zone.load = equip_load

                zones.append(new_zone)
            else:
                zones.append(None)

flows = []
for i in range(len(zones)):
    if zones[i] is None:
        continue
    coords = zones[i].coords
    if coords[0] < tick_x - 1:
        # print("looping " + str(coords[0]) + " | " + str(coords[1]))
        if mesh_mask[tick_y - coords[1] - 1][coords[0] + 1] == 1:
            targetId = retrieve_id([coords[0] + 1, coords[1], coords[2]], tick_x, tick_y)
            flows.append(Flow(len(flows), 0, targetId, i, 0, 1))
            zones[i].connect[1] = 1
            zones[i].port[1] = flows[-1].name
            zones[targetId].connect[0] = 1
            zones[targetId].port[0] = zones[i].port[1]
    if coords[1] < tick_y - 1:
        if mesh_mask[tick_y - coords[1] - 2][coords[0]] == 1:
            flows.append(Flow(len(flows), 1,
                              retrieve_id([coords[0], coords[1] + 1, coords[2]], tick_x, tick_y), i, 2, 3))
            zones[i].connect[3] = 1
            zones[i].port[3] = flows[-1].name
            zones[retrieve_id([coords[0], coords[1] + 1, coords[2]], tick_x, tick_y)].connect[2] = 1
            zones[retrieve_id([coords[0], coords[1] + 1, coords[2]], tick_x, tick_y)].port[2] = zones[i].port[3]
    if coords[2] < tick_z - 1:
        flows.append(Flow(len(flows), 2,
                          retrieve_id([coords[0], coords[1], coords[2] + 1], tick_x, tick_y), i, 4, 5))
        zones[i].connect[5] = 1
        zones[i].port[5] = flows[-1].name
        zones[retrieve_id([coords[0], coords[1], coords[2] + 1], tick_x, tick_y)].connect[4] = 1
        zones[retrieve_id([coords[0], coords[1], coords[2] + 1], tick_x, tick_y)].port[4] = zones[i].port[5]

'''
for zone in zones:
    if zone != None:
        chain = ""
        for i in range(len(zone._connect)):
            chain += ", " + str(zone._connect[i])
        print(chain)
'''

walls = []
for i in range(len(zones)):
    if zones[i] is None:
        continue
    connectors = zones[i].connect
    for j in range(len(connectors)):
        # any exposed connector (exposed) will be assigned with a wall
        if not connectors[j]:
            # stretch zone index to get the base point of the zone
            loc_pt = np.array(zones[i].coords) * dims
            # if the wall has negative normal vector
            # plus the dimension size accordingly to get its location point
            loc_pt[j // 2] += (j % 2) * dims[j // 2]
            new_wall = Wall(len(walls), j // 2, i, j, loc_pt, [dim_x, dim_y, dim_z])
            walls.append(new_wall)
            zones[i].port[j] = new_wall.name

# determine the wall types
for wall in walls:
    if wall.direction == 2:
        wall.is_adiabatic = True
        if wall.zone_i < tick_x * tick_y:
            wall.is_sunlit = True
    elif retrieve_adjacency_mask(wall.zone_i, tick_x, tick_y, mesh_mask, wall.port_i) != 0:
        wall.is_adiabatic = True

    if RADIATED_MODE:
        wall.is_radiated = True

    wall.t_0 = temp_init

    print("initial value: R-{} C-{}".format(wall.resistance, wall.capacity))
    loc_pt = np.array([wall.coords[0], wall.coords[1]])
    vec_norm = np.array([wall.normal[0], wall.normal[1]])
    # skim out horizontal surface
    if np.linalg.norm(vec_norm) == 0:
        wall.resistance = R_WALL
        wall.capacity = C_WALL
        continue
    for i in range(len(pts_outer) - 1):
        if i + 1 <= len(MASK_WWR):
            if MASK_WWR[i] > 0:
                vec_edge = pts_outer[i + 1] - pts_outer[i]
                check_vertical = np.dot(vec_edge, vec_norm)
                check_parallel = np.cross(pts_outer[i + 1] - loc_pt, pts_outer[i] - loc_pt)
                if check_vertical == 0 and check_parallel == 0:
                    wall.is_transparent = True
                    wall.resistance = R_WIN
                    wall.capacity = C_WIN
                    break
                else:
                    wall.resistance = R_WALL
                    wall.capacity = C_WALL

# ################################# VIEW MATRIX #####################################


discretization_num = 5
mrt_measure_height = 1.2

wall_normals = []
wall_basepts = []
for wall in walls:
    # wall.coords type is np.array
    wall_basepts.append(wall.coords)
    # wall.normal type is np.array
    wall_normals.append(wall.normal)

if RADIATED_MODE:
    if PATH_VIEW_FACTOR:
        view_matrix = np.genfromtxt(PATH_VIEW_FACTOR, delimiter=',')
    else:
        view_matrix = get_view_factor(wall_basepts, wall_normals, dims, discretization_num)

mrt_matrix = get_mrt_factor(wall_basepts, wall_normals, dims, mesh_mask, discretization_num, mrt_measure_height)

# ################################# SERIALIZATION ####################################


# cache the mask matrix to .csv for post-processing
np.savetxt("mask.csv", mesh_mask, fmt='%d', delimiter=',')

# cache the internal view factor matrix for checking
if not PATH_VIEW_FACTOR and RADIATED_MODE:
    np.savetxt("view.csv", view_matrix, fmt='%2f', delimiter=',')

# cache the MRT matrix to .csv for post-processing
np.savetxt("mrt.csv", mrt_matrix, fmt='%2f', delimiter=',')

# ############################## FILE RADIATION MODEL #################################


if RADIATED_MODE:
    # not sure about the solar distribution right now, assuming the solar absorption = 1
    # so there's no redistribution of solar load on inner surfaces
    alpha_solar = 0.8
    rho_solar = 0.2

    alpha_infra = 0.8
    rho_infra = 0.2
    tao_infra = 0.25

    fo = open("../VEPZO/Radbox.mo", "w")
    fo.write("within VEPZO;\n")
    fo.write("model Radbox\n")
    fo.write("  HeatPort port_s;\n")  # for diffuse solar load
    for i in range(len(walls)):
        fo.write("  HeatPort port_" + str(i) + ";\n")
        fo.write("  SI.HeatFlowRate Qin_" + str(i) + ";\n")
        fo.write("  SI.HeatFlowRate Qout_" + str(i) + ";\n")
    fo.write("equation\n")
    fo.write("  port_s.T = 293.15;\n".format(i))
    for i in range(len(walls)):
        Sigma = ""
        for j in range(len(walls)):
            if i != j:
                Sigma += "Qout_{0} * {1} + ".format(j, round(view_matrix[i][j], 4))
        Sigma = Sigma[:len(Sigma) - 3]
        if walls[i].is_transparent:
            fo.write("  Qout_{0} = {1} * Modelica.Constants.sigma * port_{0}.T ^ 4 + \
                port_s.Q_flow + {2} * Qin_{0};\n".format(i, alpha_infra, tao_infra))

        else:
            fo.write("  Qout_{0} = {1} * Modelica.Constants.sigma * port_{0}.T ^ 4 + {2} * Qin_{0};\n"
                     .format(i, alpha_infra, rho_infra))
        fo.write("  Qin_{0} = {1};\n".format(i, Sigma))
        fo.write("  port_{0}.Q_flow = Qin_{0} - Qout_{0};\n".format(i))
    fo.write("end Radbox;\n")
    fo.close()

# ################################# FILE ROOM MODEL ##################################


fo = open("../VEPZO/Samples/" + MODEL_NAME + ".mo", "w")
fo.write("within VEPZO.Samples;\n")
fo.write("model " + MODEL_NAME + "\n")
# fo.write("  extends VEPZO.VariableSelection;\n")

# # MODULE
fo.write("  Modelica.Blocks.Sources.TimeTable ambient(table = [{0}]);\n"
         .format(zip_time_table(series_time, series_temp)))
# use global equipment load schedule
# fo.write("  Modelica.Blocks.Sources.TimeTable equipLoadSchedule(table = [{0}]);\n"
#          .format(zip_time_table(timestamps, series_equip_sche * EQUIP_LOAD)))
fo.write("  Modelica.Blocks.Sources.TimeTable equipZeroSchedule(table = [{0}]);\n"
         .format(zip_time_table(timestamps, series_equip_zero)))
for zone in zones:
    if zone is not None:
        fo.write(zone.serialize())
        fo.write("  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow loadEquip{0};\n".format(zone._idx))
        if zone.load:
            fo.write("  Modelica.Blocks.Sources.TimeTable equipLoadSchedule{0}(table = [{1}]);\n"
                     .format(zone._idx, zip_time_table(timestamps, series_equip_sche * zone.load)))
for flow in flows:
    fo.write(flow.serialize())
for wall in walls:
    print("{} - {} - C:{}".format(wall.name, wall.is_transparent, wall.capacity))
    # print("Wall-" + str(walls.index(wall)) + " assigned boundary condition ", end='')
    # type 1, by default the wall is adiabatic without solar load
    fo.write(wall.serialize())

    # type 2, on the floor, adiabatic with solar stress
    if wall.is_sunlit:
        fo.write("  Modelica.Blocks.Sources.TimeTable solarBeam{0}(table = [{1}]);\n"
                 .format(wall._idx, zip_time_table(timestamps, solar_series[wall.zone_i])))
        fo.write("  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow loadBeam{0};\n".format(wall._idx))
    else:
        fo.write("  Adiabatic adiabaticS{0};\n".format(wall._idx))
    # type 3, perimeter, assigned prescribed environment temperature
    if wall.is_adiabatic:
        fo.write("  Adiabatic adiabaticT{0};\n".format(wall._idx))
    else:
        fo.write("  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature temp{};\n".format(wall._idx))
        fo.write("  Modelica.Thermal.HeatTransfer.Components.ThermalResistor wallR{}(R = {:.2f});\n"
                 .format(wall._idx, 1 / 7.6 / wall.area))
        fo.write("  Bound bound{0};\n".format(wall._idx))

# include the radiation box model
if RADIATED_MODE:
    fo.write("  Radbox radbox;\n")
    fo.write("  Modelica.Blocks.Sources.TimeTable solarDiff(table = [{0}]);\n"
             .format(zip_time_table(series_time, series_solar_diff)))
    fo.write("  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow loadDiff;\n")

# # EQUATION
fo.write("equation\n")

for flow in flows:
    axis = "x"
    if flow.direction == 1: axis = "y"
    if flow.direction == 2: axis = "z"
    fo.write("  connect({0}.port_a, {1}.port_{2}1);\n".format(flow.name, "zone" + str(flow.zone_i), axis))
    fo.write("  connect({0}.port_b, {1}.port_{2}2);\n".format(flow.name, "zone" + str(flow.zone_o), axis))
for wall in walls:
    axis = "x"
    if wall.direction == 1: axis = "y"
    if wall.direction == 2: axis = "z"
    fo.write("  connect({0}.port_a, {1}.port_{2}{3});\n".format(wall.name, "zone" + str(wall.zone_i),
                                                                axis, str(wall.port_i % 2 + 1)))
for flow in flows:
    fo.write("  connect({0}.i[1], {1}.o);\n".format(flow.name, "zone" + str(flow.zone_o)))
    fo.write("  connect({0}.i[2], {1}.o);\n".format(flow.name, "zone" + str(flow.zone_i)))

for zone in zones:
    if zone is not None:
        fo.write("  connect({0}.port_s, loadEquip{1}.port);\n".format(zone.name, zone._idx))
        if zone.load:
            fo.write("  connect(equipLoadSchedule{0}.y, loadEquip{0}.Q_flow);\n".format(zone._idx))
        else:
            fo.write("  connect(equipZeroSchedule.y, loadEquip{0}.Q_flow);\n".format(zone._idx))
        for i in range(6):
            fo.write("  connect({0}.i[{1}], {2}.o);\n".format(zone.name, i + 1, zone.port[i]))

for wall in walls:
    fo.write("  connect({0}.i, {1}.o);\n".format(wall.name, "zone" + str(wall.zone_i)))
    adjMaskValue = retrieve_adjacency_mask(wall.zone_i, tick_x, tick_y, mesh_mask, wall.port_i)
    if wall.is_sunlit:
        fo.write("  connect(wall{0}.port_s, loadBeam{0}.port);\n".format(wall._idx))
        fo.write("  connect(solarBeam{0}.y, loadBeam{0}.Q_flow);\n".format(wall._idx))
    else:
        fo.write("  connect(wall{0}.port_s, adiabaticS{0}.port);\n".format(wall._idx))
    if wall.direction != 2 and adjMaskValue == 0:
        fo.write("  connect(temp{0}.T, ambient.y);\n".format(wall._idx))
        fo.write("  connect(temp{0}.port, wallR{0}.port_b);\n".format(wall._idx))
        fo.write("  connect(wallR{0}.port_a, wall{0}.port_t);\n".format(wall._idx))
        fo.write("  connect(temp{0}.port, bound{0}.port_b);\n".format(wall._idx))
        fo.write("  connect(bound{0}.port_a, wall{0}.port_t);\n".format(wall._idx))
        # fo.write("  connect(wallC{0}.port, wallR{0}.port_b);\n".format(wall._idx))
    else:
        fo.write("  connect(adiabaticT{0}.port, wall{0}.port_t);\n".format(wall._idx))

    # to radiation box
    if RADIATED_MODE:
        fo.write("  connect(wall{0}.port_r, radbox.port_{0});\n".format(wall._idx))

if RADIATED_MODE:
    fo.write("  connect(radbox.port_s, loadDiff.port);\n")
    fo.write("  connect(solarDiff.y, loadDiff.Q_flow);\n")

# fo.write("  MatrixConverter();\n")
# experiment configuration only for Wolfram
fo.write("annotation(experiment(StopTime = {0}, __Wolfram_NumberOfIntervals = -1));\n".format(time_simulation))

fo.write("end " + MODEL_NAME + ";\n")
fo.close()
