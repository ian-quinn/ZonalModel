import math
import time
import datetime

import numpy as np
from mesh import GetMesh
from fake import FakeInput
from rad import GetViewFactor

# assistants
def RetrieveId(coord, x, y):
    return x * y * coord[2] + y * coord[1] + coord[0]

def RetrieveCoord(id, x, y):
    return [id % (x * y) % x, id % (x * y) // x, id // (x * y)]

# Retrieve the adjacent cell in the mesh mask
# -1 for indexation failure
# 0 for void region / 1 for floorplan occupied region
# 2 for void region enclosed / 3 for outer adiabatic boundary
def RetrieveAdjacencyMask(index, x, y, mask, portId):
    pt = RetrieveCoord(index, x, y)
    if portId == 0:
        if pt[0] - 1 >= 0:
            adj = [pt[0] - 1, pt[1]]
            return mask[y - adj[1] - 1, adj[0]]
    if portId == 1:
        if pt[0] + 1 < x:
            adj = [pt[0] + 1, pt[1]]
            return mask[y - adj[1] - 1, adj[0]]
    if portId == 2:
        if pt[1] - 1 >= 0:
            adj = [pt[0], pt[1] - 1]
            return mask[y - adj[1] - 1, adj[0]]
    if portId == 3:
        if pt[1] + 1 < y:
            adj = [pt[0], pt[1] + 1]
            return mask[y - adj[1] - 1, adj[0]]
    # print("WARNING retrieving outside the mesh boundary")
    return -1

def ZipTimeTable(time, value):
    content = ""
    for i in range(len(time)):
        content += "{}, {:.2f}".format(time[i], value[i])
        if i + 1 != len(time):
            content += "; "
    return content

# module class

class Zone(object):
    is_loaded = False
    t_0 = 293.15
    def __init__(self, idx, coords, dims):
        self._idx = idx
        self._coords = coords
        self._connect = [0, 0, 0, 0, 0, 0]
        self._port = ["", "", "", "", "", ""]

        self.name = "zone" + str(idx)
        self.dx = dims[0]
        self.dy = dims[1]
        self.dz = dims[2]

    def __repr__(self):
        print("zone{0} ({1}, {2}, {3})".format(
            self.idx, self.coords[0], self.coords[1], self.coords[2]))

    def Serialize(self):
        return "  Zone {}(IsSource = {}, dx = {}, dy = {}, dz = {}, T_0 = {:.2f});\n"\
            .format(self.name, str(self.is_loaded).lower(), self.dx, self.dy, self.dz, self.t_0)

class Flow(object):
    def __init__(self, idx, direction, idxA, idxB, portA, portB):
        self._idx_zone_i = idxA
        self._idx_zone_o = idxB
        self._idx_port_i = portA
        self._idx_port_o = portB

        self.name = "flow" + str(idx)
        self.direction = direction
        
    def __repr__(self):
        print("flow{0} direction{1} {2}/{3}".format(self.idxIn, self.direction, self.idxIn, self.idxOut))

    def Serialize(self):
        return "  Flow {}(Direction = {});\n".format(self.name, self.direction)

class Wall(object):
    resistance = 1
    capacity = 10000
    t_0 = 293.15
    is_solarloaded = False
    is_adiabatic = False
    is_radiated = False
    is_transparent = False
    def __init__(self, idx, direction, idxA, portA, coords, dims):
        area = 1
        if direction == 0: area = dims[1] * dims[2]
        if direction == 1: area = dims[0] * dims[2]
        if direction == 2: area = dims[0] * dims[1]
        normal = [0, 0, 0]
        normal[portA // 2] = 1 - portA % 2 * 2

        self._idx = idx
        self._idx_zone_i = idxA        # index of the connected zone
        self._idx_port_i = portA       # one of the zone port [0, 0, 0, 0, 0, 0]
        self._coords = coords           # np.array coords of the wall basepoint (min xyz)

        self.name = "wall" + str(idx)
        self.direction = direction
        self.normal = np.array(normal)
        self.area = area
        
    def __repr__(self):
        print("wall{0} direction{1} {2}".format(self.idx, self.direction, self.idxIn))

    def Serialize(self):
        return "  Wall {}(IsSource = {}, IsAdiabatic = {}, IsRadiated = {}, Direction = {}, \
            R = {:.2f}, C = {:.0f}, T_0 = {:.2f});\n"\
            .format(self.name, str(self.is_solarloaded).lower(), 
                str(self.is_adiabatic).lower(), str(self.is_radiated).lower(), 
                self.direction, self.resistance / self.area, self.capacity * self.area, self.t_0)
        
#################### USER SETUP ####################
'''
DIAMOND
[[[3, 0], [10, 0], [10, 7], [7, 10], [0, 10], [0, 3], [3, 0]], 
   [[4, 4], [6, 4], [6, 6], [4, 6], [4, 4]]]

RAND
[[[12, 0], [30, 0], [30, 15], [18, 15], [18, 30], [0, 30], [0, 15], [12, 15], [12, 0]]]

SNAKE
[[[0, 0], [33, 0], [33, 9], [6, 9], [6, 12], [33, 12], [33, 33], [0, 33], [0, 24], \
   [27, 24], [27, 21], [0, 21], [0, 0]]]

FRAME
POLYGON = [[[0, 0], [21, 0], [21, 21], [0, 21], [0, 0]], \
    [[6, 6], [15, 6], [15, 15], [6, 15], [6, 6]]]
'''
POLYGON = [[[0, 0], [9, 0], [9, 9], [0, 9], [0, 0]]]

# if you need to fake some day slices do set this as 6, 12, 18 or 24
FAKE_SLICE = 0
RADIATED_MODE = True

MODEL_NAME = "test"
PATH_VIEW_FACTOR = "view_9_3.csv"
PATH_EPW = "../Weather/Shanghai.epw"

MASK_WWR = [0.8]
MASK_ADIABATIC = [0, 0, 0, 0]

MESH_SCALE = 3

LATITUDE = 31.17
LONGITUDE = 121.43
UTC = 8

TIME_STEP = 600                 # in seconds
TIME_INIT = datetime.datetime.strptime("2022-01-01 00:00:00", "%Y-%m-%d %H:%M:%S")
TIME_START = datetime.datetime.strptime("2022-05-15 00:00:00", "%Y-%m-%d %H:%M:%S")
TIME_END = datetime.datetime.strptime("2022-06-30 23:59:59", "%Y-%m-%d %H:%M:%S")
TIME_DELTA = datetime.timedelta(seconds=TIME_STEP)

BEAM_RAD = 800
SHGC = 0.3

# general heat capacity of wall material (J/K*m2) (cp * rho * thickness)
C_WIN = 10000      
C_WALL = 20000
# general thermal resistance of wall material (m2*K/W) (thickness / conductance + 1 / convection rate)
R_WIN = 0.5
R_WALL = 1

T_DELTA = 0

##################### SOME PREPROCESSING ###################
# just for now, include several geometry stuff
pts_outer = []
offset = np.array([MESH_SCALE, MESH_SCALE])
for coords in POLYGON[0]:
    pts_outer.append(np.array(coords) + offset)


##################### MESHING WORK ###################

# good to go for the meshing work
# snapshots / 3d array in nested lists form representing solar load matrix of each step
# mesh_mask / 2d array representing solid cell in bounding box of target polygon
(snapshots, mesh_mask, cell_dim) = GetMesh(
    POLYGON, MASK_WWR, MASK_ADIABATIC, MESH_SCALE, 
    (LATITUDE, LONGITUDE, UTC), 
    (TIME_START, TIME_END, TIME_DELTA), 
    BEAM_RAD, SHGC, PATH_EPW, FAKE_SLICE)

# passed down geometry information
dimx = cell_dim[0]
dimy = cell_dim[1]
dimz = 1                        # default as always
dims = np.array([dimx, dimy, dimz])

tickx = np.size(mesh_mask, 1)
ticky = np.size(mesh_mask, 0)
tickz = 3                       # default as always

# flip the snapshots to get the beam load
solardistributions = []
for snapshot in snapshots:
    solardistribution = [item for row in snapshot for item in row]
    solardistributions.append(solardistribution)
solarseries = np.array(solardistributions).transpose()
solarseries_lumped = sum(solarseries)

time_simulation = (TIME_END - TIME_START).total_seconds()
if FAKE_SLICE:
    time_simulation = 3600 * 24 * FAKE_SLICE

# grab timestamps
# the timestamps returned by GetMesh() should correspond to the following series_time
# make sure all functions follow the TIME_STEP setting
# the question is, do we all need a small time step?
timestamps = []          # time stamps in second
for i in range(len(snapshots)):
    timestamps.append(i * TIME_STEP)


# grab temperature series from epw file
if FAKE_SLICE:
    series_temp, actual_temp = FakeInput(PATH_EPW, FAKE_SLICE, TIME_STEP)
    series_time = [TIME_STEP * i for i in range(len(series_temp))]
    series_temp = series_temp + T_DELTA
    temp_init = series_temp[0]
    series_diff = []
    roomarea = tickx * ticky * MESH_SCALE * MESH_SCALE
    for i in range(len(solarseries_lumped)):
        solarperarea = solarseries_lumped[i] / roomarea
        if solarperarea > 0:
            series_diff.append(math.sqrt(solarperarea) * 10 * SHGC)
        else:
            series_diff.append(0)
    # print(series_diff)
else:
    series_temp = []
    series_time = []
    series_diff = []
    timeticks = 0
    f = open(PATH_EPW, mode="r")
    for line in f.readlines():
        datalist = line.split(',')
        if len(datalist) == 35:
            timeticks = timeticks + 1
            if timeticks * 3600 > (TIME_START - TIME_INIT).total_seconds() and \
                (timeticks - 1) * 3600 < (TIME_END - TIME_INIT).total_seconds():
                series_temp.append(round(float(datalist[6]) + 273.15 + T_DELTA, 2))
                series_diff.append(int(datalist[14]) * SHGC)
                if len(series_temp) == 0:
                    series_time.append(0.0)
                else:
                    series_time.append((timeticks - 1) * 3600 - (TIME_START - TIME_INIT).total_seconds())
    temp_init = series_temp[0]

if len(timestamps) != len(series_time):
    print("########################################################")
    print("WARNING! Time stamps not match {} - {}".format(len(timestamps), len(series_time)))


################# MODULE GENERATION ##################

# zone will cover all over the mesh regardless of inside/outside the floorplan
# the index of each zone must projected to coordinates one and only
zonelist = [] # list of module Zone
for k in range(tickz):
    for j in range(ticky):
        for i in range(tickx):
            if mesh_mask[ticky - j - 1][i] == 1:
                new_zone = Zone(len(zonelist), [i, j, k], [dimx, dimy, dimz])
                new_zone.t_0 = temp_init
                zonelist.append(new_zone)
                # print(new_zone._coords)
            else:
                zonelist.append(None)
# zonecheck = ""
# for zone in zonelist:
#     if zone == None: 
#         zonecheck += " 0"
#     else:
#         zonecheck += " 1"
# print(zonecheck)

flowlist = [] # list of module Flow
for i in range(len(zonelist)):
    if zonelist[i] == None:
        continue
    coords = zonelist[i]._coords
    if coords[0] < tickx - 1:
        # print("looping " + str(coords[0]) + " | " + str(coords[1]))
        if mesh_mask[ticky - coords[1] - 1][coords[0] + 1] == 1:
            targetId = RetrieveId([coords[0] + 1, coords[1], coords[2]], tickx, ticky)
            flowlist.append(Flow(len(flowlist), 0, targetId, i, 0, 1))
            zonelist[i]._connect[1] = 1
            zonelist[i]._port[1] = flowlist[-1].name
            zonelist[targetId]._connect[0] = 1
            zonelist[targetId]._port[0] = zonelist[i]._port[1]
    if coords[1] < ticky - 1:
        if mesh_mask[ticky - coords[1] - 2][coords[0]] == 1:
            flowlist.append(Flow(len(flowlist), 1, 
                RetrieveId([coords[0], coords[1] + 1, coords[2]], tickx, ticky), i, 2, 3))
            zonelist[i]._connect[3] = 1
            zonelist[i]._port[3] = flowlist[-1].name
            zonelist[RetrieveId([coords[0], coords[1] + 1, coords[2]], tickx, ticky)]._connect[2] = 1
            zonelist[RetrieveId([coords[0], coords[1] + 1, coords[2]], tickx, ticky)]._port[2] = zonelist[i]._port[3]
    if coords[2] < tickz - 1:
        flowlist.append(Flow(len(flowlist), 2, 
            RetrieveId([coords[0], coords[1], coords[2] + 1], tickx, ticky), i, 4, 5))
        zonelist[i]._connect[5] = 1
        zonelist[i]._port[5] = flowlist[-1].name
        zonelist[RetrieveId([coords[0], coords[1], coords[2] + 1], tickx, ticky)]._connect[4] = 1
        zonelist[RetrieveId([coords[0], coords[1], coords[2] + 1], tickx, ticky)]._port[4] = zonelist[i]._port[5]

'''
for zone in zonelist:
    if zone != None:
        chain = ""
        for i in range(len(zone._connect)):
            chain += ", " + str(zone._connect[i])
        print(chain)
'''

walllist = [] # list of module Wall
for i in range(len(zonelist)):
    if zonelist[i] == None:
        continue
    connectors = zonelist[i]._connect
    for j in range(len(connectors)):
        # any exposed connector (exposed) will be assigned with a wall
        if not connectors[j]:
            # stretch zone index to get the base point of the zone
            basepoint = np.array(zonelist[i]._coords) * dims
            # if the wall has negative normal vector
            # plus the dimension size accordingly to get its location point
            basepoint[j // 2] += (j % 2) * dims[j // 2]
            new_wall = Wall(len(walllist), j // 2, i, j, basepoint, [dimx, dimy, dimz])
            walllist.append(new_wall)
            zonelist[i]._port[j] = new_wall.name

# determine the wall types
for wall in walllist:
    if wall.direction == 2:
        wall.is_adiabatic = True
        if wall._idx_zone_i < tickx * ticky:
            wall.is_solarloaded = True
    elif RetrieveAdjacencyMask(wall._idx_zone_i, tickx, ticky, mesh_mask, wall._idx_port_i) != 0:
        wall.is_adiabatic = True

    if RADIATED_MODE:
        wall.is_radiated = True

    wall.t_0 = temp_init

    print("initial value: R-{} C-{}".format(wall.resistance, wall.capacity))
    basept = np.array([wall._coords[0], wall._coords[1]])
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
                check_verticle = np.dot(vec_edge, vec_norm)
                check_parallel = np.cross(pts_outer[i + 1] - basept, pts_outer[i] - basept)
                if check_verticle == 0 and check_parallel == 0:
                    wall.is_transparent = True
                    wall.resistance = R_WIN
                    wall.capacity = C_WIN
                    break
                else:
                    wall.resistance = R_WALL
                    wall.capacity = C_WALL

#################### VIEW MATRIX #####################

if RADIATED_MODE:
    wallnormals = []
    wallbasepts = []
    for wall in walllist:
        wallbasepts.append(wall._coords)
        wallnormals.append(wall.normal)

    if PATH_VIEW_FACTOR:
        view_matrix = np.genfromtxt(PATH_VIEW_FACTOR, delimiter=',')
    else:
        view_matrix = GetViewFactor(wallbasepts, wallnormals, dims, 0.5)


######################################################################################
######################################################################################
################################### SERIALIZATION ####################################


# dump the mask matrix to csv for post-processing
np.savetxt("mask.csv", mesh_mask, fmt='%d', delimiter=',')
if not PATH_VIEW_FACTOR and RADIATED_MODE:
    np.savetxt("view.csv", view_matrix, fmt='%2f', delimiter=',')


##################### FILE RADIATION MODEL ###################
if RADIATED_MODE:
    # not sure about the solar distribution right now, assuming the solar absorption = 1
    # so there's no redistribution of solar load on inner surfaces
    alpha_solar = 0.8
    rho_solar = 0.2

    alpha_infra = 0.8
    rho_infra = 0.2
    tao_infra = 0.25

    fo = open("../Radbox.mo", "w")
    fo.write("within VEPZO;\n")
    fo.write("model Radbox\n")
    fo.write("  HeatPort port_s;\n") # for diffuse solar load
    for i in range(len(walllist)):
        fo.write("  HeatPort port_" + str(i) + ";\n")
        fo.write("  SI.HeatFlowRate Qin_" + str(i) + ";\n")
        fo.write("  SI.HeatFlowRate Qout_" + str(i) + ";\n")
    fo.write("equation\n")
    fo.write("  port_s.T = 293.15;\n".format(i));
    for i in range(len(walllist)):
        Sigma = ""
        for j in range(len(walllist)):
            if i != j:
                Sigma += "Qout_{0} * {1} + ".format(j, round(view_matrix[i][j], 4))
        Sigma = Sigma[:len(Sigma) - 3]
        if walllist[i].is_transparent:
            fo.write("  Qout_{0} = {1} * Modelica.Constants.sigma * port_{0}.T ^ 4 + \
                port_s.Q_flow + {2} * Qin_{0};\n".format(i, alpha_infra, tao_infra))

        else:
            fo.write("  Qout_{0} = {1} * Modelica.Constants.sigma * port_{0}.T ^ 4 + {2} * Qin_{0};\n"
                .format(i, alpha_infra, rho_infra))
        fo.write("  Qin_{0} = {1};\n".format(i, Sigma))
        fo.write("  port_{0}.Q_flow = Qin_{0} - Qout_{0};\n".format(i))
    fo.write("end Radbox;\n")
    fo.close()


##################### FILE ROOM MODEL ###################
fo = open("../Samples/" + MODEL_NAME + ".mo", "w")
fo.write("within VEPZO.Samples;\n")
fo.write("model " + MODEL_NAME + "\n")

### MODULE
fo.write("  Modelica.Blocks.Sources.TimeTable ambient(table = [{0}]);\n"
    .format(ZipTimeTable(series_time, series_temp)))
for zone in zonelist:
    if zone != None:
        fo.write(zone.Serialize())
for flow in flowlist:
    fo.write(flow.Serialize())
for wall in walllist:
    print("{} - {} - C:{}".format(wall.name, wall.is_transparent, wall.capacity))
    # print("Wall-" + str(walllist.index(wall)) + " assigned boundary condition ", end='')
    # type 1, by default the wall is adiabatic without solarload
    fo.write(wall.Serialize())
    
    # type 2, on the floor, adiabatic with solar stress
    if wall.is_solarloaded:
        fo.write("  Modelica.Blocks.Sources.TimeTable solarBeam{0}(table = [{1}]);\n"
            .format(wall._idx, ZipTimeTable(timestamps, solarseries[wall._idx_zone_i])))
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
            .format(ZipTimeTable(series_time, series_diff)))
    fo.write("  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow loadDiff;\n")

### EQUATION
fo.write("equation\n")

for flow in flowlist:
    axis = "x"
    if flow.direction == 1: axis = "y"
    if flow.direction == 2: axis = "z"
    fo.write("  connect({0}.port_a, {1}.port_{2}1);\n".format(flow.name, "zone" + str(flow._idx_zone_i), axis))
    fo.write("  connect({0}.port_b, {1}.port_{2}2);\n".format(flow.name, "zone" + str(flow._idx_zone_o), axis))
for wall in walllist:
    axis = "x"
    if wall.direction == 1: axis = "y"
    if wall.direction == 2: axis = "z"
    fo.write("  connect({0}.port_a, {1}.port_{2}{3});\n".format(wall.name, "zone" + str(wall._idx_zone_i), 
        axis, str(wall._idx_port_i % 2 + 1)))
for flow in flowlist:
    fo.write("  connect({0}.i[1], {1}.o);\n".format(flow.name, "zone" + str(flow._idx_zone_o)))
    fo.write("  connect({0}.i[2], {1}.o);\n".format(flow.name, "zone" + str(flow._idx_zone_i)))

for zone in zonelist:
    if zone != None:
        if zone._idx < tickx * ticky:
            pass
            # fo.write("  connect({0}.port_s, heatFlow{1}.port);\n".format(zone.name, zone.idx))
            # fo.write("  connect(solarBeam{0}.y, heatFlow{1}.Q_flow);\n".format(zone.idx, zone.idx))
        for i in range(6):
            fo.write("  connect({0}.i[{1}], {2}.o);\n".format(zone.name, i+1, zone._port[i]))

for wall in walllist:
    fo.write("  connect({0}.i, {1}.o);\n".format(wall.name, "zone" + str(wall._idx_zone_i)))
    adjMaskValue = RetrieveAdjacencyMask(wall._idx_zone_i, tickx, ticky, mesh_mask, wall._idx_port_i)
    if wall.is_solarloaded:
        fo.write("  connect(wall{0}.port_s, loadBeam{1}.port);\n".format(wall._idx, wall._idx))
        fo.write("  connect(solarBeam{0}.y, loadBeam{1}.Q_flow);\n".format(wall._idx, wall._idx))
    else:
        fo.write("  connect(wall{0}.port_s, adiabaticS{1}.port);\n".format(wall._idx, wall._idx))
    if wall.direction != 2 and adjMaskValue == 0:
        fo.write("  connect(temp{0}.T, ambient.y);\n".format(wall._idx))
        fo.write("  connect(temp{0}.port, wallR{1}.port_b);\n".format(wall._idx, wall._idx))
        fo.write("  connect(wallR{0}.port_a, wall{1}.port_t);\n".format(wall._idx, wall._idx))
        fo.write("  connect(temp{0}.port, bound{1}.port_b);\n".format(wall._idx, wall._idx))
        fo.write("  connect(bound{0}.port_a, wall{1}.port_t);\n".format(wall._idx, wall._idx))
        # fo.write("  connect(wallC{0}.port, wallR{1}.port_b);\n".format(wall._idx, wall._idx))
    else:
        fo.write("  connect(adiabaticT{0}.port, wall{1}.port_t);\n".format(wall._idx, wall._idx))

    # to radiation box
    if RADIATED_MODE:
        fo.write("  connect(wall{0}.port_r, radbox.port_{0});\n".format(wall._idx))

if RADIATED_MODE:
    fo.write("  connect(radbox.port_s, loadDiff.port);\n")
    fo.write("  connect(solarDiff.y, loadDiff.Q_flow);\n")

# expriment configuration only for Wolfram
fo.write("annotation(experiment(StopTime = {0}, __Wolfram_NumberOfIntervals = -1));\n".format(time_simulation))
        
fo.write("end " + MODEL_NAME + ";\n")
fo.close()
