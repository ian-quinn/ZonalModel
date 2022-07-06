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
class Point(object):
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
    def __repr__(self):
        return "Point({0}, {1}, {2})".format(self.x, self.y, self.z)

class Zone(object):
    def __init__(self, idx, x, y, z):
        self.idx = idx
        self.coords = [x, y, z]
        self.name = "zone" + str(idx)
        self.connect = [0, 0, 0, 0, 0, 0]
        self.port = ["", "", "", "", "", ""]
    def __repr__(self):
        print("zone{0} ({1}, {2}, {3})".format(
            self.idx, self.coords[0], self.coords[1], self.coords[2]))

class Flow(object):
    def __init__(self, idx, direction, idxA, idxB, portA, portB):
        self.idxIn = idxA
        self.idxOut = idxB
        self.portIn = portA
        self.portOut = portB
        self.name = "flow" + str(idx)
        self.direction = direction
    def __repr__(self):
        print("flow{0} direction{1} {2}/{3}".format(self.idxIn, self.direction, self.idxIn, self.idxOut))

class  Wall(object):
    def __init__(self, idx, direction, idxA, portA, pt):
        self.idx = idx          # index of this wall
        self.name = "wall" + str(idx)
        self.direction = direction
        self.idxIn = idxA       # index of the connected zone
        self.portIn = portA     # one of the zone port [0, 0, 0, 0, 0, 0]
        self.basepoint = pt     # np.array coords of the wall basepoint (min xyz)
    def __repr__(self):
        print("wall{0} direction{1} {2}".format(self.idx, self.direction, self.idxIn))
        
#################### USER SETUP ####################

# if you need to fake some day slices do set this as 6, 12, 18 or 24
fakeSlice = 0
isRadiated = True

model_name = "test"
vfm_path = "view_21_3.csv"
epw_path = "../Weather/Shanghai.epw"


# DIAMOND
# vertexloops = [[[3, 0], [10, 0], [10, 7], [7, 10], [0, 10], [0, 3], [3, 0]], 
#    [[4, 4], [6, 4], [6, 6], [4, 6], [4, 4]]]

# RAND
# vertexloops = [[[12, 0], [30, 0], [30, 15], [18, 15], [18, 30], [0, 30], [0, 15], [12, 15], [12, 0]]]

# SNAKE
# vertexloops = [[[0, 0], [33, 0], [33, 9], [6, 9], [6, 12], [33, 12], [33, 33], [0, 33], [0, 24], \
#    [27, 24], [27, 21], [0, 21], [0, 0]]]


# FRAME
# vertexloops = [[[0, 0], [21, 0], [21, 21], [0, 21], [0, 0]], \
#     [[6, 6], [15, 6], [15, 15], [6, 15], [6, 6]]]
# vertexloops = [[[0, 0], [7, 0], [7, 7], [0, 7], [0, 0]], \
#     [[2, 2], [5, 2], [5, 5], [2, 5], [2, 2]]]
vertexloops = [[[0, 0], [21, 0], [21, 21], [0, 21], [0, 0]]]

mask_wwr = [0.8]
# mask_adia = [0]
# mask_wwr = [0.8, 0.8, 0.8, 0.8]
mask_adia = [0, 0, 0, 0]

scale_factor = 3

latitude = 31.17
longitude = 121.43
utc = 8      # time zone

time_step = 600                 # in seconds
time_init = datetime.datetime.strptime("2022-01-01 00:00:00", "%Y-%m-%d %H:%M:%S")
time_start = datetime.datetime.strptime("2022-05-15 00:00:00", "%Y-%m-%d %H:%M:%S")
time_end = datetime.datetime.strptime("2022-06-30 23:59:59", "%Y-%m-%d %H:%M:%S")
time_delta = datetime.timedelta(seconds=time_step)

solarload = 800
SHGC = 0.3

# thermal property settings
C_win = 10000      # general heat capacity of wall material (J/K*m2) (cp * rho * thickness)
R_win = 1   # general thermal resistance of wall material (m2*K/W) (thickness / conductance + 1 / convection rate)
# lumped thermal property of perimeter enclosing?
C_wall = 20000
R_wall = 1

delta_t = 0

##################### SOME PREPROCESSING ###################
# just for now, include several geometry stuff
outerpts = []
offset = np.array([scale_factor, scale_factor])
for coords in vertexloops[0]:
    outerpts.append(np.array(coords) + offset)


##################### MESHING WORK ###################

# good to go for the meshing work
# snapshots / 3d array in nested lists form representing solar load matrix of each step
# mesh_mask / 2d array representing solid cell in bounding box of target polygon
(snapshots, mesh_mask, cell_dim) = GetMesh(
    vertexloops, mask_wwr, mask_adia, scale_factor, 
    (latitude, longitude, utc), 
    (time_start, time_end, time_delta), 
    solarload, SHGC, epw_path, fakeSlice)

# passed down geometry information
dimx = cell_dim[0]
dimy = cell_dim[1]
dimz = 1                        # default as always
dims = np.array([dimx, dimy, dimz])
# print(dims)
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

time_simulation = (time_end - time_start).total_seconds()
if fakeSlice:
    time_simulation = 3600 * 24 * fakeSlice

# grab timestamps
# the timestamps returned by GetMesh() should correspond to the following series_time
# make sure all functions follow the time_step setting
# the question is, do we all need a small time step?
timestamps = []          # time stamps in second
for i in range(len(snapshots)):
    timestamps.append(i * time_step)


# grab temperature series from epw file
if fakeSlice:
    series_temp, actual_temp = FakeInput(epw_path, fakeSlice, time_step)
    series_time = [time_step * i for i in range(len(series_temp))]
    series_temp = series_temp + delta_t
    temp_init = series_temp[0]
    series_diff = []
    roomarea = tickx * ticky * scale_factor *scale_factor
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
    f = open(epw_path, mode="r")
    for line in f.readlines():
        datalist = line.split(',')
        if len(datalist) == 35:
            timeticks = timeticks + 1
            if timeticks * 3600 > (time_start - time_init).total_seconds() and \
                (timeticks - 1) * 3600 < (time_end - time_init).total_seconds():
                series_temp.append(round(float(datalist[6]) + 273.15 + delta_t, 2))
                series_diff.append(int(datalist[14]) * SHGC)
                if len(series_temp) == 0:
                    series_time.append(0.0)
                else:
                    series_time.append((timeticks - 1) * 3600 - (time_start - time_init).total_seconds())
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
            #zonelist.append(Zone(len(zonelist), i, j, k))
            if mesh_mask[ticky - j - 1][i] == 1:
                zonelist.append(Zone(len(zonelist), i, j, k))
                # print(zonelist[-1].coords)
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
    coords = zonelist[i].coords
    if coords[0] < tickx - 1:
        # print("looping " + str(coords[0]) + " | " + str(coords[1]))
        if mesh_mask[ticky - coords[1] - 1][coords[0] + 1] == 1:
            targetId = RetrieveId([coords[0] + 1, coords[1], coords[2]], tickx, ticky)
            flowlist.append(Flow(len(flowlist), 0, targetId, i, 0, 1))
            zonelist[i].connect[1] = 1
            zonelist[i].port[1] = flowlist[-1].name
            zonelist[targetId].connect[0] = 1
            zonelist[targetId].port[0] = zonelist[i].port[1]
    if coords[1] < ticky - 1:
        if mesh_mask[ticky - coords[1] - 2][coords[0]] == 1:
            flowlist.append(Flow(len(flowlist), 1, 
                RetrieveId([coords[0], coords[1] + 1, coords[2]], tickx, ticky), i, 2, 3))
            zonelist[i].connect[3] = 1
            zonelist[i].port[3] = flowlist[-1].name
            zonelist[RetrieveId([coords[0], coords[1] + 1, coords[2]], tickx, ticky)].connect[2] = 1
            zonelist[RetrieveId([coords[0], coords[1] + 1, coords[2]], tickx, ticky)].port[2] = zonelist[i].port[3]
    if coords[2] < tickz - 1:
        flowlist.append(Flow(len(flowlist), 2, 
            RetrieveId([coords[0], coords[1], coords[2] + 1], tickx, ticky), i, 4, 5))
        zonelist[i].connect[5] = 1
        zonelist[i].port[5] = flowlist[-1].name
        zonelist[RetrieveId([coords[0], coords[1], coords[2] + 1], tickx, ticky)].connect[4] = 1
        zonelist[RetrieveId([coords[0], coords[1], coords[2] + 1], tickx, ticky)].port[4] = zonelist[i].port[5]

for zone in zonelist:
    if zone != None:
        chain = ""
        for i in range(len(zone.connect)):
            chain += ", " + str(zone.connect[i])
        # print(chain)

walllist = [] # list of module Wall
for i in range(len(zonelist)):
    if zonelist[i] == None:
        continue
    connectors = zonelist[i].connect
    for j in range(len(connectors)):
        # any exposed connector (exposed) will be assigned with a wall
        if not connectors[j]:
            wallname = "wall" + str(len(walllist))
            # stretch zone index to get the base point of the zone
            basepoint = np.array(zonelist[i].coords) * dims
            # if the wall has negative normal vector
            # plus the dimension size accordingly to get its location point
            basepoint[j // 2] += (j % 2) * dims[j // 2]
            walllist.append(Wall(len(walllist), j // 2, i, j, basepoint))
            zonelist[i].port[j] = wallname

#################### VIEW MATRIX #####################

if isRadiated:
    wallnormals = []
    wallbasepts = []
    for wall in walllist:
        wallbasepts.append(wall.basepoint)
        wallnormal = [0, 0, 0]
        wallnormal[wall.portIn // 2] = 1 - wall.portIn % 2 * 2
        wallnormals.append(np.array(wallnormal))

    if vfm_path:
        view_matrix = np.genfromtxt(vfm_path, delimiter=',')
    else:
        view_matrix = GetViewFactor(wallbasepts, wallnormals, dims, 0.5)

    # locate transparent walls
    wallistransparent = []
    for i in range(len(walllist)):
        basept = np.array([walllist[i].basepoint[0], walllist[i].basepoint[1]])
        vec_norm = np.array([wallnormals[i][0], wallnormals[i][1]])
        # skim out horizontal surface
        if np.linalg.norm(vec_norm) == 0:
            wallistransparent.append(False)
            continue
        for j in range(len(outerpts) - 1):
            if j + 1 <= len(mask_wwr):
                if mask_wwr[j] > 0:
                    vec_edge = outerpts[j + 1] - outerpts[j]
                    check_verticle = np.dot(vec_edge, vec_norm)
                    check_parallel = np.cross(outerpts[j + 1] - basept, outerpts[j] - basept)
                    if check_verticle == 0 and check_parallel == 0:
                        wallistransparent.append(True)
                    else:
                        wallistransparent.append(False)

# locate the zone at the bottom and the top
adiabaticZonelist = [] # index list for top and bottom zones
for zone in zonelist:
    if zone != None:
        if zone.idx < tickx * ticky or zone.idx >= (tickz - 1) * tickx * ticky:
            adiabaticZonelist.append(zone.idx)


######################################################################################
######################################################################################
################################### SERIALIZATION ####################################


# dump the mask matrix to csv for post-processing
np.savetxt("mask.csv", mesh_mask, fmt='%d', delimiter=',')
if not vfm_path and isRadiated:
    np.savetxt("view.csv", view_matrix, fmt='%2f', delimiter=',')


##################### FILE RADIATION MODEL ###################
if isRadiated:
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
        if wallistransparent[i]:
            fo.write("  Qout_{0} = {1} * Modelica.Constants.sigma * port_{0}.T ^ 4 + port_s.Q_flow + {2} * Qin_{0};\n"
                .format(i, alpha_infra, tao_infra))

        else:
            fo.write("  Qout_{0} = {1} * Modelica.Constants.sigma * port_{0}.T ^ 4 + {2} * Qin_{0};\n"
                .format(i, alpha_infra, rho_infra))
        fo.write("  Qin_{0} = {1};\n".format(i, Sigma))
        fo.write("  port_{0}.Q_flow = Qin_{0} - Qout_{0};\n".format(i))
    fo.write("end Radbox;\n")
    fo.close()


##################### FILE ROOM MODEL ###################
fo = open("../Samples/" + model_name + ".mo", "w")
fo.write("within VEPZO.Samples;\n")
fo.write("model " + model_name + "\n")

##################### MODULE DECLARATION ###################
fo.write("  Modelica.Blocks.Sources.TimeTable ambient(table = [{0}]);\n"
    .format(ZipTimeTable(series_time, series_temp)))
for zone in zonelist:
    if zone != None:
        if zone.idx < tickx * ticky:
            fo.write("  Zone {}(dx = {}, dy = {}, dz = {}, T_0 = {:.2f});\n"
                .format(zone.name, dimx, dimy, dimz, temp_init))

            # fo.write("  Zone {0}(IsSource = true, dx = {1}, dy = {2}, dz = {3}, T_0 = {4});\n"
            #     .format(zone.name, dimx, dimy, dimz, temp_init))
            # fo.write("  Modelica.Blocks.Sources.TimeTable solar{0}(table = [{1}]);\n"
            #     .format(zone.idx, ZipTimeTable(timestamps, solarseries[zone.idx])))
            # fo.write("  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow heatFlow{0};\n".format(zone.idx))
        else:
            fo.write("  Zone {}(dx = {}, dy = {}, dz = {}, T_0 = {:.2f});\n"
                .format(zone.name, dimx, dimy, dimz, temp_init))

for flow in flowlist:
    fo.write("  Flow {0}(Direction = {1});\n".format(flow.name, flow.direction))
for wall in walllist:
    # print("Wall-" + str(walllist.index(wall)) + " assigned boundary condition ", end='')
    area = dimy * dimz
    if wall.direction == 1:
        area = dimx * dimz
    else:
        area = dimx * dimy

    # type 1, on the floor, assigned with port_s from solar load
    # type 2, on outer boundary, assigned with port_t from prescribed env temperature
    # type 3, elsewhere, adiabatic with capacity only
    if wall.idxIn < tickx * ticky and wall.direction == 2:
        fo.write("  Wall {}(IsSource = true, IsRadiated = true, Direction = {}, R = {:.2f}, C = {}, T_0 = {:.2f});\n"
            .format(wall.name, wall.direction, R_wall / area, C_wall * area, temp_init))
        fo.write("  Modelica.Blocks.Sources.TimeTable solarBeam{0}(table = [{1}]);\n"
            .format(wall.idx, ZipTimeTable(timestamps, solarseries[wall.idxIn])))
        fo.write("  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow loadBeam{0};\n".format(wall.idx))
    elif wall.direction != 2 and RetrieveAdjacencyMask(wall.idxIn, tickx, ticky, mesh_mask, wall.portIn) == 0:
        if walllist.index(wall) in wallistransparent:
            fo.write("  Wall {}(IsAdiabatic = false, IsRadiated = true, Direction = {}, R = {:.2f}, C = {}, T_0 = {:.2f});\n"
                .format(wall.name, wall.direction, R_win / area, C_win * area, temp_init))
        else:
            fo.write("  Wall {}(IsAdiabatic = false, IsRadiated = true, Direction = {}, R = {:.2f}, C = {}, T_0 = {:.2f});\n"
                .format(wall.name, wall.direction, R_wall / area, C_wall * area, temp_init))
        fo.write("  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature temp{};\n".format(wall.idx))
        fo.write("  Modelica.Thermal.HeatTransfer.Components.ThermalResistor wallR{}(R = {:.2f});\n"
            .format(wall.idx, 1 / 7.6 / area))
        fo.write("  Bound bound{0};\n".format(wall.idx))
        # fo.write("  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor wallC{0}(C = {1}, T.start = {2});\n"
        #     .format(wall.idx, 1000000, temp_init))
    else:
        if walllist.index(wall) in wallistransparent:
            fo.write("  Wall {}(IsRadiated = true, Direction = {}, R = {:.2f}, C = {}, T_0 = {:.2f});\n"
                .format(wall.name, wall.direction, R_win / area, C_win * area, temp_init))
        else:
            fo.write("  Wall {}(IsRadiated = true, Direction = {}, R = {:.2f}, C = {}, T_0 = {:.2f});\n"
                .format(wall.name, wall.direction, R_wall / area, C_wall * area, temp_init))

    # if wall.idxIn in adiabaticZonelist and wall.direction == 2:
    #     print("Adiabatic (top/bottom)")
    #     fo.write("  Wall {0}(Direction = {1}, R = {2}, C = {3});\n".format(wall.name, wall.direction, 
    #         resistance / area, capacity * area))
    #     continue
    # if RetrieveAdjacencyMask(wall.idxIn, tickx, ticky, mesh_mask, wall.portIn) == 2:
    #     print("Adiabatic (hole)")
    #     fo.write("  Wall {0}(Direction = {1}, R = {2}, C = {3});\n".format(wall.name, wall.direction, 
    #         resistance / area, capacity * area))
    #     continue
    # if RetrieveAdjacencyMask(wall.idxIn, tickx, ticky, mesh_mask, wall.portIn) == 3:
    #     print("Adiabatic (manual)")
    #     fo.write("  Wall {0}(Direction = {1}, R = {2}, C = {3});\n".format(wall.name, wall.direction, 
    #         resistance / area, capacity * area))
    #     continue
    # print("Temperature (epw)")
    # fo.write("  Wall {0}(Direction = {1}, IsSource = true);\n".format(wall.name, wall.direction))
    # fo.write("  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature temp{0};\n".format(wall.idx))
    # fo.write("  Modelica.Thermal.HeatTransfer.Components.ThermalResistor wallR{0}(R = {1});\n"
    #     .format(wall.idx, resistance / area))
    # fo.write("  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor wallC{0}(C = {1}, T.start = {2});\n"
    #     .format(wall.idx, capacity * area, temp_init))

# include the radiation box model
if isRadiated:
    fo.write("  Radbox radbox;\n")
    fo.write("  Modelica.Blocks.Sources.TimeTable solarDiff(table = [{0}]);\n"
            .format(ZipTimeTable(series_time, series_diff)))
    fo.write("  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow loadDiff;\n")

##################### EQUATION DECLARATION ###################
fo.write("equation\n")
if isRadiated:
    fo.write("  connect(radbox.port_s, loadDiff.port);\n")
    fo.write("  connect(solarDiff.y, loadDiff.Q_flow);\n")

for flow in flowlist:
    axis = "x"
    if flow.direction == 1: axis = "y"
    if flow.direction == 2: axis = "z"
    fo.write("  connect({0}.port_a, {1}.port_{2}1);\n".format(flow.name, "zone" + str(flow.idxIn), axis))
    fo.write("  connect({0}.port_b, {1}.port_{2}2);\n".format(flow.name, "zone" + str(flow.idxOut), axis))
for wall in walllist:
    axis = "x"
    if wall.direction == 1: axis = "y"
    if wall.direction == 2: axis = "z"
    fo.write("  connect({0}.port_a, {1}.port_{2}{3});\n".format(wall.name, "zone" + str(wall.idxIn), 
        axis, str(wall.portIn % 2 + 1)))
for flow in flowlist:
    fo.write("  connect({0}.i[1], {1}.o);\n".format(flow.name, "zone" + str(flow.idxOut)))
    fo.write("  connect({0}.i[2], {1}.o);\n".format(flow.name, "zone" + str(flow.idxIn)))
for zone in zonelist:
    if zone != None:
        if zonelist.index(zone) < tickx * ticky:
            pass
            # fo.write("  connect({0}.port_s, heatFlow{1}.port);\n".format(zone.name, zone.idx))
            # fo.write("  connect(solarBeam{0}.y, heatFlow{1}.Q_flow);\n".format(zone.idx, zone.idx))
        for i in range(6):
            fo.write("  connect({0}.i[{1}], {2}.o);\n".format(zone.name, i+1, zone.port[i]))
for wall in walllist:
    fo.write("  connect({0}.i, {1}.o);\n".format(wall.name, "zone" + str(wall.idxIn)))
    adjMaskValue = RetrieveAdjacencyMask(wall.idxIn, tickx, ticky, mesh_mask, wall.portIn)
    if wall.idxIn < tickx * ticky and wall.direction == 2:
        fo.write("  connect(wall{0}.port_s, loadBeam{1}.port);\n".format(wall.idx, wall.idx))
        fo.write("  connect(solarBeam{0}.y, loadBeam{1}.Q_flow);\n".format(wall.idx, wall.idx))
    elif wall.direction != 2 and adjMaskValue == 0:
        fo.write("  connect(temp{0}.T, ambient.y);\n".format(wall.idx))
        fo.write("  connect(temp{0}.port, wallR{1}.port_b);\n".format(wall.idx, wall.idx))
        fo.write("  connect(wallR{0}.port_a, wall{1}.port_t);\n".format(wall.idx, wall.idx))
        fo.write("  connect(temp{0}.port, bound{1}.port_b);\n".format(wall.idx, wall.idx))
        fo.write("  connect(bound{0}.port_a, wall{1}.port_t);\n".format(wall.idx, wall.idx))
        # fo.write("  connect(wallC{0}.port, wallR{1}.port_b);\n".format(wall.idx, wall.idx))
    else:
        pass

    # to radiation box
    if isRadiated:
        fo.write("  connect(wall{0}.port_r, radbox.port_{0});\n".format(wall.idx))

    # if not (wall.idxIn in adiabaticZonelist and wall.direction == 2) and adjMaskValue == 0:
    #     fo.write("  connect(temp{0}.T, ambient.y);\n".format(wall.idx))
    #     fo.write("  connect(temp{0}.port, wallR{1}.port_b);\n".format(wall.idx, wall.idx))
    #     fo.write("  connect(wallR{0}.port_a, wall{1}.port_s);\n".format(wall.idx, wall.idx))
    #     fo.write("  connect(wallC{0}.port, wallR{1}.port_a);\n".format(wall.idx, wall.idx))

# expriment configuration only for Wolfram
fo.write("annotation(experiment(StopTime = {0}, __Wolfram_NumberOfIntervals = -1));\n".format(time_simulation))
        
fo.write("end " + model_name + ";\n")
fo.close()
