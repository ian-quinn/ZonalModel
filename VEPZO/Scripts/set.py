import math
import time
import datetime

import numpy as np
from mesh import GetMesh
from fake import FakeInput

# assitants
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
        content += str(time[i]) + ", " + str(value[i])
        if i + 1 != len(time):
            content += "; "
    return content

# module class
class Zone(object):
    def __init__(self, idx, x, y, z):
        self.coord = [x, y, z]
        self.idx = idx
        self.name = "zone" + str(idx)
        self.connect = [0, 0, 0, 0, 0, 0]
        self.port = ["", "", "", "", "", ""]
    def serialize(self):
        print("zone{0} ({1}, {2}, {3})".format(
            self.idx, self.coord[0], self.coord[1], self.coord[2]))

class Flow(object):
    def __init__(self, idx, direction, idxA, idxB, portA, portB):
        self.idxIn = idxA
        self.idxOut = idxB
        self.portIn = portA
        self.portOut = portB
        self.name = "flow" + str(idx)
        self.direction = direction
    def serialize(self):
        print("flow{0} direction{1} {2}/{3}".format(self.idxIn, self.direction, self.idxIn, self.idxOut))

class  Wall(object):
    def __init__(self, idx, direction, idxA, portA):
        self.idxIn = idxA       # index of the connected zone
        self.portIn = portA     # one of the zone port [0, 0, 0, 0, 0, 0]
        self.idx = idx          # index of this wall
        self.name = "wall" + str(idx)
        self.direction = direction
    def serialize(self):
        print("wall{0} direction{1} {2}".format(self.idx, self.direction, self.idxIn))
        
#################### USER SETUP ####################

model_name = "test"
# if you need to fake some day slices do set this as 6, 12, 18 or 24
fakeSlice = 0

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
time_start = datetime.datetime.strptime("2022-01-25 00:00:00", "%Y-%m-%d %H:%M:%S")
time_end = datetime.datetime.strptime("2022-01-25 23:59:59", "%Y-%m-%d %H:%M:%S")
time_delta = datetime.timedelta(seconds=time_step)

solarload = 250

# thermal property settings
capacity = 30000      # general heat capacity of wall material (J/K*m2) (cp * rho * thickness)
resistance = 0.5   # general thermal resistance of wall material (m2*K/W) (thickness / conductance + 1 / convection rate)
# lumped thermal property of perimeter enclosing?

epw_path = "Shanghai.epw"
delta_t = 0

######################################################

# good to go for the meshing work
# snapshots / 3d array in nested lists form representing solar load matrix of each step
# mesh_mask / 2d array representing solid cell in bounding box of target polygon
(snapshots, mesh_mask, cell_dim) = GetMesh(
    vertexloops, mask_wwr, mask_adia, scale_factor, 
    (latitude, longitude, utc), 
    (time_start, time_end, time_delta), 
    solarload, fakeSlice)

# passed down geometry information
dimx = cell_dim[0]
dimy = cell_dim[1]
dimz = 1                        # default as always
tickx = np.size(mesh_mask, 1)
ticky = np.size(mesh_mask, 0)
tickz = 3                       # default as always

time_simulation = (time_end - time_start).total_seconds()
if fakeSlice:
    time_simulation = 3600 * 24 * fakeSlice

# grab timestamps
timestamps = []          # time stamps in second
for i in range(len(snapshots)):
    timestamps.append(i * time_step)

# grab temperature seris from epw file
temp_ambient = ""
if fakeSlice:
    series_temp, actual_temp = FakeInput(epw_path, fakeSlice)
    series_time = [3600 * i for i in range(len(series_temp))]
    series_temp = series_temp + delta_t
    temp_init = series_temp[0]
    for i in range(len(series_temp)):
        temp_ambient = temp_ambient + str(series_time[i]) + ", " + str(series_temp[i])
        if i < len(series_temp) - 1:
            temp_ambient = temp_ambient + "; "
else:
    series_temp = []
    series_time = []
    timeticks = 0
    f = open(epw_path, mode="r")
    for line in f.readlines():
        datalist = line.split(',')
        if len(datalist) == 35:
            timeticks = timeticks + 1
            if timeticks * 3600 > (time_start - time_init).total_seconds() and \
                (timeticks - 1) * 3600 < (time_end - time_init).total_seconds():
                if len(series_temp) == 0:
                    series_time.append(0.0)
                    series_temp.append(round(float(datalist[6]) + 273.15 + delta_t, 2))
                else:
                    series_time.append((timeticks - 1) * 3600 - (time_start - time_init).total_seconds())
                    series_temp.append(round(float(datalist[6]) + 273.15 + delta_t, 2))
    temp_init = series_temp[0]

    for i in range(len(series_temp)):
        temp_ambient = temp_ambient + str(series_time[i]) + ", " + str(series_temp[i])
        if i < len(series_temp) - 1:
            temp_ambient = temp_ambient + "; "

print(temp_ambient)

# grab solar load
solardistributions = []
for snapshot in snapshots:
    solardistribution = [item for row in snapshot for item in row]
    solardistributions.append(solardistribution)
solarseries = np.array(solardistributions).transpose()

# indoor load


# module generation
# zone will cover all over the mesh regardelss of inside/outside the floorplan
# the index of each zone must projected to coordinates one and only
zonelist = [] # list of module Zone
for k in range(tickz):
    for j in range(ticky):
        for i in range(tickx):
            #zonelist.append(Zone(len(zonelist), i, j, k))
            if mesh_mask[ticky - j - 1][i] == 1:
                zonelist.append(Zone(len(zonelist), i, j, k))
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
    coords = zonelist[i].coord
    if coords[0] < tickx - 1:
        print("looping " + str(coords[0]) + " | " + str(coords[1]))
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
        print(chain)

walllist = [] # list of module Wall
for i in range(len(zonelist)):
    if zonelist[i] == None:
        continue
    connectors = zonelist[i].connect
    for j in range(len(connectors)):
        # any exposed connector (exposed) will be assigned with a wall
        if not connectors[j]:
            wallname = "wall" + str(len(walllist))
            walllist.append(Wall(len(walllist), j // 2, i, j))
            zonelist[i].port[j] = wallname

# locate the zone at the bottom and the top
adiabaticZonelist = [] # index list for top and bottom zones
for zone in zonelist:
    if zone != None:
        if zone.idx < tickx * ticky or zone.idx >= (tickz - 1) * tickx * ticky:
            adiabaticZonelist.append(zone.idx)


############################# module serialization #################################

# dump the mask matrix to csv for post-processing
print("#############" + str(type(mesh_mask[0][0])))
np.savetxt("mask.csv", mesh_mask, fmt='%d', delimiter=',')

fo = open(model_name + ".mo", "w")
fo.write("within VEPZO.Samples;\n")
fo.write("model " + model_name + "\n")
# fo.write("  import VEPZO;\n")
fo.write("  Modelica.Blocks.Sources.TimeTable ambient(table = [{0}]);\n".format(temp_ambient))
for zone in zonelist:
    if zone != None:
        if zone.idx < tickx * ticky:
            fo.write("  Zone {0}(dx = {1}, dy = {2}, dz = {3}, T_0 = {4});\n"
                .format(zone.name, dimx, dimy, dimz, temp_init))

            # fo.write("  Zone {0}(IsSource = true, dx = {1}, dy = {2}, dz = {3}, T_0 = {4});\n"
            #     .format(zone.name, dimx, dimy, dimz, temp_init))
            # fo.write("  Modelica.Blocks.Sources.TimeTable solar{0}(table = [{1}]);\n"
            #     .format(zone.idx, ZipTimeTable(timestamps, solarseries[zone.idx])))
            # fo.write("  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow heatFlow{0};\n".format(zone.idx))
        else:
            fo.write("  Zone {0}(dx = {1}, dy = {2}, dz = {3}, T_0 = {4});\n"
                .format(zone.name, dimx, dimy, dimz, temp_init))

for flow in flowlist:
    fo.write("  Flow {0}(Direction = {1});\n".format(flow.name, flow.direction))
for wall in walllist:
    print("Wall-" + str(walllist.index(wall)) + " assigned boundary condition ", end='')
    area = dimy * dimz
    if wall.direction == 1:
        area = dimx * dimz
    else:
        area = dimx * dimy

    # type 1, on the floor, assigned with port_s from solar load
    # type 2, on outer boundary, assigned with port_t from prescribed env temperature
    # type 3, elsewhere, adiabatic with capacity only
    if wall.idxIn < tickx * ticky and wall.direction == 2:
        fo.write("  Wall {0}(IsSource = true, Direction = {1}, R = {2}, C = {3}, T_0 = {4});\n"
            .format(wall.name, wall.direction, resistance / area, capacity * area, temp_init))
        fo.write("  Modelica.Blocks.Sources.TimeTable solar{0}(table = [{1}]);\n"
            .format(wall.idx, ZipTimeTable(timestamps, solarseries[wall.idxIn])))
        fo.write("  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow heatFlow{0};\n".format(wall.idx))
    elif wall.direction != 2 and RetrieveAdjacencyMask(wall.idxIn, tickx, ticky, mesh_mask, wall.portIn) == 0:
        fo.write("  Wall {0}(IsAdiabatic = false, Direction = {1}, R = {2}, C = {3}, T_0 = {4});\n"
            .format(wall.name, wall.direction, resistance / area, capacity * area, temp_init))
        fo.write("  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature temp{0};\n".format(wall.idx))
        fo.write("  Modelica.Thermal.HeatTransfer.Components.ThermalResistor wallR{0}(R = {1});\n"
            .format(wall.idx, 0.001))
        fo.write("  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor wallC{0}(C = {1}, T.start = {2});\n"
            .format(wall.idx, 100000000, temp_init))
    else:
        fo.write("  Wall {0}(Direction = {1}, R = {2}, C = {3}, T_0 = {4});\n"
            .format(wall.name, wall.direction, resistance / area, capacity * area, temp_init))

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

fo.write("equation\n")
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
            # fo.write("  connect(solar{0}.y, heatFlow{1}.Q_flow);\n".format(zone.idx, zone.idx))
        for i in range(6):
            fo.write("  connect({0}.i[{1}], {2}.o);\n".format(zone.name, i+1, zone.port[i]))
for wall in walllist:
    fo.write("  connect({0}.i, {1}.o);\n".format(wall.name, "zone" + str(wall.idxIn)))
    adjMaskValue = RetrieveAdjacencyMask(wall.idxIn, tickx, ticky, mesh_mask, wall.portIn)
    if wall.idxIn < tickx * ticky and wall.direction == 2:
        fo.write("  connect(wall{0}.port_s, heatFlow{1}.port);\n".format(wall.idx, wall.idx))
        fo.write("  connect(solar{0}.y, heatFlow{1}.Q_flow);\n".format(wall.idx, wall.idx))
    elif wall.direction != 2 and adjMaskValue == 0:
        fo.write("  connect(temp{0}.T, ambient.y);\n".format(wall.idx))
        fo.write("  connect(temp{0}.port, wallR{1}.port_b);\n".format(wall.idx, wall.idx))
        fo.write("  connect(wallR{0}.port_a, wall{1}.port_t);\n".format(wall.idx, wall.idx))
        fo.write("  connect(wallC{0}.port, wallR{1}.port_a);\n".format(wall.idx, wall.idx))
    else:
        pass
    # if not (wall.idxIn in adiabaticZonelist and wall.direction == 2) and adjMaskValue == 0:
    #     fo.write("  connect(temp{0}.T, ambient.y);\n".format(wall.idx))
    #     fo.write("  connect(temp{0}.port, wallR{1}.port_b);\n".format(wall.idx, wall.idx))
    #     fo.write("  connect(wallR{0}.port_a, wall{1}.port_s);\n".format(wall.idx, wall.idx))
    #     fo.write("  connect(wallC{0}.port, wallR{1}.port_a);\n".format(wall.idx, wall.idx))
# expriment configuration only for Wolfram
fo.write("annotation(experiment(StopTime = {0}, __Wolfram_NumberOfIntervals = -1));\n".format(time_simulation))
        
fo.write("end " + model_name + ";\n")
fo.close()
