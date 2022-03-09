import math
import time
import datetime

import numpy as np
from shapely.geometry import Polygon
from matplotlib import pyplot as plt
import matplotlib.animation as animation

# solar calculator
def getSEA(struct_time, latitude, longitude, utc_offset):
    date = struct_time
    # print('date', date)
    hour = date[3]
    minute = date[4]
    # Check your timezone to add the offset
    hour_minute = (hour + minute / 60) - utc_offset
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

def getAZ(struct_time, latitude, longitude, utc_offset):
    date = struct_time
    hour = date[3]
    minute = date[4]
    # Check your timezone to add the offset
    hour_minute = (hour + minute / 60) - utc_offset
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


# assitants
def RetrieveId(coord):
    return x * y * coord[2] + y * coord[1] + coord[0]
def RetrieveCoord(id):
    return [id % (x * y) % x, id % (x * y) // x, id // (x * y)]
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
        self.idxOut = idxA
        self.portIn = portA
        self.portOut = portB
        self.name = "flow" + str(idx)
        self.direction = direction
    def serialize(self):
        print("flow{0} direction{1} {2}/{3}".format(self.idxIn, self.direction, self.idxIn, self.idxOut))
class  Wall(object):
    def __init__(self, idx, direction, idxA, portA):
        self.idxIn = idxA
        self.portIn = portA
        self.idx = idx
        self.name = "wall" + str(idx)
        self.direction = direction
    def serialize(self):
        print("wall{0} direction{1} {2}".format(self.idx, self.direction, self.idxIn))
        

# global settings for user input
visualizeShadow = False      # generate a movie flip for solar distribution
useActuralAmbient = True    # use ambient temperature from EPW file

x = 10      # grid on abscissa
y = 10      # grid on ordinate
z = 3       # grid on vertical

dim_x = 30  # abscissa dimension of the box
dim_y = 30  # ordinate dimension of the box
dim_z = 3   # vertical dimension of the box

gap = 0.5   # gap between the sub-surface window to its hosting surface
sill = 1    # height of the window sill
h = 1.5     # window height
Q = 1200    # general beam radiation intensity

capacity = 800      # general heat capacity of wall material
resistance = 0.26   # general thermal resistance of wall material
#ambientTemp = "0, 278.25; 3600, 278.95; 7200, 280.45; 10800, 281.55; 14400, 282.35; 18000, 282.75; 21600, 282.95; 25200, 283.15; 28800, 283.35; 32400, 283.55; 36000, 283.35; 39600, 282.75; 43200, 281.95"
ambientTemp = "0, 276.95; 3600, 276.65; 7200, 276.35; 10800, 275.95; 14400, 275.55; 18000, 275.25; 21600, 275.15; 25200, 275.25; 28800, 275.35; 32400, 275.55; 36000, 275.85; 39600, 276.25; 43200, 276.85; 46800, 277.65; 50400, 278.65; 54000, 279.75; 57600, 280.35; 61200, 280.15; 64800, 279.55; 68400, 279.55; 68400, 278.85; 72000, 278.15; 75600, 277.65; 79200, 277.35; 82800, 277.15; 86400, 277.15"

time_start = datetime.datetime.strptime("2022-03-07 00:00:00", "%Y-%m-%d %H:%M:%S")
time_end = datetime.datetime.strptime("2022-03-07 23:59:59", "%Y-%m-%d %H:%M:%S")
latitude = 33
longitude = 122
utc_offset = 8      # time zone
timestep = 10       # interval to update the solar radiation (in minute)


# shadow calculation
delta = datetime.timedelta(minutes=timestep)
zeltas = []         # in degree
thetas = []         # in degree
times = []          # time stamps in second
time_ellapse = 0    # accumulated time in second
time_current = time_start
while time_current < time_end:
    times.append(time_ellapse)
    time_current += delta
    time_ellapse += timestep * 60
    thetas.append(getSEA(time_current.timetuple(), latitude, longitude, utc_offset))
    zeltas.append(getAZ(time_current.timetuple(), latitude, longitude, utc_offset))

grids = [] # list of shapely Polygon for radiation calculation
solars = [ [] for i in range(x * y) ]
for i in range(x * y):
    row = i % x
    col = i // x
    grids.append(Polygon([
        [row * dim_x / x, col * dim_y / y], 
        [(row + 1) * dim_x / x, col * dim_y / y], 
        [(row + 1) * dim_x / x, (col + 1) * dim_y / y], 
        [row * dim_x / x, (col + 1) * dim_y / y]
        ]))
for i in range(len(zeltas)):
    # initialize to zero
    for j in range(len(grids)):
            solars[j].append(0)
    if zeltas[i] > 90 and zeltas[i] < 270 and thetas[i] > 0:
        project_x = math.tan(math.radians(thetas[i])) * math.tan(math.radians(zeltas[i] - 90))
        project_y = math.tan(math.radians(thetas[i]))
        shadow = Polygon([
            [gap - sill / project_x, sill / project_y], 
            [dim_x - gap - sill / project_x, sill / project_y], 
            [dim_x - gap - (sill + h) / project_x, (sill + h) / project_y], 
            [gap - (sill + h) / project_x, (sill + h) / project_y]
            ])
        area = (dim_x - 2 * gap) * h / math.tan(math.radians(thetas[i]))
        # print(list(shadow.exterior.coords))
        for j in range(len(grids)):
            sect = shadow.intersection(grids[j])
            # print(list(sect.exterior.coords))
            # print(sect.area)
            solars[j][i] += sect.area / area * Q * math.sin(math.radians(thetas[i])) * math.sin(math.radians(zeltas[i] - 90))
    if zeltas[i] > 0 and zeltas[i] < 180 and thetas[i] > 0:
        project_x = math.tan(math.radians(thetas[i])) * math.tan(math.radians(zeltas[i]))
        project_y = math.tan(math.radians(thetas[i]))
        shadow = Polygon([
            [dim_x - sill / project_y, gap - sill / project_x], 
            [dim_y - sill / project_y, dim_y - gap - sill / project_x], 
            [dim_y - (sill + h) / project_y, dim_y - gap - (sill + h) / project_x], 
            [dim_x - (sill + h) / project_y, gap - (sill + h) / project_x]
            ])
        area = (dim_y - 2 * gap) * h / math.tan(math.radians(thetas[i]))
        for j in range(len(grids)):
            sect = shadow.intersection(grids[j])
            solars[j][i] += sect.area / area * Q * math.sin(math.radians(thetas[i])) * math.sin(math.radians(zeltas[i]))
    if zeltas[i] > 0 and zeltas[i] < 90 and thetas[i] > 0 or zeltas[i] > 270 and zeltas[i] < 360 and thetas[i] > 0:
        zelta = zeltas[i] + 90
        if zeltas[i] > 270 and zeltas[i] < 360:
            zelta = zeltas[i] - 270 
        project_x = math.tan(math.radians(thetas[i])) * math.tan(math.radians(zelta))
        project_y = math.tan(math.radians(thetas[i]))
        shadow = Polygon([
            [dim_x - gap - sill / project_x, dim_y - sill / project_y], 
            [gap - sill / project_x, dim_y - sill / project_y], 
            [gap - (sill + h) / project_x, dim_y - (sill + h) / project_y], 
            [dim_x - gap - (sill + h) / project_x, dim_y - (sill + h) / project_y]
            ])
        area = (dim_x - 2 * gap) * h / math.tan(math.radians(thetas[i]))
        for j in range(len(grids)):
            sect = shadow.intersection(grids[j])
            solars[j][i] += sect.area / area * Q * math.sin(math.radians(thetas[i])) * math.sin(math.radians(zelta))
    if zeltas[i] > 180 and zeltas[i] < 360 and thetas[i] > 0:
        project_x = math.tan(math.radians(thetas[i])) * math.tan(math.radians(zeltas[i] - 180))
        project_y = math.tan(math.radians(thetas[i]))
        shadow = Polygon([
            [sill / project_y, dim_y - gap - sill / project_x], 
            [sill / project_y, gap - sill / project_x], 
            [(sill + h) / project_y, gap - (sill + h) / project_x], 
            [(sill + h) / project_y, dim_y - gap - (sill + h) / project_x]
            ])
        area = (dim_y - 2 * gap) * h / math.tan(math.radians(thetas[i]))
        for j in range(len(grids)):
            sect = shadow.intersection(grids[j])
            solars[j][i] += sect.area / area * Q * math.sin(math.radians(thetas[i])) * math.sin(math.radians(zeltas[i] - 180))


# visualize radiation intensity on the floorplan
if visualizeShadow:
    snapshots = [] # list of 2D numpy matrix for time series of solar distribution
    maxRadiation = 0
    for i in range(len(zeltas)):
        temp = []
        for j in range(len(solars)):
            temp.append(solars[j][i])
            if solars[j][i] > maxRadiation:
                maxRadiation = solars[j][i]
        snapshot = np.array(temp).reshape((x, y))
        snapshots.append(snapshot)

    fig = plt.figure()
    xlattice = np.arange(0, x + 1)
    ylattice = np.arange(0, y + 1)
    imgs = []
    for snapshot in snapshots:
        imgs.append((plt.pcolor(xlattice, ylattice, snapshot, norm=plt.Normalize(0, maxRadiation)),))
    img_ani = animation.ArtistAnimation(fig, imgs, interval = 50, repeat_delay=0, blit=True)
    plt.show()

    # pip install pillow
    # https://stackoverflow.com/questions/51512141/how-to-make-matplotlib-saved-gif-looping
    class LoopingPillowWriter(animation.PillowWriter):
        def finish(self):
            self._frames[0].save(
                self._outfile, save_all=True, append_images=self._frames[1:],
                duration=int(1000 / self.fps), loop=0)
    img_ani.save(r"animation.gif", writer=LoopingPillowWriter(fps=20)) 
    # img_ani.save(r"animation.gif", writer=animation.PillowWriter(fps=30))
    # img_ani.save(r"animation.gif", writer='imagemagick')


# module generation
zonelist = [] # list of module Zone
for k in range(z):
    for j in range(y):
        for i in range(x):
            zonelist.append(
                Zone(len(zonelist), i, j, k))

flowlist = [] # list of module Flow
for i in range(len(zonelist)):
    coords = RetrieveCoord(i)
    if coords[0] < x - 1:
        targetId = RetrieveId([coords[0] + 1, coords[1], coords[2]])
        flowlist.append(Flow(len(flowlist), 0, targetId, i, 0, 1))
        zonelist[i].connect[1] = 1
        zonelist[i].port[1] = flowlist[-1].name
        zonelist[targetId].connect[0] = 1
        zonelist[targetId].port[0] = zonelist[i].port[1]
    if coords[1] < y - 1:
        flowlist.append(Flow(len(flowlist), 1, RetrieveId([coords[0], coords[1] + 1, coords[2]]), i, 2, 3))
        zonelist[i].connect[3] = 1
        zonelist[i].port[3] = flowlist[-1].name
        zonelist[RetrieveId([coords[0], coords[1] + 1, coords[2]])].connect[2] = 1
        zonelist[RetrieveId([coords[0], coords[1] + 1, coords[2]])].port[2] = zonelist[i].port[3]
    if coords[2] < z - 1:
        flowlist.append(Flow(len(flowlist), 2, RetrieveId([coords[0], coords[1], coords[2] + 1]), i, 4, 5))
        zonelist[i].connect[5] = 1
        zonelist[i].port[5] = flowlist[-1].name
        zonelist[RetrieveId([coords[0], coords[1], coords[2] + 1])].connect[4] = 1
        zonelist[RetrieveId([coords[0], coords[1], coords[2] + 1])].port[4] = zonelist[i].port[5]
for zone in zonelist:
    chain = ""
    for i in range(len(zone.connect)):
        chain += ", " + str(zone.connect[i])
    # print(chain)

walllist = [] # list of module Wall
for i in range(len(zonelist)):
    connectors = zonelist[i].connect
    for j in range(len(connectors)):
        if not connectors[j]:
            wallname = "wall" + str(len(walllist))
            walllist.append(Wall(len(walllist), j // 2, i, j))
            zonelist[i].port[j] = wallname

# locate the zone at the bottom and the top
adiabaticZonelist = [] # index list for top and bottom zones
for zone in zonelist:
    if zone.idx < x * y or zone.idx >= (z - 1) * x * y:
        adiabaticZonelist.append(zone.idx)

# module serialization
fo = open("Box.mo", "w")
fo.write("model Box\n")
if useActuralAmbient:
    fo.write("  Modelica.Blocks.Sources.TimeTable ambient(table = [{0}]);\n".format(ambientTemp))
for zone in zonelist:
    if zone.idx < x * y:
        fo.write("  Zone {0}(IsSource = true, dx = {1}, dy = {2}, dz = {3});\n".
            format(zone.name, dim_x/x, dim_y/y, dim_z/z))
        fo.write("  Modelica.Blocks.Sources.TimeTable solar{0}(table = [{1}]);\n".
            format(zone.idx, ZipTimeTable(times, solars[zone.idx])))
        fo.write("  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow heatFlow{0};\n".format(zone.idx))
    else:
        fo.write("  Zone {0}(dx = {1}, dy = {2}, dz = {3});\n".format(zone.name, dim_x/x, dim_y/y, dim_z/z))

for flow in flowlist:
    fo.write("  Flow {0}(Direction = {1});\n".format(flow.name, flow.direction))
for wall in walllist:
    if wall.idxIn in adiabaticZonelist and wall.direction == 2:
        fo.write("  Wall {0}(Direction = {1});\n".format(wall.name, wall.direction))
    else:
        fo.write("  Wall {0}(Direction = {1}, IsSource = true);\n".format(wall.name, wall.direction))
        fo.write("  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature temp{0};".format(wall.idx))
        fo.write("  Modelica.Thermal.HeatTransfer.Components.ThermalResistor wallR{0}(R = {1});".format(wall.idx, resistance))
        fo.write("  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor wallC{0}(C = {1});".format(wall.idx, capacity))

fo.write("equation\n")
for flow in flowlist:
    axis = "x"
    if flow.direction == 1: axis = "y"
    if flow.direction == 2: axis = "z"
    fo.write("  connect({0}.port_a, {1}.port_{2}1);\n".format(flow.name, "zone" + str(flow.idxIn), axis))
    fo.write("  connect({0}.port_b, {1}.port_{2}2);\n".format(flow.name, "zone" + str(flow.idxIn), axis))
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
    if zonelist.index(zone) < x * y:
        fo.write("  connect({0}.port_s, heatFlow{1}.port);\n".format(zone.name, zone.idx))
        fo.write("  connect(solar{0}.y, heatFlow{1}.Q_flow);\n".format(zone.idx, zone.idx))
    for i in range(6):
        fo.write("  connect({0}.i[{1}], {2}.o);\n".format(zone.name, i+1, zone.port[i]))
for wall in walllist:
    fo.write("  connect({0}.i, {1}.o);\n".format(wall.name, "zone" + str(wall.idxIn)))
    if not (wall.idxIn in adiabaticZonelist and wall.direction == 2):
        fo.write("  connect(temp{0}.T, ambient.y);\n".format(wall.idx))
        fo.write("  connect(temp{0}.port, wallR{1}.port_b);\n".format(wall.idx, wall.idx))
        fo.write("  connect(wallR{0}.port_a, wall{1}.port_s);\n".format(wall.idx, wall.idx))
        fo.write("  connect(wallC{0}.port, wallR{1}.port_a);\n".format(wall.idx, wall.idx))
        
        
fo.write("end Box;\n")
fo.close()