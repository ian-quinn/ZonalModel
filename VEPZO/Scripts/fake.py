import math

import numpy as np
from matplotlib import pyplot as plt

def bubbleSort(arr):
	n = len(arr)
	for i in range(n):
		for j in range(0, n-i-1):
			if arr[j] > arr[j+1]:
				arr[j], arr[j+1] = arr[j+1], arr[j]

def FakeInput(epw_path, slice_num):
	# grab temperature seris from epw file
	series_temp = []		# contains all temperature series
	slices_temp = [[] for i in range(slice_num)]	# tempearture data of this slice only
	timeticks = 0
	f = open(epw_path, mode="r")
	for line in f.readlines():
		datalist = line.split(',')
		# the standard epw file has 35 columns each time step
		# and covers all 8760 hours of one typical year
		if len(datalist) == 35 and timeticks <= 8760:
			slices_temp[timeticks//int(8760/slice_num)]\
				.append(round(float(datalist[6]) + 273.15, 2))
			series_temp.append(round(float(datalist[6]) + 273.15, 2))
			timeticks = timeticks + 1

	amps = []
	for temps in slices_temp:
		bubbleSort(temps)
		amps.append(round(temps[-1] - temps[0], 1) / 2)

	print(amps)

	col_1 = []
	col_2 = []
	col_3 = []
	for i in range(8760):
		col_1.append(math.sin(2 * math.pi / 8760 * i))
		col_2.append(math.cos(2 * math.pi / 8760 * i))
		col_3.append(1)

	params = np.matrix([col_1, col_2, col_3]).T
	values = np.matrix([series_temp]).T
	fits = np.matmul(np.linalg.pinv(params), values)

	x = np.arange(0, 8760, 1)
	# y = fits[0, 0] * np.sin(2 * math.pi / 8760 * x) + \
	# 	fits[1, 0] * np.cos(2 * math.pi / 8760 * x) + \
	# 	fits[2, 0]
	# fig, ax = plt.subplots()
	# ax.plot(x, y, linewidth=2)
	# ax.scatter(x, series_temp, s=1)
	# plt.show()

	t = np.arange(0, 24 * slice_num, 1)
	temps = fits[0, 0] * np.sin(2 * math.pi / 24 / slice_num * t) + \
		fits[1, 0] * np.cos(2 * math.pi / 24 / slice_num * t) + fits[2, 0]
	overlay = [0 for i in range(len(t))]
	for i in range(len(t)):
		overlay[i] += -amps[i//24] * np.sin(2 * math.pi / 24 * (i%24))
	for i in range(len(t)):
		if i + 6 >= len(t):
			temps[i] += overlay[i + 6 - len(t)]
		else:
			temps[i] += overlay[i + 6]

	return temps, series_temp

if __name__ =='__main__':

	slice_num = 24
	fake_temp, actual_temp = FakeInput("Shanghai.epw", slice_num)

	x = np.arange(0, 8760, 1)
	t = np.arange(0, 24 * slice_num, 1)
	fig, ax = plt.subplots()
	ax.plot(t, fake_temp, linewidth=2)
	ax.scatter(x * 24 * slice_num / 8760, actual_temp, s=1, c="#88c999")
	plt.show()