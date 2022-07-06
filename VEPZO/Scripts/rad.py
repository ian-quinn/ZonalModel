import math
import numpy as np

def AssignProduct(arr, mask):
	result = [i for i in mask]
	idx = 0
	for i in range(len(mask)):
		if mask[i] != 0:
			result[i] = arr[idx]
			idx += 1
		else:
			result[i] = 0
	return np.array(result)

def GetViewFactor(basepts, normals, dims, d):
	view_matrix = np.zeros((len(basepts), len(basepts)))
	for i in range(len(basepts)):
		revnormal_A = 1 - abs(normals[i])
		dims_A = revnormal_A * dims
		basearea = dims_A[dims_A != 0][0] * dims_A[dims_A != 0][1]
		for j in range(len(basepts)):
			print("Pairing {0}-{1}".format(i, j))
			# only consider 2d rectangle for now
			# the self-viewed factor must be 0
			if i == j:
				view_matrix[i, j] = 0
			else:
				# discretization process
				base_A = basepts[i]
				base_B = basepts[j]
				revnormal_B = 1 - abs(normals[j])
				mesh_A = []
				mesh_B = []

				u = 0
				vec = dims_A[dims_A != 0]
				while u < vec[0]:
					v = 0
					while v < vec[1]:
						mesh_A.append(base_A + 0.5 * d * revnormal_A + \
							AssignProduct(np.array([u, v]), revnormal_A))
						v += d
					u += d

				u = 0
				dims_B = revnormal_B * dims
				vec = dims_B[dims_B != 0]
				while u < vec[0]:
					v = 0
					while v < vec[1]:
						mesh_B.append(base_B + 0.5 * d * revnormal_B + \
							AssignProduct(np.array([u, v]), revnormal_B))
						v += d
					u += d

				# if i == 0:
				# 	print(mesh_A)
				# 	print(mesh_B)

				# if i < 2:
				#	 for pt in mesh_B:
				#		 print("({0}, {1}, {2})".format(pt.x, pt.y, pt.z))
				summation = 0
				for ptA in mesh_A:
					for ptB in mesh_B:
						view_d = np.linalg.norm(ptB - ptA)
						if view_d < 0.000001:
							continue
						cos_theta_i = sum(normals[i] * (ptB - ptA))
						cos_theta_j = sum(normals[j] * (ptA - ptB))
						summation += cos_theta_i * cos_theta_j * d ** 4 / math.pi / view_d ** 4
						if i == 0:
							print(normals[i], normals[j])
							print(ptA, ptB)
							print(cos_theta_i, cos_theta_j)
				view_matrix[i, j] = summation / basearea
	return view_matrix


# if __name__ =='__main__':