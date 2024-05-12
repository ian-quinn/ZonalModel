import math
import time
import numpy as np


def assign_product(arr, mask):
    result = [i for i in mask]
    idx = 0
    for i in range(len(mask)):
        if mask[i] != 0:
            result[i] = arr[idx]
            idx += 1
        else:
            result[i] = 0
    return np.array(result)


def area_on_normal_direction(dims, normal):
    normal_rev = 1 - abs(normal)
    area = 1
    for coord in dims * normal_rev:
        if coord != 0:
            area = area * coord
    return area


def get_view_factor(locations, normals, dims, size: int):
    """
    :param locations: list of np.array representing the location of wall surface (closest point to (0, 0, 0))
    :param normals: list of np.array representing the face normal
    :param dims: np.array for the size of the zone cell (all zones are the same)
    :param size: number of mesh along the edge. large number will increase the calculation time
    :return: the matrix of view factors with dimension len(locations)^2
    """
    start_time = time.time()
    surface_num = len(locations)
    view_matrix = np.zeros((surface_num, surface_num))
    for i in range(surface_num):
        normal_rev_a = 1 - abs(normals[i])
        dims_a = normal_rev_a * dims
        area_a = dims_a[dims_a != 0][0] * dims_a[dims_a != 0][1]
        for j in range(surface_num):
            # print("Pairing {0}-{1}".format(i, j))
            # only consider 2d rectangle for now
            if i == j:
                view_matrix[i, j] = 0  # the self-viewed factor must be 0
            else:
                # discretization process
                origin_a = locations[i]
                origin_b = locations[j]
                normal_rev_b = 1 - abs(normals[j])
                mesh_a = []
                mesh_b = []

                u = 0
                vec = dims_a[dims_a != 0]  # take dimensions that are not 0, XYZ -> UV system
                while u < size:
                    v = 0
                    while v < size:
                        mesh_a.append(origin_a + 0.5 * dims_a / size +
                                      assign_product(np.array([u * vec[0] / size, v * vec[1] / size]), normal_rev_a))
                        v += 1
                    u += 1

                u = 0
                dims_b = normal_rev_b * dims
                vec = dims_b[dims_b != 0]
                while u < size:
                    v = 0
                    while v < size:
                        mesh_b.append(origin_b + 0.5 * dims_b / size +
                                      assign_product(np.array([u * vec[0] / size, v * vec[1] / size]), normal_rev_b))
                        v += 1
                    u += 1

                summation = 0
                for ptA in mesh_a:
                    for ptB in mesh_b:
                        view_d = np.linalg.norm(ptB - ptA)
                        if view_d < 0.000001:
                            continue
                        cos_theta_i = sum(normals[i] * (ptB - ptA)) / view_d
                        cos_theta_j = sum(normals[j] * (ptA - ptB)) / view_d
                        # if the grid system is 3*3, then the mesh area is 1/9 of the base area
                        # the z direction sizing corresponds to the dim_x / size
                        # thus we can take this ratio as 1 / (dims[0]/size)^2
                        summation += (cos_theta_i * cos_theta_j *
                                      (area_on_normal_direction(dims, normals[i]) / size ** 2) *
                                      (area_on_normal_direction(dims, normals[j]) / size ** 2) /
                                      math.pi / view_d ** 2)
                view_matrix[i, j] = summation / area_a

    end_time = time.time()
    print("Radiation box view factor matrix calculation time: ", end_time - start_time, " s")

    # a simple box for example:
    # the view factor between two surfaces facing each other: 0.19982489569838746
    # the view factor between two adjacent surfaces: 0.20004377607540316

    return view_matrix


def get_mrt_factor(locations, normals, dims, mask, size: int, height=0):
    """
    :param locations: list of np.array representing the base point of wall (closest point to origin)
    :param normals: list of np.array representing the face normal
    :param dims: np.array for the size of the zone cell (all zones are the same)
    :param mask: 2D np.array representing the floorplan mesh that will be the MRT contour
    :param size: number of mesh grid along the edge. large number will increase the calculation time
    :param height: operative height to measure the MRT. iteration through all zones by default
    :return: the matrix of view factors with dimension (num_bottom_faces * num_other_faces)
    """
    start_time = time.time()
    surface_num = len(locations)
    # the mesh cell iteration starts from the left-bottom corner
    # first go row then go column
    # its sequence should be in line with the zone
    mask_row = np.size(mask, 0)
    mask_col = np.size(mask, 1)
    mesh_num = 0
    mesh_centroids = []
    for z in range(3):
        for i in range(mask_row - 1, -1, -1):
            for j in range(mask_col):
                if mask[i, j] == 1:
                    mesh_num += 1
                    z_coord = dims[2] * z + dims[2] / 2
                    if height >= 0:
                        z_coord = height
                    else:
                        z_coord = 0
                    mesh_coord = [j * dims[1] + dims[1] / 2, (mask_row - i - 1) * dims[0] + dims[0] / 2, z_coord]
                    mesh_centroids.append(np.array(mesh_coord))
        if height >= 0:
            break

    # each mesh cell will reach out to all faces
    view_matrix = np.zeros((mesh_num, surface_num))

    # iterate through each MRT mesh grid
    for i in range(mesh_num):

        # iterate through each wall face
        for j in range(surface_num):

            normal_rev_b = 1 - abs(normals[j])
            dims_b = normal_rev_b * dims

            # discretize the face A
            origin_b = locations[j]
            mesh_b = []

            # ptA is the observation point, the differential area
            ptA = mesh_centroids[i]

            u = 0
            vec = dims_b[dims_b != 0]
            while u < size:
                v = 0
                while v < size:
                    mesh_b.append(origin_b + 0.5 * dims_b / size +
                                  assign_product(np.array([u * vec[0] / size, v * vec[1] / size]), normal_rev_b))
                    v += 1
                u += 1

            summation = 0
            for ptB in mesh_b:
                view_d = np.linalg.norm(ptB - ptA)
                if view_d < 0.000001:
                    continue
                # normal of A is unit Z constantly
                cos_theta_i = sum(np.array([0, 0, 1]) * (ptB - ptA)) / view_d
                cos_theta_j = sum(normals[j] * (ptA - ptB)) / view_d
                summation += (cos_theta_i * cos_theta_j *
                              area_on_normal_direction(dims, normals[i]) / size ** 2 /
                              math.pi / view_d ** 2)
            view_matrix[i, j] = summation

    end_time = time.time()
    print("MRT solution matrix calculation time: ", end_time - start_time, " s")

    # a simple box for example:
    # the view factor from the centroid of the bottom to the vertical wall: 0.19013588238480664
    # the view factor from the centroid of the bottom to the top face: 0.239456470460773536

    return view_matrix


# test this module
if __name__ == '__main__':
    """
    cube, each face 3x3m
    0,0,0
    0,1,0
    0,0,0
    """
    mesh_str = "0,0,0,0,1,0,0,0,0"
    # sequence: bottom, top, left, right, front, back
    basepts = [np.array([3., 3., 0.]), np.array([3., 3., 3.]), np.array([3., 3., 0.]),
               np.array([6., 3., 0.]), np.array([3., 3., 0.]), np.array([3., 6., 0.])]
    normals = [np.array([0, 0, 1]), np.array([0, 0, -1]), np.array([1, 0, 0]),
               np.array([-1, 0, 0]), np.array([0, 1, 0]), np.array([0, -1, 0])]
    dims = np.array([3, 3, 3])
    mesh = np.fromstring(mesh_str, sep=',').reshape(3, 3)
    factors = get_mrt_factor(basepts, normals, dims, mesh, 10, 0)
    print(factors)
