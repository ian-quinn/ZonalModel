import math
import time
import numpy as np


# #################### Geometry transformations #########################

def trace_vec(vec):
    return vec[0] + vec[1] + vec[2]


def unit_vec(vec):
    return vec / np.linalg.norm(vec)


def angle_between(v1, v2):
    v1_u = unit_vec(v1)
    v2_u = unit_vec(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


def angle_from(v1, v2):
    angle_radians = angle_between(v1, v2)
    if np.cross(v1, v2)[2] < 0:
        angle_radians = 2 * np.pi - angle_radians
    return angle_radians


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


def rotate_point(pt, axis: str, angle):
    rad = angle / 180 * math.pi
    if axis == 'x':
        rotation_matrix = np.array([
            [1, 0, 0],
            [0, np.cos(rad), -np.sin(rad)],
            [0, np.sin(rad), np.cos(rad)]
        ])
    elif axis == 'y':
        rotation_matrix = np.array([
            [np.cos(rad), 0, np.sin(rad)],
            [0, 1, 0],
            [-np.sin(rad), 0, np.cos(rad)]
        ])
    elif axis == 'z':
        rotation_matrix = np.array([
            [np.cos(rad), -np.sin(rad), 0],
            [np.sin(rad), np.cos(rad), 0],
            [0, 0, 1]
        ])

    # Perform the rotation
    rotated_pt = np.dot(rotation_matrix, pt)

    return np.round(rotated_pt, decimals=2)


def rotate_vector(vec, axis, angle):
    rad = angle / 180 * math.pi
    cross_product_matrix = np.array([
        [0, -axis[2], axis[1]],
        [axis[2], 0, -axis[0]],
        [-axis[1], axis[0], 0]
    ])
    rotation_matrix = (np.eye(3) * np.cos(rad) + np.sin(rad) * cross_product_matrix +
                       (1 - np.cos(rad)) * np.outer(axis, axis))

    # Perform the rotation
    rotated_vec = np.dot(rotation_matrix, vec)

    return np.round(rotated_vec, decimals=2)


def replicate_face(origin, normal, dims):
    vertices = [origin]
    if normal[2] == 1 or normal[2] == -1:
        vertices.append(np.array([origin[0] + dims[0], origin[1], origin[2]]))
        vertices.append(np.array([origin[0] + dims[0], origin[1] + dims[1], origin[2]]))
        vertices.append(np.array([origin[0], origin[1] + dims[1], origin[2]]))
    if normal[1] == 1 or normal[1] == -1:
        vertices.append(np.array([origin[0] + dims[0], origin[1], origin[2]]))
        vertices.append(np.array([origin[0] + dims[0], origin[1], origin[2] + dims[2]]))
        vertices.append(np.array([origin[0], origin[1], origin[2] + dims[2]]))
    if normal[0] == 1 or normal[0] == -1:
        vertices.append(np.array([origin[0], origin[1] + dims[1], origin[2]]))
        vertices.append(np.array([origin[0], origin[1] + dims[1], origin[2] + dims[2]]))
        vertices.append(np.array([origin[0], origin[1], origin[2] + dims[2]]))
    return vertices


def rotate_face(face, axis, angle):
    face_rotated = []
    for vertex in face:
        face_rotated.append(rotate_point(vertex, axis, angle))
    return face_rotated


def locate_bounding_box_pts(pts):
    min_x = 10000.0
    min_y = 10000.0
    min_z = 10000.0
    max_x = -10000.0
    max_y = -10000.0
    max_z = -10000.0
    for pt in pts:
        if pt[0] < min_x: min_x = pt[0]
        if pt[1] < min_y: min_y = pt[1]
        if pt[2] < min_z: min_z = pt[2]
        if pt[0] > max_x: max_x = pt[0]
        if pt[1] > max_y: max_y = pt[1]
        if pt[2] > max_z: max_z = pt[2]
    return [[min_x, max_x], [min_y, max_y], [min_z, max_z]]


def area_on_normal_direction(dims, normal):
    normal_rev = 1 - abs(normal)
    area = 1
    for coord in dims * normal_rev:
        if coord != 0:
            area = area * coord
    return area


# #################### VF formulas between rectangles ##################


def get_vertical_rectangle_vf(x: tuple[float, float], y: tuple[float, float],
                              xi: tuple[float, float], eta: tuple[float, float]) -> float:
    def calc_b(x, y, xi, eta):
        c = math.sqrt(pow(x, 2) + pow(xi, 2))
        d = (y - eta) / c
        b = ((y - eta) * c * math.atan(d) - pow(c, 2) / 4 * (1 - pow(d, 2))
             * math.log(pow(c, 2) * (1 + pow(d, 2))))
        return b

    # safe lock
    if x[0] == 0: x = (0.000001, x[1])
    if xi[0] == 0: xi = (0.000001, xi[1])

    sigma = 0
    for l in range(2):
        for k in range(2):
            for j in range(2):
                for i in range(2):
                    sigma += pow(-1, i + j + k + l) * calc_b(x[i], y[j], xi[k], eta[l])
    return sigma / 2 / math.pi / (x[1] - x[0]) / (y[1] - y[0])


def get_parallel_rectangle_vf(x: tuple[float, float], y: tuple[float, float],
                              xi: tuple[float, float], eta: tuple[float, float],
                              z: float) -> float:
    def calc_b(x, y, xi, eta, z):
        u = x - xi
        v = y - eta
        p = math.sqrt(pow(u, 2) + pow(z, 2))
        q = math.sqrt(pow(v, 2) + pow(z, 2))
        b = (v * p * math.atan(v / p) + u * q * math.atan(u / q) -
             pow(z, 2) / 2 * math.log(pow(u, 2) + pow(v, 2) + pow(z, 2)))
        return b

    sigma = 0
    for l in range(2):
        for k in range(2):
            for j in range(2):
                for i in range(2):
                    sigma += pow(-1, i + j + k + l) * calc_b(x[i], y[j], xi[k], eta[l], z)
    return sigma / 2 / math.pi / (x[1] - x[0]) / (y[1] - y[0])


def switch_dimension(origin, normal, dims):
    # with limited usage
    tension = dims * (1 - abs(normal))
    terminal = origin + tension
    # list for tuples representing the value range of each dimension, max=2
    ranges = []
    for i in range(3):
        if normal[i] == 0:
            ranges.append((origin[i], terminal[i]))
    return ranges[0], ranges[1]


# #################### Main operations ####################


def get_view_factor(locations, normals, dims, size: int):
    """
    The numerical way for view factor calculation (update for any shape in future)
    :param locations: list of np.array representing the location of wall surface (closest point to (0, 0, 0))
    :param normals: list of np.array representing the face normal
    :param dims: np.array for the size of the zone cell (all zones are the same)
    :param size: number of mesh along the edge. large number will increase the calculation time
    :return: the matrix of view factors with dimension len(locations)^2 (numerical method)
    """
    start_time = time.time()
    surface_num = len(locations)
    view_matrix = np.zeros((surface_num, surface_num))
    for i in range(surface_num):
        normal_rev_a = 1 - abs(normals[i])
        dims_a = normal_rev_a * dims
        area_a = dims_a[dims_a != 0][0] * dims_a[dims_a != 0][1]
        for j in range(surface_num):
            # note that view_matrix is symmetric. copy its value if view_matrix[j,i] exists
            if view_matrix[j, i] > 0:
                view_matrix[i, j] = view_matrix[j, i]
                continue
            # if two faces are on the same plane
            if np.array_equal(normals[i], normals[j]):
                view_matrix[i, j] = 0
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


def calc_view_factor(locations, normals, dims):
    """
    The analytical way for view factor calculation on parallel or vertical rectangles
    :param locations: list of np.array representing the location of wall surface (closest point to (0, 0, 0))
    :param normals: list of np.array representing the face normal
    :param dims: np.array for the size of the zone cell (all zones are the same)
    :param size: number of mesh along the edge. large number will increase the calculation time
    :return: the matrix of view factors with dimension len(locations)^2 (analysis method)
    """
    start_time = time.time()
    surface_num = len(locations)
    view_matrix = np.zeros((surface_num, surface_num))
    for i in range(surface_num):
        for j in range(surface_num):
            # note that view_matrix is symmetric. copy its value if view_matrix[j,i] exists
            if view_matrix[j, i] > 0:
                view_matrix[i, j] = view_matrix[j, i]
                continue
            # if two faces are on the same plane
            if np.array_equal(normals[i], normals[j]):
                view_matrix[i, j] = 0
            else:
                # if two faces are parallel and frontal
                if np.dot(normals[i], normals[j]) != 0:
                    z = np.linalg.norm(locations[i] * normals[i] + locations[j] * normals[j])
                    x, y = switch_dimension(locations[i], normals[i], dims)
                    u, v = switch_dimension(locations[j], normals[j], dims)
                    view_matrix[i, j] = get_parallel_rectangle_vf(x, y, u, v, z)
                    print(f"{i}-{j} para: {view_matrix[i, j]}")
                if np.dot(normals[i], normals[j]) == 0:
                    a_face = replicate_face(locations[i], normals[i], dims)
                    b_face = replicate_face(locations[j], normals[j], dims)
                    a_vec = normals[i]
                    b_vec = normals[j]

                    hinge_axis = np.cross(normals[i], normals[j])
                    if abs(hinge_axis)[0] > 0:
                        a_face = rotate_face(a_face, "z", -90)
                        b_face = rotate_face(b_face, "z", -90)
                        a_vec = rotate_vector(a_vec, np.array([0, 0, 1]), -90)
                        b_vec = rotate_vector(b_vec, np.array([0, 0, 1]), -90)
                    if abs(hinge_axis)[2] > 0:
                        a_face = rotate_face(a_face, "x", -270)
                        b_face = rotate_face(b_face, "x", -270)
                        a_vec = rotate_vector(a_vec, np.array([1, 0, 0]), -270)
                        b_vec = rotate_vector(b_vec, np.array([1, 0, 0]), -270)
                    while trace_vec(a_vec) < 0 or trace_vec(b_vec) < 0:
                        a_face = rotate_face(a_face, "y", 90)
                        b_face = rotate_face(b_face, "y", 90)
                        a_vec = rotate_vector(a_vec, np.array([0, 1, 0]), 90)
                        b_vec = rotate_vector(b_vec, np.array([0, 1, 0]), 90)
                    # move faces based on box corner to global origin
                    box = locate_bounding_box_pts(a_face + b_face)
                    origin = np.array([box[0][0], box[1][0], box[2][0]])
                    a_face_ = [a_face[i] - origin for i in range(4)]
                    b_face_ = [b_face[i] - origin for i in range(4)]
                    a_box = locate_bounding_box_pts(a_face_)
                    b_box = locate_bounding_box_pts(b_face_)
                    a_range = []
                    b_range = []
                    if a_vec[0] == 1:
                        a_range = [a_box[2], a_box[1]]
                        b_range = [b_box[0], b_box[1]]
                    else:
                        a_range = [a_box[0], a_box[1]]
                        b_range = [b_box[2], b_box[1]]
                    x = a_range[0]
                    y = a_range[1]
                    u = b_range[0]
                    v = b_range[1]

                    view_matrix[i, j] = get_vertical_rectangle_vf(x, y, u, v)
                    print(f"{i}-{j} vert: {view_matrix[i, j]}")

    end_time = time.time()
    print("Radiation box view factor matrix calculation time: ", end_time - start_time, " s")

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

    # test case 1x1 cube with 6 faces
    # basepts = [np.array([3., 3., 0.]), np.array([3., 3., 3.]), np.array([3., 3., 0.]),
    #            np.array([6., 3., 0.]), np.array([3., 3., 0.]), np.array([3., 6., 0.])]
    # normals = [np.array([0, 0, 1]), np.array([0, 0, -1]), np.array([1, 0, 0]),
    #            np.array([-1, 0, 0]), np.array([0, 1, 0]), np.array([0, -1, 0])]

    # test case 2x2 cube with 24 faces
    basepts = [np.array([3., 3., 0.]), np.array([3., 3., 0.]), np.array([3., 3., 0.]),
               np.array([9., 3., 0.]), np.array([6., 3., 0.]), np.array([6., 3., 0.]),
               np.array([3., 6., 0.]), np.array([3., 9., 0.]), np.array([3., 6., 0.]),
               np.array([9., 6., 0.]), np.array([6., 9., 0.]), np.array([6., 6., 0.]),
               np.array([3., 3., 3.]), np.array([3., 3., 3.]), np.array([3., 3., 6.]),
               np.array([9., 3., 3.]), np.array([6., 3., 3.]), np.array([6., 3., 6.]),
               np.array([3., 6., 3.]), np.array([3., 9., 3.]), np.array([3., 6., 6.]),
               np.array([9., 6., 3.]), np.array([6., 9., 3.]), np.array([6., 6., 6.])]

    normals = [np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1]),
               np.array([-1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1]),
               np.array([1, 0, 0]), np.array([0, -1, 0]), np.array([0, 0, 1]),
               np.array([-1, 0, 0]), np.array([0, -1, 0]), np.array([0, 0, 1]),
               np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, -1]),
               np.array([-1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, -1]),
               np.array([1, 0, 0]), np.array([0, -1, 0]), np.array([0, 0, -1]),
               np.array([-1, 0, 0]), np.array([0, -1, 0]), np.array([0, 0, -1])]

    dims = np.array([3, 3, 3])
    mesh = np.fromstring(mesh_str, sep=',').reshape(3, 3)
    # factors = get_mrt_factor(basepts, normals, dims, mesh, 10, 0)
    factors = calc_view_factor(basepts, normals, dims)
    print(factors)
