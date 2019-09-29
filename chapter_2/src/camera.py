import numpy as np

from src.geometry import Point, Rotation, Translation, Transform


class Pixel(Point):
    def __init(self, x, y, depth):
        super(Point, self).__init__()
        self.x = x
        self.y = y
        self.z = depth

    @property
    def depth(self):
        return self.z

    @property
    def uv(self):
        u = int(self.x / self.depth)
        v = int(self.y / self.depth)
        return u, v


class Distortion(object):
    pass


class IntrinsicParam(object):
    def __init__(self, fx, fy, cx, cy):
        self.__fx = fx
        self.__fy = fy
        self.__cx = cx
        self.__cy = cy
        self.__K = np.array([self.__fx, 0.0, self.__cx,
                             0.0, self.__fy, self.__cy,
                             0.0, 0.0, 1.0])
        self.__K = self.__K.reshape((3, 3))

    @property
    def K(self):
        return self.__K

    def __repr__(self):
        return "{}".format(self.K)


class Camera(object):
    def __init__(self, fx, fy, cx, cy, rotation, translation):
        self.__intrinsic_param = IntrinsicParam(fx, fy, cx, cy)
        self.__rotation = rotation
        self.__translation = translation
        self.__world_to_camera = Transform(self.__rotation, self.__translation)

    def project(self, point):
        """
        camera pinhole mode
        :param point: point in world coordinates
        :type point: Point
        :return: point in image coordinates
        :rtype: Pixel
        """
        # TODO: implement camera model
        ### your code starts ###
        M = self.__world_to_camera.homogeneous
        print(M)
        k = self.__intrinsic_param.K
        points_homo = point.to_homo_vector()
        Pc = np.dot(M, points_homo)
        InnerCo = np.dot(k, Pc[0:3, :])
        x = InnerCo[0]
        y = InnerCo[1]

        ### your code ends   ###
        return Pixel(x, y, 1)


def create_test_camera_model():
    # TODO: try different camera parameters
    fx = 300.0
    fy = 400.0
    cx = 450.0
    cy = 445.0

    rotation = Rotation()
    translation = Translation(0.0, 0.0, 0.0)
    return Camera(fx, fy, cx, cy, rotation, translation)
