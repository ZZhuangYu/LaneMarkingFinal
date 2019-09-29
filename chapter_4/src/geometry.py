import numpy as np


def convert_to_homogeneous(rotation, translation):
    return np.block([[rotation.rotation, translation.to_vector()],
                     [0.0, 0.0, 0.0, 1.0]])


def rotation_to_quaternion(r):
    return Quaternion()


def rotation_to_rpy(r):
    return 0.0, 0.0, 0.0


def rpy_to_rotation(r, p, y):
    yaw_matrix = np.matrix([
        [np.cos(y), -np.sin(y), 0.0],
        [np.sin(y), np.cos(y), 0.0],
        [0.0, 0.0, 1.0]
    ])
    pitch_matrix = np.matrix([
        [np.cos(p), 0.0, np.sin(p)],
        [0.0, 1.0, 0.0],
        [-np.sin(p), 0.0, np.cos(p)]
    ])
    roll_matrix = np.matrix([
        [1.0, 0.0, 0.0],
        [0.0, np.cos(r), -np.sin(r)],
        [0.0, np.sin(r), np.cos(r)]
    ])
    rotation = yaw_matrix * pitch_matrix * roll_matrix
    return rotation


def rpy_to_quaternion(r, p, y):
    return Quaternion()


def quaternion_to_rpy(q):
    return 0.0, 0.0, 0.0


def quaternion_to_rotation(q):
    return np.eye(3)


def is_rotation_valid(r):
    return True


class Point(object):
    def __init__(self, x, y, z):
        self.__x = x
        self.__y = y
        self.__z = z

    @property
    def x(self):
        return self.__x

    @x.setter
    def x(self, x):
        self.__x = x

    @property
    def y(self):
        return self.__y

    @y.setter
    def y(self, y):
        self.__y = y

    @property
    def z(self):
        return self.__z

    @z.setter
    def z(self, z):
        self.__z = z

    def to_vector(self):
        vec = np.array([self.x, self.y, self.z])
        vec = vec.reshape((3, 1))
        return vec

    def to_homo_vector(self):
        vec = np.array([self.x, self.y, self.z, 1.0])
        vec = vec.reshape((4, 1))
        return vec

    def __repr__(self):
        return "Point({}, {}, {})".format(self.x, self.y, self.z)


class Translation(Point):
    def __init__(self, x, y, z):
        super(Point, self).__init__()
        self.x = x
        self.y = y
        self.z = z

    def __repr__(self):
        return "Translation({}, {}, {})".format(self.x, self.y, self.z)

    def apply_to_point(self, point):
        translated = point.to_vector() + self.to_vector()
        return Point(translated[0, 0], translated[1, 0], translated[2, 0])


class Quaternion(object):
    def __init__(self, w=1.0, x=0.0, y=0.0, z=0.0):
        self.__w = w
        self.__x = x
        self.__y = y
        self.__z = z

    def __repr__(self):
        return "Quaternion({}, {}, {}, {})".format(self.__w, self.__x, self.__y, self.__z)


class Rotation(object):
    def __init__(self):
        self.__roll = 0.0
        self.__pitch = 0.0
        self.__yaw = 0.0
        self.__rotation = np.eye(3)
        self.__quaternion = Quaternion()

    @staticmethod
    def from_rpy(roll=0.0, pitch=0.0, yaw=0.0):
        rotation = Rotation()
        rotation.__roll = roll
        rotation.__pitch = pitch
        rotation.__yaw = yaw
        rotation.__rotation = rpy_to_rotation(roll, pitch, yaw)
        rotation.__quaternion = rpy_to_quaternion(roll, pitch, yaw)
        return rotation

    @staticmethod
    def from_rotation_matrix(r=np.eye(3)):
        if not is_rotation_valid(r):
            raise Exception("Invalid rotation matrix")
        rotation = Rotation()
        rotation.__rotation = r
        rotation.__roll, rotation.__pitch, rotation.__yaw = rotation_to_quaternion(r)
        return rotation

    @staticmethod
    def from_quaternion(q=Quaternion()):
        rotation = Rotation()
        rotation.__quaternion = q
        rotation.__roll, rotation.__pitch, rotation.__yaw = quaternion_to_rpy(q)
        rotation.__rotation = quaternion_to_rotation(q)
        return rotation

    @property
    def rotation(self):
        return self.__rotation

    @rotation.setter
    def rotation(self, r):
        self.__rotation = r
        self.__roll, self.__pitch, self.__yaw = rotation_to_rpy(r)
        self.__quaternion = rotation_to_quaternion(r)
    
    @property
    def rpy(self):
        return self.__roll, self.__pitch, self.__yaw

    @rpy.setter
    def rpy(self, r, p, y):
        self.__roll = r
        self.__pitch = p
        self.__yaw = y
        self.__rotation = rpy_to_rotation(r, p, y)
        self.__quaternion = rpy_to_quaternion(r, p, y)

    @property
    def quaternion(self):
        return self.__quaternion

    @quaternion.setter
    def quaternion(self, q):
        self.__quaternion = q
        self.__roll, self.__pitch, self.__yaw = quaternion_to_rpy(q)
        self.__rotation = quaternion_to_rotation(q)

    def apply_to_point(self, point):
        rotated = np.dot(self.rotation, point.to_vector())
        return Point(rotated[0, 0], rotated[1, 0], rotated[2, 0])

    def __repr__(self):
        text = "Rotation:\n"
        text += "{}".format(self.rotation)
        text += "\nRPY: {}".format(self.rpy)
        text += "\n{}".format(self.quaternion)
        return text


class Transform(object):
    def __init__(self, rotation, translation):
        self.rotation = rotation
        self.translation = translation
        self.homogeneous = convert_to_homogeneous(self.rotation, self.translation)

    def apply_to_point(self, point):
        rotated = self.rotation.apply_to_point(point)
        translated = self.translation.apply_to_point(rotated)
        return translated

    def __repr__(self):
        return "{}\n{}".format(self.rotation, self.translation)


if __name__ == "__main__":
    # point
    p = Point(1.0, 2.0, 3.0)
    print(p)
    # translation
    trans = Translation(2, 3, 4)
    print(trans)
    p_trans = trans.apply_to_point(p)
    print(p_trans)
    # rotation
    r = Rotation.from_rpy(0.0, 1.0, 2.0)
    print(r)
    p_rotated = r.apply_to_point(p)
    print(p_rotated)
    transform = Transform(r, trans)
    print(transform)
    transformed = transform.apply_to_point(p)
    print(transformed)
