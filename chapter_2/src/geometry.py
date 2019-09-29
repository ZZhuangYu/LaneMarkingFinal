import numpy as np


def convert_to_homogeneous(rotation, translation):
    return np.block([[rotation.rotation, translation.to_vector()],
                     [0.0, 0.0, 0.0, 1.0]])


def rotation_to_quaternion(r):
    # TODO: rotation matrix to quaternion
    ### your code starts ###
    w = 1/2 * (np.sqrt(1+r[0][0]+r[1][1])+r[2][2])
    x = (r[2][1]-r[1][2])/4 * w
    y = (r[0][2]-r[2][0])/4 * w
    z = (r[1][0]-r[0][1])/4 * w
    ### your code ends   ###
    return Quaternion(w, x, y, z)


def rotation_to_rpy(r):
    # TODO: rotation matrix to euler angles
    ### your code starts ###
    # Base on the rotation matrix below, we could solve Euler Angle 'r', 'p'  'y' from a given rotation matrix.
    # r = arctan(r32/r33)
    # p = arcsin(-r31)
    # y = arctan(r21/r11)

    r11 = r[0][0]
    r21 = r[1][0]
    r31 = r[2][0]
    r32 = r[2][1]
    r33 = r[2][2]
    r = (np.arctan2(r32, r33))*180/np.pi
    p = np.arcsin(-r31)*180/np.pi
    y = np.arctan2(r21, r11)*180/np.pi
    ### your code ends   ###
    return r, p, y


def rpy_to_rotation(r, p, y):
    # TODO: euler angles to rotation matrix
    ### your code starts ###
    # We assume that the rotation is executed in order: First rotate angle 'r' in direction 'x',
    # then 'p' in direction 'y', and finally 'y' in direction 'z',
    # which determines the order of matrices of each direction as: Rz * Ry * Rx

    # To be more specific, the matrix of their product would be:

    # [ cos(y)cos(z)       -cos(x)sin(z)+sin(x)sin(y)cos(z)       sin(x)sin(z)+cos(x)sin(y)cos(z)
    #   cos(y)sin(z)       cos(x)cos(z)+sin(x)sin(y)sin(z)        -sin(x)cos(z)+cos(x)sin(y)sin(z)
    #   -sin(y)            sin(x)cos(y)                           cos(x)cos(y)                     ]

    #In the following program 'x','y','z' would be replaced by 'r','p','y', where they have the same meanings.
    #By setting each rotation angle to 0, we could check the correctness of the matrix above.
    Cx = np.cos(r)
    Cy = np.cos(p)
    Cz = np.cos(y)
    Sx = np.sin(r)
    Sy = np.sin(p)
    Sz = np.sin(y)

    r11 = Cz*Cy
    r12 = Sx*Sy*Cz - Cx*Sz
    r13 = Sx*Sz + Cx*Sy*Cz
    r21 = Cy*Sz
    r22 = Cx*Cz + Sx*Sy*Sz
    r23 = Cx*Sy*Sz - Sx*Cz
    r31 = -Sy
    r32 = Sx*Sy
    r33 = Cx*Cy
    rotation_mat = np.array([[r11, r12, r13], [r21, r22, r23], [r31, r32, r33]])
    ### your code ends   ###
    return rotation_mat


def rpy_to_quaternion(r, p, y):
    # TODO: euler angles to quaternion
    ### your code starts ###
    w = np.cos(p/2)*np.cos(y/2)*np.cos(r/2)-np.sin(p/2)*np.sin(y/2)*np.sin(r/2)
    x = np.sin(p/2)*np.sin(y/2)*np.cos(r/2)+np.cos(p/2)*np.cos(y/2)*np.sin(r/2)
    y = np.sin(p/2)*np.cos(y/2)*np.cos(r/2)+np.cos(p/2)*np.sin(y/2)*np.sin(r/2)
    z = np.cos(p/2)*np.sin(y/2)*np.sin(y/2)-np.sin(p/2)*np.cos(y/2)*np.sin(r/2)
    ### your code ends   ###
    return Quaternion(w, x, y, z)


def quaternion_to_rpy(q):
    # TODO: quaternion to euler angles
    ### your code starts ###
    wq = q[0]
    xq = q[1]
    yq = q[2]
    zq = q[3]
    y = np.arctan2(2*(xq*yq-wq*zq), (pow(wq, 2)+pow(xq, 2)-pow(yq, 2)-pow(zq, 2)))
    r = np.arctan2(2*(yq*zq-wq*xq), (pow(wq, 2)-pow(xq, 2)-pow(yq, 2)-pow(zq, 2)))
    p = np.arcsin(-2*(wq*yq + xq*zq))
    ### your code ends   ###
    return r, p, y


def quaternion_to_rotation(q):
    # TODO: quaternion to ration_matrix
    ### your code starts ###
    w = q[0]
    x = q[1]
    y = q[2]
    z = q[3]
    r = np.array([[pow(w, 2)+pow(x, 2)-pow(y, 2)-pow(z, 2), 2*(w*z-x*y), 2*(x*z-w*y)],
                  [2*(x*y-w*z), pow(w, 2)-pow(x, 2)+pow(y, 2)-pow(z, 2), 2*(w*x+y*z)],
                  [2*(w*y+x*z), 2*(y*z-w*x), pow(w, 2)-pow(x, 2)-pow(y, 2)+pow(z, 2)]])

    ### your code ends   ###
    return r


def is_rotation_valid(r):
    # TODO: validate rotation matrix
    ### your code starts ###
    #Rotation Matrix requires its columns are orthogonal, inner product equals 1 and
    #determinant equals 1
    #We assume the matrix inputted is size of 3*3.
    r12 = np.dot(r[:, 0], r[:, 1])
    r23 = np.dot(r[:, 1], r[:, 2])
    r31 = np.dot(r[:, 2], r[:, 0])
    det = np.linalg.det(r)

    if r12 == 0 and r23 == 0 and r31 == 0 and det == 1:
    ### your code ends   ###
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
        return rotation

    @staticmethod
    def from_rotation_matrix(r=np.eye(3)):
        if not is_rotation_valid(r):
            raise Exception("Invalid rotation matrix")
        rotation = Rotation()
        rotation.__rotation = r
        return rotation

    @staticmethod
    def from_quaternion(q=Quaternion()):
        rotation = Rotation()
        rotation.__quaternion = q
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
