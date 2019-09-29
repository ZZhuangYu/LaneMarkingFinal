import numpy as np

from src.utils import Pixel
from src.geometry import Point


class ParametricPolyCurve(object):
    """
    Three order parametric polynomial curve fitting
    x = cx0 + cx1 * s + cx2 * s^2 + cx3 * s^3
    y = cy0 + cy1 * s + cy2 * s^2 + cy3 * s^3
    """
    def __init__(self, points, weights):
        if len(points) < 4:
            raise Exception("not enough points for fitting")
        self.points = points
        self.start = points[0]
        self.end = points[-1]
        self.arc_lengths = self.calc_arc_length(self.points)
        self.coeff_x, self.coeff_y = self.fit(self.points, weights)

    @property
    def arc_length_range(self):
        return self.arc_lengths[0], self.arc_lengths[-1]

    def calc_arc_length(self, points):
        arc_lengths = [0.0]
        for i in range(len(points) - 1):
            arc_length = np.sqrt(np.power(points[i].x - points[i+1].x, 2) + np.power(points[i].y - points[i+1].y, 2))
            arc_lengths.append(arc_length + arc_lengths[-1])
        return arc_lengths

    def fit(self, points, weights):
        coeff_x = np.empty((4, 1))
        coeff_y = np.empty((4, 1))
        # TODO: implement parametric polynomial fitting
        ################ your code starts #################
        xs = np.zeros((len(points), 1))
        ys = np.zeros((len(points), 1))
        S = np.ones((len(points), 4))
        s = self.calc_arc_length(points)
        w = weights
        for i in range(len(points) - 1):
            xs[i] = points[i].x
            ys[i] = points[i].y
            S[i][1] = s[i]
            S[i][2] = pow(s[i], 2)
            S[i][3] = pow(s[i], 3)

        St = S.transpose()
        Stw = St.dot(w)
        StwS = Stw.dot(S)
        StwS_inv = np.linalg.inv(StwS)
        StwS_invSt = StwS_inv.dot(St)
        StwS_invStw = StwS_invSt.dot(w)
        coeff_x = StwS_invStw.dot(xs)
        coeff_y = StwS_invStw.dot(ys)
        ################# your code ends ##################
        return coeff_x, coeff_y

    def get_point(self, s):
        s_vector = np.array([1.0, s, s * s, s * s * s])
        x = np.dot(s_vector, self.coeff_x)
        y = np.dot(s_vector, self.coeff_y)
        return Point(x, y, 0.0)


class PolynomialCurve(object):
    """
    Three order polynomial curve fitting for lane marking.
    u = c0 + c1 * v + c2 * v^2 + c3 * v^3

    """
    def __init__(self, start_v, end_v, points):
        self.start_v = start_v
        self.end_v = end_v
        self.c0, self.c1, self.c2, self.c3 = self.fit(points)

    def fit(self, points):
        num_of_points = len(points)
        v_matrix = np.empty((num_of_points, 4))
        u_vector = np.empty((num_of_points, 1))
        for i, p in enumerate(points):
            v_matrix[i, 0] = 1.0
            v_matrix[i, 1] = p.v
            v_matrix[i, 2] = p.v * p.v
            v_matrix[i, 3] = p.v * p.v * p.v
            u_vector[i, 0] = p.u
        v_matrix_transpose = np.transpose(v_matrix)
        v_matrix_square = np.dot(v_matrix_transpose, v_matrix)
        v_t_u = np.dot(v_matrix_transpose, u_vector)
        coeffs = np.dot(np.linalg.inv(v_matrix_square), v_t_u)
        return coeffs[0], coeffs[1], coeffs[2], coeffs[3]

    def get_points(self):
        step = 1 if self.end_v > self.start_v else -1
        points = []
        for v in range(self.start_v, self.end_v + step, step):
            u = self.c0 + self.c1 * v + self.c2 * v * v + self.c3 * v * v * v
            u = int(u)
            points.append(Pixel(u, v))
        return points
