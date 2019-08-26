import numpy as np

from src.utils import Pixel


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
        c0 = 0.0
        c1 = 0.0
        c2 = 0.0
        c3 = 0.0
        ######################################################
        # TODO: implement your polynomial curve fitting here #
        ################ your code starts ####################


        ################# your code ends #####################

        return c0, c1, c2, c3

    def get_points(self):
        step = 1 if self.end_v > self.start_v else -1
        points = []
        for v in range(self.start_v, self.end_v + step, step):
            u = self.c0 + self.c1 * v + self.c2 * v * v + self.c3 * v * v * v
            u = int(u)
            points.append(Pixel(u, v))
        return points
