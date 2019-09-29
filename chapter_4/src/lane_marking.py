import numpy as np

from src.polynomial_curve import ParametricPolyCurve
from src.geometry import Rotation, Translation, Transform


class LaneMarking(object):
    def __init__(self, points, weights=None, pose=None):
        if weights is None:
            weights = np.eye(len(points))
        if pose is None:
            pose = Transform(Rotation.from_rpy(0.0, 0.0, 0.0), Translation(0.0, 0.0, 0.0))
        self.pose = pose
        self.polynomial = ParametricPolyCurve(points, weights)
        self.start_s, self.end_s = self.polynomial.arc_length_range

    def get_points(self, start_s, end_s, step=0.1):
        if start_s < self.start_s or end_s > self.end_s:
            raise Exception("s out of range")
        points = []
        arc_lengths = []
        s = start_s
        while s <= end_s:
            arc_lengths.append(s)
            point = self.polynomial.get_point(s)
            # transform point to world coordinates
            point_in_world = self.pose.apply_to_point(point)
            points.append(point_in_world)
            s += step
        return points, arc_lengths

    def get_all_points(self, step=0.1):
        return self.get_points(self.start_s, self.end_s, step)

    def get_weights(self, arc_lengths):
        weights = np.eye(len(arc_lengths))
        # TODO: implement information matrix
        for i in range(len(arc_lengths) - 1):
            weights[i, i] = 1 - arc_lengths[i]/arc_lengths[-1]
        # trying the -1 power of arclength:
        #weights[0, 0] = 1.0
        #weights[1, 1] = 1.0
        #for i in range(2, len(arc_lengths)):
           # weights[i, i] = weights[i-1, i-1] * arc_lengths[i-1]/arc_lengths[i]
        return weights
