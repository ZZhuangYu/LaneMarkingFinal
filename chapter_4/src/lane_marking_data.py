import random

import numpy as np

from src.lane_marking import LaneMarking
from src.geometry import Point, Transform, Rotation, Translation


def load_lane_markings():
    xs = [float(i) for i in range(100)]

    def left_y(x):
        return 1.5 + 0.025 * x + 0.001 * x * x + 0.1 * np.sin(x * 0.1) + random.random() * 0.5

    def right_y(x):
        return -1.5 + 0.01 * x + 0.001 * x * x + 0.1 * np.sin(x * 0.2) + random.random() * 0.5

    left_ys = [left_y(x) for x in xs]
    right_ys = [right_y(x) for x in xs]

    raw_points = [Point(x, y, 0.0) for x, y in zip(xs + xs, left_ys + right_ys)]

    first_left = [Point(x, y, 0.0) for x, y in zip(xs[:60], left_ys[:60])]
    first_right = [Point(x, y, 0.0) for x, y in zip(xs[:60], right_ys[:60])]
    first_left_lane_marking = LaneMarking(first_left)
    first_right_lane_marking = LaneMarking(first_right)

    trans_x = 50.0
    trans_y = 0.5 * (left_y(trans_x) + right_y(trans_x))
    yaw = 1.46
    pose = Transform(Rotation.from_rpy(0.0, 0.0, yaw), Translation(trans_x, trans_y, 0.0))
    pose_inv = np.linalg.inv(pose.homogeneous)
    second_left = []
    second_right = []
    for x, left_y, right_y in zip(xs[40:], left_ys[40:], right_ys[40:]):
        left_p = np.array([x, left_y, 0.0, 1.0])
        left_p = left_p.reshape((4, 1))
        right_p = np.array([x, right_y, 0.0, 1.0])
        right_p = right_p.reshape((4, 1))
        left_p_inv = np.dot(pose_inv, left_p)
        right_p_inv = np.dot(pose_inv, right_p)
        second_left.append(Point(left_p_inv[0, 0], left_p_inv[1, 0] + 0.5 * random.random(), 0.0))
        second_right.append(Point(right_p_inv[0, 0], right_p_inv[1, 0] + 0.5 * random.random(), 0.0))
    second_left_lane_marking = LaneMarking(second_left, pose=pose)
    second_right_lane_marking = LaneMarking(second_right, pose=pose)

    return (first_left_lane_marking, first_right_lane_marking), (second_left_lane_marking, second_right_lane_marking), raw_points
