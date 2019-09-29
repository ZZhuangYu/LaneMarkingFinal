import numpy as np

from src.lane_marking import LaneMarking
from src.lane_marking_data import load_lane_markings


class LaneMarkingFuser(object):
    def __init__(self):
        self.left_lane_markings = []
        self.right_lane_markings = []

    def process_observation(self, left_lane_marking, right_lane_marking):
        # initialize fuser
        if not self.left_lane_markings or not self.right_lane_markings:
            self.left_lane_markings.append(left_lane_marking)
            self.right_lane_markings.append(right_lane_marking)
            return
        # sample observation lane markings
        left_points, left_s_values = left_lane_marking.get_all_points()
        right_points, right_s_values = right_lane_marking.get_all_points()
        left_weights = left_lane_marking.get_weights(left_s_values)
        right_weights = right_lane_marking.get_weights(right_s_values)
        # get first points of observation lane markings
        left_start_point = left_points[0]
        right_start_point = right_points[0]
        # find starting arc length on previous lane markings
        left_start_s = self.find_start_s(self.left_lane_markings[-1], left_start_point)
        right_start_s = self.find_start_s(self.right_lane_markings[-1], right_start_point)
        # calculate previous lane marking points and weights
        prev_left_points, prev_left_s_values = self.left_lane_markings[-1].get_points(left_start_s, self.left_lane_markings[-1].end_s)
        prev_right_points, prev_right_s_values = self.right_lane_markings[-1].get_points(right_start_s, self.right_lane_markings[-1].end_s)
        prev_left_weights = self.left_lane_markings[-1].get_weights(prev_left_s_values)
        prev_right_weights = self.right_lane_markings[-1].get_weights(prev_right_s_values)
        # fuse previous lane marking and observation
        fused_left_lane_marking = self.fuse_observation(prev_left_points, prev_left_weights, left_points, left_weights)
        fused_right_lane_marking = self.fuse_observation(prev_right_points, prev_right_weights, right_points, right_weights)
        # add fused lane markings to list
        self.left_lane_markings.append(fused_left_lane_marking)
        self.right_lane_markings.append(fused_right_lane_marking)

    def find_start_s(self, lane_marking, point):
        start_s = 0.0
        # TODO: find arc length on lane marking with respect to given point
        ################### your code starts ####################
        prev_points, prev_ss = lane_marking.get_all_points()  #Using E-distance as judgement of fusing
        distance = np.zeros((len(prev_points), ))
        for i in range(len(prev_points) - 1):
            distance[i] = np.sqrt(pow((prev_points[i].y - point.y), 2) + pow((prev_points[i].x - point.x), 2))

        t = int(np.argwhere(distance == np.max(distance)))
        start_s = prev_ss[t]
        ################### your code ends ######################
        return start_s

    def fuse_observation(self, prev_points, prev_weights, observation_points, observation_weights):
        num_of_prev = len(prev_points)
        num_of_observation = len(observation_points)
        assert(num_of_prev < num_of_observation)
        for i in range(num_of_prev):
            x = prev_points[i].x * prev_weights[i, i] + observation_points[i].x * observation_weights[i, i]
            y = prev_points[i].y * prev_weights[i, i] + observation_points[i].y * observation_weights[i, i]
            normalized_x = x / (prev_weights[i, i] + observation_weights[i, i])
            normalized_y = y / (prev_weights[i, i] + observation_weights[i, i])
            observation_points[i].x = normalized_x
            observation_points[i].y = normalized_y
        return LaneMarking(observation_points, observation_weights)

    def get_vis_points(self):
        left_points, _ = self.left_lane_markings[-1].get_all_points()
        right_points, _ = self.right_lane_markings[-1].get_all_points()
        return left_points, right_points


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    fuser = LaneMarkingFuser()

    (first_left, first_right), (second_left, second_right), raw_points = load_lane_markings()
    plt.plot([p.y for p in raw_points], [p.x for p in raw_points], "y.")

    second_left_points, _ = second_left.get_all_points()
    second_right_points, _ = second_right.get_all_points()

    plt.plot([p.y for p in second_left_points], [p.x for p in second_left_points], "g-.")
    plt.plot([p.y for p in second_right_points], [p.x for p in second_right_points], "g-.")

    fuser.process_observation(first_left, first_right)
    left_points, right_points = fuser.get_vis_points()

    plt.plot([p.y for p in left_points], [p.x for p in left_points], "b--")
    plt.plot([p.y for p in right_points], [p.x for p in right_points], "b--")

    fuser.process_observation(second_left, second_right)
    left_points, right_points = fuser.get_vis_points()

    plt.plot([p.y for p in left_points], [p.x for p in left_points], "r-")
    plt.plot([p.y for p in right_points], [p.x for p in right_points], "r-")

    plt.grid(True)
    plt.show()
