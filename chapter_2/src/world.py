from src.geometry import Point
from src.camera import create_test_camera_model

import cv2 as cv
import numpy as np


class ColoredLine(object):
    def __init__(self, start, end, color=(0, 255, 0), width=5):
        self.start = start
        self.end = end
        self.color = color
        self.width = width

    def project_to_image(self, image, camera):
        start_pixel = camera.project(self.start)
        end_pixel = camera.project(self.end)
        start_uv = start_pixel.uv
        end_uv = end_pixel.uv
        print(start_uv)
        print(end_uv)
        image = cv.line(image, start_uv, end_uv, self.color, thickness=self.width)
        return image


class ColoredPolygon(object):
    def __init__(self, points, color=(255, 255, 0)):
        self.points = points
        self.color = color
        self.num_of_points = len(self.points)
        if self.num_of_points < 3:
            raise Exception("Polygon needs at lease 3 points")

    def project_to_image(self, image, camera):
        pixels = np.empty((self.num_of_points, 2))
        for i, point in enumerate(self.points):
            pixel = camera.project(point)
            pixels[i, 0] = pixel.uv[0]
            pixels[i, 1] = pixel.uv[1]
        pixels = pixels.astype(int)
        print(pixels)
        image = cv.fillConvexPoly(image, pixels, self.color)
        return image


class Circle(object):
    def __init__(self, center, radius, color=(0, 0, 255)):
        self.center = center
        self.radius = radius
        self.color = color

    def project_to_image(self, image, camera):
        # TODO: project circle to image
        ### your code starts  ###
        w_center = self.center.to_vector()
        x0 = int(w_center[0][0])
        y0 = int(w_center[1][0])
        z0 = int(w_center[2][0])
        p0 = np.array([x0, y0, z0])
        p1 = np.array([x0 + self.radius, y0, z0])
        p2=  np.array([x0, y0 + self.radius, z0])  # Need to trans P1, P2(type:ndarray) to Point
        print('Center Point of the circle:' + str(p0))
        print('One of the perpendicular radius cross circle-rim point is:' + str(p1))
        print('Another perpendicular radius cross circle-rim point is:' + str(p2))
        P0 = Point(x0, y0, z0)
        P1 = Point(x0 + self.radius, y0, z0)
        P2 = Point(x0, y0 + self.radius, z0)
        r0 = camera.project(P0)
        r1 = camera.project(P1)
        r2 = camera.project(P2)
        R0 = r0.to_vector()
        R00 = np.array([R0[0][0][0], R0[1][0][0], R0[2][0]])
        R1 = r1.to_vector()
        R11 = np.array([R1[0][0][0], R1[1][0][0], R1[2][0]])
        #R1 = np.array([r1.to_vector()[0][0][0], r1.to_vector()[1][0][0]], r1.to_vector()[2][0]])
        R2 = r2.to_vector()
        #R2 = np.array([r2.to_vector()[0][0][0], r2.to_vector()[1][0][0]], r2.to_vector()[2][0]])
        R22 = np.array([R2[0][0][0], R2[1][0][0], R2[2][0]])
        print('Center Point of the projected ellipse:' + str(R00))
        print('One of the perpendicular radius cross ellipse-rim point after projection is:' + str(R11))
        print('Another perpendicular radius cross ellipse-rim point after projection is:' + str(R22))
        e_center = (int(R00[0]), int(R00[1]))
        d1 = int(np.linalg.norm(R11-R00))
        d2 = int(np.linalg.norm(R22-R00))
        image = cv.ellipse(image, e_center, (d1, d2), 0, 0, 360, self.color, 5)

        ### your code ends    ###
        return image


if __name__ == "__main__":
    camera = create_test_camera_model()

    line = ColoredLine(Point(-1.0, -0.5, 1.5), Point(-0.5, -1.0, 3.5))

    points = [Point(-0.5, 0.0, 2.0), Point(-0.5, 0.0, 4.0), Point(0.5, 0.5, 4.0), Point(0.5, 0.5, 2.0)]
    polygon = ColoredPolygon(points)

    circle = Circle(Point(0.0, 0.0, 0.0), 0.5)

    image = np.ones((768, 1024, 3)) * 255
    image = polygon.project_to_image(image, camera)
    image = line.project_to_image(image, camera)
    image = circle.project_to_image(image, camera)

    cv.imshow("image", image)
    cv.waitKey(0)
