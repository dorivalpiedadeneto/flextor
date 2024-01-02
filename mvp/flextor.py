# Minimal class (no setters/getters, no typing, no tests, minimal comments...)

from math import sin, cos, radians, sqrt

class Point(object):
    '''
    Point in 2D: y vertical (positive from top to bottom), z horizontal
    (positive from right to left). Rotation positive in clockwise direction.
    '''
    tolerance = 1.0e-6

    def __init__(self, z, y):
        self.z = z
        self.y = y
    
    def __repr__(self):
        return "point at z = {:.3f}, y = {:.3f}".format(self.z, self.y)
    
    def __str__(self):
        return "point at z = {:.3f}, y = {:.3f}".format(self.z, self.y)

    def __eq__(self, other_point):
        test_z = abs(self.z - other_point.z) < self.tolerance
        test_y = abs(self.y - other_point.y) < self.tolerance

        return test_z and test_y

    def coord(self):
        return self.z, self.y

    def move(self, dz = 0.0, dy = 0.0):
        self.z += dz
        self.y += dy
        
    def rotate(self, rpoint, angle, deg = True):
        if deg:
            angle = radians(angle)
        angle = -angle
        rz, ry = rpoint.coord()
        z, y = self.coord()
        dz = z - rz
        dy = y - ry
        self.z = rz + cos(angle) * dz - sin(angle) * dy
        self.y = ry + sin(angle) * dz + cos(angle) * dy

    def distance(self, other_point):
        oz, oy = other_point.coord()
        z, y = self.coord()
        d = sqrt((z - oz) ** 2 + (y -oy) ** 2)
        return d

class Line(object):
    '''
    Straight line in 2D space. Natural coordinate qsi from -1.0 to 1.0.
    '''

    def __init__(self, pi, pj):
        self.pi = pi
        self.pj = pj

    def __repr__(self):
        return 'line from {} to {}'.format(self.pi, self.pj)

    def __str__(self):
        return 'line from {} to {}'.format(self.pi, self.pj)

    def move(self, dz = 0.0, dy = 0.0):
        self.pi.move(dz, dy)
        self.pj.move(dz, dy)

    def rotate(self, rpoint, angle, deg = True):
        self.pi.rotate(rpoint, angle, deg)
        self.pj.rotate(rpoint, angle, deg)

    def length(self):
        zi, yi = self.pi.coord()
        zj, yj = self.pj.coord()
        return sqrt((zj - zi) ** 2 + (yj - yi) ** 2)

    def unit_tangent(self):
        zi, yi = self.pi.coord()
        zj, yj = self.pj.coord()
        l = sqrt((zj - zi) ** 2 + (yj - yi) ** 2)
        return ((zj - zi) / l, (yj - yi) / l)

    def unit_normal(self):
        tz, ty = self.unit_tangent()
        nz = -ty
        ny = tz
        return (nz, ny)       

    def projection(self, point):
        '''
        Returns the dimensionless coordinate qsi of the orthogonal projection of
        point over the line.
        qsi = -1.0 -> pi
        qsi = 1.0 -> pj
        qsi = 0.0 -> middle of the line
        qsi < -1.0 or qsi > 1.0 -> out of the line
        '''
        # point coordinate
        zp, yp = point.coord()
        # middle point (point at the middle of the line segment)
        zi, yi = self.pi.coord()
        zj, yj = self.pj.coord()
        zm = 0.5 * (zi + zj)
        ym = 0.5 * (yi + yj)
        # vector mp - from middle point to point
        zmp = zp - zm
        ymp = yp - ym
        # unit tanget vector components
        tz, ty = self.unit_tangent()
        # lenght from middle pointo to pj
        l = sqrt((zj - zm) ** 2 + (yj - ym) ** 2)
        qsi = (zmp * tz + ymp * ty) / l
        return qsi        

    def coordinate(self, qsi):
        '''
        Returns the coordinate of a point over the line in the given qsi.
        If qsi < -1.0 or qsi > 1.0, it is over the line's extension (i.e.,
        it is out of the segment).
        '''
        zi, yi = self.pi.coord()
        zj, yj = self.pj.coord()
        phii = 0.5 * (1.0 - qsi)
        phij = 0.5 * (1.0 + qsi)
        z = zi * phii + zj * phij
        y = yi * phii + yj * phij
        return (z, y)

    def distance(self, other_point):
        '''
        Returns the distance from other_point and the line. If the normal
        projection qsi is between -1.0 and 1.0, the computed distance is
        between other_point and its normal projection. If qsi < -1.0,
        returns the distance from other_point and pi; if qsi > 1.0,
        returns the distance from other_point and pj.
        '''
        qsi = self.projection(other_point)
        if qsi < -1.0:
            z, y = self.pi().coord()
        elif qsi > 1.0:
            z, y = self.pj().coord()
        else:
            z, y = self.coordinate(qsi)
        zp, yp = other_point.coord()
        d = sqrt((zp - z) ** 2 + (yp - y) ** 2)
        return d

    def bounding_box(self):
        zi, yi = self.pi.coord()
        zj, yj = self.pj.coord()
        return((max(zi, zj),max(yi, yj),(min(zi, zj), min(yi, yj))))


def test_point_dev():
    # Only during first development stage
    p = Point(1.0, -2.0)
    o = Point(1.0, -2.0)
    print(p)
    print(o)
    print("Is 'o' equal to point 'p'? {}".format(o == p))
    o.y += 1.0e-7; o.z -= 1.0e-7
    print("And now? {}".format(o == p))
    o.y += 0.1
    print("Is it false now? {}".format(o == p))
    o.y = p.y
    Point.tolerance = 1.0e-9
    print("Changed the tolerance, but z still greater -> {}".format(o==p))
    o.z = p.z
    print("Now they are equal: {}".format(o==p))
    # test move, rotate and distance
    o.move(-1.0, 2.0)
    print(o)
    p.rotate(o, 90.0)
    print(p)
    p.rotate(o, 90.0)
    print(p)
    p.rotate(o, 90.0)
    print(p)
    print('Is distance right? -> {}'.format(abs(p.distance(o) - sqrt(5)) < Point.tolerance))

def test_line_dev():
    # Only during first development stage
    p = Point(1.5, 1.5)
    l = Line(Point(1.0, 1.0), Point(1.0, 2.0))
    print(p)
    print(l)
    print(l.unit_tangent())
    print(l.unit_normal())
    print(l.projection(p))
    print(l.coordinate(l.projection(p)))


if __name__ == "__main__":
    test_point_dev()
    test_line_dev()
    #Just superficial tests... next step: write a file with more complete
    #tests