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
            z, y = self.pi.coord()
        elif qsi > 1.0:
            z, y = self.pj.coord()
        else:
            z, y = self.coordinate(qsi)
        zp, yp = other_point.coord()
        d = sqrt((zp - z) ** 2 + (yp - y) ** 2)
        return d

    def bounding_box(self):
        '''
        The bounding box (tuple of two tuples) return the bottom left and the
        top right coordinate of the line's bounding box (z horizontal, from 
        right to left, and y vertical, from top to bottom).
        '''
        zi, yi = self.pi.coord()
        zj, yj = self.pj.coord()
        return((max(zi, zj),max(yi, yj)),(min(zi, zj), min(yi, yj)))

    def is_inside_bbox(self, point, tolerance = Point.tolerance):
        z_, y_ = point.coord()
        zi, yi = self.pi.coord()
        zj, yj = self.pj.coord()
        zmax = max(zi, zj) + tolerance; ymax = max(yi, yj) + tolerance
        zmin = min(zi, zj) - tolerance; ymin = min(yi, yj) - tolerance
        return (zmin < z_ < zmax) and (ymin < y_ < ymax)

class Vertex(Point):

    def __init__(self, name, z, y):
        self.name = name        # planning to use string to make reference
        super().__init__(z, y)

    def __repr__(self):
        return "vertex {} at z = {:.3f}, y = {:.3f}".format(self.name, self.z, self.y)
    
    def __str__(self):
        return "vertex {} at z = {:.3f}, y = {:.3f}".format(self.name, self.z, self.y)

class Segment(Line):

    def __init__(self, name, pi, pj):
        self.name = name        # planning to use string to make reference
        super().__init__(pi, pj)
        self.thickness = None
        self.first_vertex = None
        self.last_vertex = None
        self.properties = {}

    def __repr__(self):
        return 'straight segment {} from {} to {}'.format(self.name, self.pi, self.pj)

    def __str__(self):
        return 'straight segment {} from {} to {}'.format(self.name, self.pi, self.pj)

    def has_property(self, name):
        return name in self.properties

    def get_property(self, name):
        if name in self.properties:
            return self.properties[name]
        else:
            return None
    

class Cross_section(object):

    def __init__(self):
        self.segments = {}      # dict of segment instances
        self.vertices = {}      # dict of vertices instances
        self.CG = None          # None or a point instance
        self.D = None           # None or a point instance 
        self.paa = None         # principal axes angle (None or float)
        self.properties = {}    # dict for storing the results


if __name__ == "__main__":
    pass