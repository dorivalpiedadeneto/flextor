import unittest

from flextor import Point
from flextor import Line
from flextor import Vertex
from flextor import Segment

class TestPoint(unittest.TestCase):

    def testConstructor(self):
        pts = ((0.0, 0.0), (1.0, 1.0), (1.0, -1.0), (-1.0, 1.0))
        for pt in pts:
            z, y = pt
            p = Point(z, y)
            self.assertIsInstance(p, Point)
            self.assertAlmostEqual(z, p.z)
            self.assertAlmostEqual(y, p.y)
    
    def testString(self):
        _str = "point at z = {:.3f}, y = {:.3f}"
        pts = ((0.0, 0.0), (1.001, -2.023), (-0.12344, 32.3245))
        for pt in pts:
            z, y = pt
            p = Point(z, y)
            self.assertEqual(_str.format(z, y), str(p))

    def testEqual(self):
        pts = ((2.0, 4.0), (3.5, -2.0), (2.987, -1.234), (9.879, -41.23))
        ds = ((0.0, Point.tolerance/2.0), (Point.tolerance/12.0, 0.0),
             (Point.tolerance/1.1, Point.tolerance/1.1))
        for d in ds:
            dz, dy = d
            for pt in pts:
                z, y = pt
                p = Point(z+dz, y+dy)
                self.assertTrue(p == Point(z, y))
        ds = ((0.0, Point.tolerance*1.1), (Point.tolerance*1.1, 0.0),
             (Point.tolerance*1.01, Point.tolerance*1.01))
        for d in ds:
            dz, dy = d
            for pt in pts:
                z, y = pt
                p = Point(z+dz, y+dy)
                self.assertFalse(p == Point(z, y))

    def testMove(self):
        pts = ((0.0, 0.0), (-2.3, 3.4), (1.27, 2.34), (-0.79, -12.98))
        ds = ((0.0, 2.0), (2.34, 0.0), (-1.238, -0.897), (1923.093, -123.9321))
        for d in ds:
            dz, dy = d
            for pt in pts:
                z, y = pt
                p = Point(z, y)
                p_ = Point(z+dz, y+dy)
                p.move(dz, dy)
                self.assertTrue(p, p_)
    
    def testRotate(self):
        cs = ((0.0, 0.0), (1.0, 2.0), (-2.1, -1.2), (-1.0, 3.4))
        for cz, cy in cs:
            pc = Point(cz, cy)
            p = Point(cz + 1.0, cy)
            p.rotate(pc, 90.0, True)
            self.assertAlmostEqual(p.z - pc.z, 0.0)
            self.assertAlmostEqual(p.y - pc.y, -1.0)
            p.rotate(pc, 90.0, True)
            self.assertAlmostEqual(p.z - pc.z, -1.0)
            self.assertAlmostEqual(p.y - pc.y, 0.0)
            p.rotate(pc, 90.0, True)
            self.assertAlmostEqual(p.z - pc.z, 0.0)
            self.assertAlmostEqual(p.y - pc.y, 1.0)
            p.rotate(pc, 90.0, True)
            self.assertAlmostEqual(p.z - pc.z, 1.0)
            self.assertAlmostEqual(p.y - pc.y, 0.0)
            sqrt2_2 = 0.7071067811865476
            sqrt3_2 = 0.8660254037844386

            p.rotate(pc, -45.0, True)
            self.assertAlmostEqual(p.z - pc.z, sqrt2_2)
            self.assertAlmostEqual(p.y - pc.y, sqrt2_2)
            p.rotate(pc, -45.0, True)
            self.assertAlmostEqual(p.z - pc.z, 0.0)
            self.assertAlmostEqual(p.y - pc.y, 1.0)
            p.rotate(pc, -45.0, True)
            self.assertAlmostEqual(p.z - pc.z,-sqrt2_2)
            self.assertAlmostEqual(p.y - pc.y, sqrt2_2)
            p.rotate(pc, -45.0, True)
            self.assertAlmostEqual(p.z - pc.z,-1.0)
            self.assertAlmostEqual(p.y - pc.y, 0.0)
            p.rotate(pc, -30.0, True)
            self.assertAlmostEqual(p.z - pc.z, -sqrt3_2)
            self.assertAlmostEqual(p.y - pc.y, -0.5)
            p.rotate(pc, -30.0, True)
            self.assertAlmostEqual(p.z - pc.z, -0.5)
            self.assertAlmostEqual(p.y - pc.y, -sqrt3_2)            
            p.rotate(pc, -30.0, True)
            self.assertAlmostEqual(p.z - pc.z, 0.0)
            self.assertAlmostEqual(p.y - pc.y, -1.0)

    def testDistance(self):
        ps = ((0.0, 0.0), (1.0, 2.0), (-2.1, -1.2), (-1.0, 3.4))
        ds = ((2.0, 0.0), (-2.3, 4.5), (10.32, -9.876))
        for pz, py in ps:
            p = Point(pz, py)
            for dz, dy in ds:
                p_ = Point(pz + dz, py + dy)
                dist = (dz**2+dy**2)**0.5
                ang = 0.0
                for _ in range(36):
                    self.assertAlmostEqual(p.distance(p_), dist)
                    ang += 10.0
                    p_.rotate(p, ang)

class testLine(unittest.TestCase):

    def testConstructor(self):
        lns = (((0.0, 0.0), (1.0, 1.0)),
               ((1.0, 1.0), (1.0, -1.0)),
               ((1.0, -1.0), (-1.0, 1.0)),
               ((-1.0, 1.0), (0.0, 0.0)))
        for ln_ in lns:
            ((zi, yi),(zj, yj)) = ln_
            pi = Point(zi, yi)
            pj = Point(zj, yj)
            l = Line(pi, pj)
            self.assertIsInstance(l.pi, Point)
            self.assertIsInstance(l.pj, Point)
            self.assertIsInstance(l, Line)
            self.assertIs(l.pi, pi)
            self.assertIs(l.pj, pj)
        
    def testString(self):
        _str = "line from {} to {}"
        lns = (((0.0, 0.0), (1.0321, 1.0213)),
               ((1.04143, 1.0214), (1.03213, -1.012)),
               ((1.0321, -1.0123), (-1.0, 1.0)),
               ((-1.0321, 1.0321), (0.0, 0.0)))
        for ln_ in lns:
            ((zi, yi),(zj, yj)) = ln_
            pi = Point(zi, yi)
            pj = Point(zj, yj)
            l = Line(pi, pj)
            self.assertEqual(_str.format(pi, pj), str(l))

    def testLength(self):
        sqrt2 = 2.0 ** 0.5
        lns = (((0.0,0.0),(0.0, 1.0)),
               ((0.0,0.0),(2.0, 0.0)),
               ((0.0,0.0),(0.0,-3.0)),
               ((0.0,0.0),(-2.0, 0.0)),
               ((0.0,0.0),(2.0, 2.0)),
               ((0.0,0.0),(-3.0, -3.0)))
        ls = (1.0, 2.0, 3.0, 2.0, 2.0*sqrt2, 3.0*sqrt2)
        rp = Point(0.0, 0.0)
        for ln_,l_ in zip(lns,ls):
            ((zi, yi),(zj, yj)) = ln_
            pi = Point(zi, yi)
            pj = Point(zj, yj)
            l = Line(pi, pj)
            self.assertAlmostEqual(l.length(), l_)
            ds = ((1.0, 1.0),(-2.0, 0.0),(0.0, -2.0),(0.0, 3.0))
            for dz, dy in ds:
                l.move(dz, dy)
                rp.move(dz, dy)
                self.assertAlmostEqual(l.length(), l_)
                ang = 0.0
                for _ in range(36):
                    self.assertAlmostEqual(l.length(), l_)
                    l.rotate(rp, 10.0)
                    ang += 10.0

    def testUnitTangent(self):
        from math import radians, sin, cos
        l = Line(Point(0.0, 0.0), Point(3.0, 0.0))
        rp = Point(0.0, 0.0)
        ds = ((1.0, 1.0),(-2.0, 0.0), (0.0, -2.0), (2.0, 0.0))
        for dz, dy in ds:
            l.move(dz, dy)
            rp.move(dz, dy)
            ang = 0.0
            for _ in range(36):
                tz, ty = l.unit_tangent()
                #print(ang, tz, cos(radians(-ang)), ty, sin(radians(-ang)))
                self.assertAlmostEqual(tz, cos(radians(-ang)))
                self.assertAlmostEqual(ty, sin(radians(-ang)))
                ang += 10.0
                l.rotate(rp, 10.0)

    def testUnitNormal(self):
        from math import radians, sin, cos
        l = Line(Point(0.0, 0.0), Point(3.0, 0.0))
        rp = Point(0.0, 0.0)
        ds = ((1.0, 1.0),(-2.0, 0.0), (0.0, -2.0), (2.0, 0.0))
        for dz, dy in ds:
            l.move(dz, dy)
            rp.move(dz, dy)
            ang = 0.0
            for _ in range(36):
                nz, ny = l.unit_normal()
                #print(ang, tz, cos(radians(-ang)), ty, sin(radians(-ang)))
                self.assertAlmostEqual(nz, cos(radians(-ang+90.0)))
                self.assertAlmostEqual(ny, sin(radians(-ang+90.0)))
                l.rotate(rp, 10.0)
                ang += 10.0

    def testProjection(self):
        l = Line(Point(-1.0, 0.0), Point(1.0, 0.0))
        p = Point(-2.0, 0.0)
        qsi = -2.0
        while p.z < 2.0:
            self.assertAlmostEqual(l.projection(p), qsi)
            p.move(0.1, 0.0)
            qsi += 0.1
        p.move(0.0, 1.0)
        while p.z < 2.0:
            self.assertAlmostEqual(l.projection(p), qsi)
            p.move(-0.1, 0.0)
            qsi -= 0.1
        l.move(-10.0, 2.0)
        p.move(-10.0, 0.0)
        while p.z < 4.0:
            self.assertAlmostEqual(l.projection(p), qsi)
            p.move(0.1, 0.0)
            qsi += 0.1

    def testCoordinate(self):
        l = Line(Point(-1.0,0.0), Point(1.0,0.0))
        z_ = -1.0; y_ = 0.0; qsi_ = -1.0
        for _ in range(21):
            z, y = l.coordinate(qsi_)
            self.assertAlmostEqual(z, z_)
            self.assertAlmostEqual(y, y_)
            qsi_ += 0.1
            z_ += 0.1
            y_ += 0.0
        l.pi.move(10.0, 1.0)
        l.pj.move(10.0,-1.0)
        z_ = 9.0; y_ = 1.0; qsi_ = -1.0
        for _ in range(21):
            z, y = l.coordinate(qsi_)
            self.assertAlmostEqual(z, z_)
            self.assertAlmostEqual(y, y_)
            qsi_ += 0.1
            z_ += 0.1
            y_ -= 0.1
        l.pi.move(-20.0,-2.0)
        l.pj.move(-20.0,-2.0)
        z_ =-11.0; y_ = -1.0; qsi_ = -1.0
        for _ in range(21):
            z, y = l.coordinate(qsi_)
            self.assertAlmostEqual(z, z_)
            self.assertAlmostEqual(y, y_)
            qsi_ += 0.1
            z_ += 0.1
            y_ -= 0.1
    
    def testDistance(self):
        l = Line(Point(1.0, 1.0), Point(2.0, 1.0))
        self.assertAlmostEqual(l.distance(Point(0.0, 0.0)), 2.0**0.5)
        self.assertAlmostEqual(l.distance(Point(0.5, 0.5)), 0.5*2.0**0.5)
        self.assertAlmostEqual(l.distance(Point(0.5, 1.0)), 0.5)
        self.assertAlmostEqual(l.distance(Point(0.5, 1.5)), 0.5*2.0**0.5)
        self.assertAlmostEqual(l.distance(Point(1.0, 0.5)), 0.5)
        self.assertAlmostEqual(l.distance(Point(1.0, 1.5)), 0.5)
        self.assertAlmostEqual(l.distance(Point(1.0, 1.0)), 0.0)
        self.assertAlmostEqual(l.distance(Point(2.0, 0.5)), 0.5)
        self.assertAlmostEqual(l.distance(Point(2.0, 1.0)), 0.0)
        l = Line(Point(-1.0, 0.0), Point(0.0, -1.0))
        self.assertAlmostEqual(l.distance(Point(0.0, 0.0)), 0.5*2.0**0.5)
        self.assertAlmostEqual(l.distance(Point(-1.0, -1.0)), 0.5*2.0**0.5)
        self.assertAlmostEqual(l.distance(Point(-0.5, -0.5)), 0.0)  

    def testBoundingBox(self):
        l = Line(Point(1.0, 0.0), Point(0.0, 0.0))
        ((zmax_, ymax_),(zmin_, ymin_)) = l.bounding_box()
        self.assertAlmostEqual(zmax_, 1.0)
        self.assertAlmostEqual(zmin_, 0.0)
        self.assertAlmostEqual(ymax_, 0.0)
        self.assertAlmostEqual(ymin_, 0.0)
        l = Line(Point(-1.0, 0.0), Point(-3.5, -0.1))
        ((zmax_, ymax_),(zmin_, ymin_)) = l.bounding_box()
        self.assertAlmostEqual(zmax_, -1.0)
        self.assertAlmostEqual(zmin_, -3.5)
        self.assertAlmostEqual(ymax_, 0.0)
        self.assertAlmostEqual(ymin_, -0.1)
        l = Line(Point(1.0, 0.23), Point(-7.5, -1.1))
        ((zmax_, ymax_),(zmin_, ymin_)) = l.bounding_box()
        self.assertAlmostEqual(zmax_, 1.0)
        self.assertAlmostEqual(zmin_, -7.5)
        self.assertAlmostEqual(ymax_, 0.23)
        self.assertAlmostEqual(ymin_, -1.1)
    
    def testIsInsideBbox(self):
        l = Line(Point(3.0,2.0),Point(-1.0, -0.5))
        self.assertTrue(l.is_inside_bbox(Point(0.0, 0.0)))
        self.assertTrue(l.is_inside_bbox(Point(1.0, 1.0)))
        self.assertTrue(l.is_inside_bbox(Point(3.0, 2.0)))
        self.assertTrue(l.is_inside_bbox(Point(-1.0,-0.5)))
        self.assertTrue(l.is_inside_bbox(Point(3.0, -0.5)))
        self.assertFalse(l.is_inside_bbox(Point(3.01, 0.0)))
        self.assertFalse(l.is_inside_bbox(Point(3.00, -0.51)))
        self.assertFalse(l.is_inside_bbox(Point(3.0, 3.0)))


class TestVertex(unittest.TestCase):

    def testConstructor(self):
        vts = (('a',(0.0, 0.0)), ('12',(1.0, 1.0)),
               ('ii',(1.0, -1.0)), ('ca',(-1.0, 1.0)))
        for vt in vts:
            nm, (z, y) = vt
            v = Vertex(nm, z, y)
            self.assertIsInstance(v, Vertex)
            self.assertIsInstance(v, Point)
            self.assertEqual(nm, v.name)
            self.assertAlmostEqual(z, v.z)
            self.assertAlmostEqual(y, v.y)
    
    def testString(self):
        _str = "vertex {} at z = {:.3f}, y = {:.3f}"
        vts = (('a',(0.0, 0.0)), ('12',(1.0, 1.0)),
               ('ii',(1.0, -1.0)), ('ca',(-1.0, 1.0)))
        for vt in vts:
            nm, (z, y) = vt
            v = Vertex(nm, z, y)
            self.assertEqual(_str.format(nm, z, y), str(v))

class TestSegment(unittest.TestCase):

    def testConstructor(self):
        lns = (('a',(0.0, 0.0), (1.0, 1.0)),
               ('3',(1.0, 1.0), (1.0, -1.0)),
               ('II',(1.0, -1.0), (-1.0, 1.0)),
               ('iv',(-1.0, 1.0), (0.0, 0.0)))
        for l in lns:
            nm, (zi, yi), (zj, yj) = l
            s = Segment(nm, Point(zi, yi), Point(zj, yj))
            self.assertIsInstance(s, Segment)
            self.assertIsInstance(s, Line)
            self.assertIsInstance(s.pi, Point)
            self.assertIsInstance(s.pj, Point)

    def testString(self):
        _str = "straight segment {} from {} to {}"
        lns = (('a',(0.0, 0.0), (1.0, 1.0)),
               ('3',(1.0, 1.0), (1.0, -1.0)),
               ('II',(1.0, -1.0), (-1.0, 1.0)),
               ('iv',(-1.0, 1.0), (0.0, 0.0)))
        for l in lns:
            nm, (zi, yi), (zj, yj) = l
            s = Segment(nm, Point(zi, yi), Point(zj, yj))
            self.assertEqual(_str.format(nm, Point(zi, yi), Point(zj, yj)), str(s))

    def testThickness(self):
        s = Segment('A', Point(0.0, 0.0), Point(1.0, 1.0))
        self.assertIsNone(s.thickness)
        s.thickness = 1.0
        self.assertIsInstance(s.thickness, float)
        self.assertAlmostEqual(s.thickness, 1.0)
        s.thickness = None
        self.assertIsNone(s.thickness)

    def testVertices(self):
        s = Segment('A', Point(0.0, 0.0), Point(1.0, 1.0))
        self.assertIsNone(s.first_vertex)
        self.assertIsNone(s.last_vertex)
        A = Vertex('A', 0.0, 0.0)
        B = Vertex('B', 1.0, 1.0)
        s.first_vertex = A
        s.last_vertex = B
        self.assertIsInstance(s.first_vertex, Vertex)
        self.assertIsInstance(s.last_vertex, Vertex)
        self.assertIs(s.first_vertex, A)
        self.assertIs(s.last_vertex, B)

    def testProperties(self):
        s = Segment('A', Point(0.0, 6.0), Point(0.0, -6.0))
        s.thickness = 2.0
        self.assertFalse(s.has_property('area'))
        self.assertFalse(s.has_property('Izl'))
        self.assertFalse(s.has_property('Iyl'))
        self.assertFalse(s.has_property('Iz'))
        self.assertFalse(s.has_property('Iy'))
        self.assertFalse(s.has_property('Izy'))
        self.assertFalse(s.has_property('CG'))
        self.assertIsNone(s.get_property('area'))
        self.assertIsNone(s.get_property('Izl'))
        self.assertIsNone(s.get_property('Iyl'))
        self.assertIsNone(s.get_property('Iz'))
        self.assertIsNone(s.get_property('Iy'))
        self.assertIsNone(s.get_property('Izy'))
        self.assertIsNone(s.get_property('CG'))
        s.compute_properties()
        self.assertTrue(s.has_property('area'))
        self.assertTrue(s.has_property('Izl'))
        self.assertTrue(s.has_property('Iyl'))
        self.assertTrue(s.has_property('Iz'))
        self.assertTrue(s.has_property('Iy'))
        self.assertTrue(s.has_property('Izy'))
        self.assertTrue(s.has_property('CG'))
        self.assertIsNotNone(s.get_property('area'))
        self.assertIsNotNone(s.get_property('Izl'))
        self.assertIsNotNone(s.get_property('Iyl'))
        self.assertIsNotNone(s.get_property('Iz'))
        self.assertIsNotNone(s.get_property('Iy'))
        self.assertIsNotNone(s.get_property('Izy'))
        self.assertIsNotNone(s.get_property('CG'))
        # Testing basic computations
        self.assertAlmostEqual(s.get_property('area'), 24.0)
        self.assertAlmostEqual(s.get_property('Izl'), 288.0)
        self.assertAlmostEqual(s.get_property('Iyl'), 8.0)
        self.assertAlmostEqual(s.get_property('Iz'), 288.0)
        self.assertAlmostEqual(s.get_property('Iy'), 8.0)
        self.assertAlmostEqual(s.get_property('Izy'), 0.0)
        self.assertEqual(s.get_property('CG'), Point(0.0, 0.0))

    def testComputeProperties(self):
        s = Segment('A', Point(0.0,6.0), Point(0.0, -6.0))
        s.thickness = 2.0
        ctr = Point(0.0, 0.0)
        # Rotating 45 degrees 8 times (theta = 0, 45, 90, 135, 180, 225, ...)
        Izs = (288.0, 148.0, 8.0, 148.0, 288.0, 148.0, 8.0, 148.0, 288.0)
        Iys = (8.0, 148.0, 288.0, 148.0, 8.0, 148.0, 288.0, 148.0, 8.0)
        Izys = (0.0, 140.0, 0.0, -140.0, 0.0, 140.0, 0.0, -140.0, 0.0)
        s.rotate(ctr,-45.0, True)
        for Iz, Iy, Izy in zip(Izs, Iys, Izys):
            s.rotate(ctr, 45.0, True)
            s.compute_properties()
            self.assertAlmostEqual(s.get_property('area'), 24.0)
            self.assertAlmostEqual(s.get_property('Izl'), 288.0)
            self.assertAlmostEqual(s.get_property('Iyl'), 8.0)
            self.assertAlmostEqual(s.get_property('Iz'), Iz)
            self.assertAlmostEqual(s.get_property('Iy'), Iy)
            self.assertAlmostEqual(s.get_property('Izy'), Izy)
            self.assertEqual(s.get_property('CG'),ctr)
        # Now testing also moving the position
        dzs = (2.0,-4.0, 0.0, 4.0)
        dys = (2.0, 0.0,-4.0, 0.0)
        for dz, dy in zip(dzs, dys):
            s.move(dz,dy)
            ctr.move(dz,dy)
            Izs = (288.0, 148.0, 8.0, 148.0, 288.0, 148.0, 8.0, 148.0, 288.0)
            Iys = (8.0, 148.0, 288.0, 148.0, 8.0, 148.0, 288.0, 148.0, 8.0)
            Izys = (0.0, 140.0, 0.0, -140.0, 0.0, 140.0, 0.0, -140.0, 0.0)
            s.rotate(ctr,-45.0, True)
            for Iz, Iy, Izy in zip(Izs, Iys, Izys):
                s.rotate(ctr, 45.0, True)
                s.compute_properties()
                self.assertAlmostEqual(s.get_property('area'), 24.0)
                self.assertAlmostEqual(s.get_property('Izl'), 288.0)
                self.assertAlmostEqual(s.get_property('Iyl'), 8.0)
                self.assertAlmostEqual(s.get_property('Iz'), Iz)
                self.assertAlmostEqual(s.get_property('Iy'), Iy)
                self.assertAlmostEqual(s.get_property('Izy'), Izy)
                self.assertEqual(s.get_property('CG'), ctr)


if __name__ == "__main__":
    unittest.main()
