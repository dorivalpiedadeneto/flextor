import unittest

from flextor import Point
from flextor import Line
from flextor import Vertex
from flextor import Segment
from flextor import Cross_section

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

    def testDataMethods(self):
        p = Point(10.0, 2.0)
        self.assertEqual(p.data, {'z': 10.0, 'y':2.0})
        p.data = {'z': -11.3, 'y': -0.1}
        self.assertEqual(p.data, {'z': -11.3, 'y': -0.1})

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

    def testDataMethods(self):
        pi = Point(0.0, 0.0); pj = Point(1.0, -1.0)
        l = Line(pi, pj)
        self.assertEqual(l.data, {'pi': pi.data, 'pj': pj.data})
        pi_data = {'z': -1.0, 'y': 2.0}
        pj_data = {'z': -0.1, 'y': 3.2}
        l.pi.data = pi_data
        l.pj.data = pj_data
        self.assertEqual(l.data, {'pi':pi_data, 'pj': pj_data})

    def testIntersection(self):
        l = Line(Point(0.0, 0.0), Point(1.0, 1.0))
        ol = Line(Point(0.0, 1.0), Point(1.0, 0.0))
        res = l.intersection(ol)
        self.assertIsInstance(res, tuple)
        self.assertAlmostEqual(res[0], 0.0)
        self.assertAlmostEqual(res[1], 0.0)
        # Parallel lines
        l = Line(Point(0.0, 0.0), Point(0.0, 1.0))
        ol = Line(Point(1.0, 0.0), Point(1.0, 1.0))
        self.assertIsNone(l.intersection(ol))
        # Aligned lines
        l = Line(Point(0.0, 0.0), Point(0.0, 1.0))
        ol = Line(Point(0.0, 1.0), Point(0.0, 0.0))
        self.assertIsNone(l.intersection(ol))
        # Testing qsi sign
        l = Line(Point(0.0, 0.0), Point(1.0, 0.0))
        ol = Line(Point(0.25, -1.0), Point(0.25, 1.0))
        res = l.intersection(ol)
        self.assertIsInstance(res, tuple)
        self.assertAlmostEqual(res[0], -0.5)
        self.assertAlmostEqual(res[1], 0.0)
        ol.move(dz=0.25,dy=0.0)
        res = l.intersection(ol)
        self.assertIsInstance(res, tuple)
        self.assertAlmostEqual(res[0], 0.0)
        self.assertAlmostEqual(res[1], 0.0)
        ol.move(dz=0.25,dy=0.0)
        res = l.intersection(ol)
        self.assertIsInstance(res, tuple)
        self.assertAlmostEqual(res[0], 0.5)
        self.assertAlmostEqual(res[1], 0.0)
        ol.move(dz=0.25,dy=0.0)
        res = l.intersection(ol)
        self.assertIsInstance(res, tuple)
        self.assertAlmostEqual(res[0], 1.0)
        self.assertAlmostEqual(res[1], 0.0)
        ol.move(dz=0.25,dy=0.0)
        res = l.intersection(ol)
        self.assertIsInstance(res, tuple)
        self.assertAlmostEqual(res[0], 1.5)
        self.assertAlmostEqual(res[1], 0.0)
        l = Line(Point(1.0, 0.0), Point(0.0, 0.0))
        res = l.intersection(ol)
        self.assertIsInstance(res, tuple)
        self.assertAlmostEqual(res[0], -1.5)
        self.assertAlmostEqual(res[1], 0.0)
        ol.move(dz = 0.0, dy = 0.5)
        res = l.intersection(ol)
        self.assertIsInstance(res, tuple)
        self.assertAlmostEqual(res[0], -1.5)
        self.assertAlmostEqual(res[1], -0.5)
        ol = Line(Point(1.25, 1.5), Point(1.25, -0.5))
        res = l.intersection(ol)
        self.assertIsInstance(res, tuple)
        self.assertAlmostEqual(res[0], -1.5)
        self.assertAlmostEqual(res[1], 0.5)

    def testIntersects(self):
        l = Line(Point(0.0, 0.0), Point(1.0, 1.0))
        ol = Line(Point(0.0, 1.0), Point(1.0, 0.0))
        self.assertTrue(l.intersects(ol))
        # Parallel lines
        l = Line(Point(0.0, 0.0), Point(0.0, 1.0))
        ol = Line(Point(1.0, 0.0), Point(1.0, 1.0))
        self.assertFalse(l.intersects(ol))
        # Aligned lines
        l = Line(Point(0.0, 0.0), Point(0.0, 1.0))
        ol = Line(Point(0.0, 1.0), Point(0.0, 0.0))
        self.assertFalse(l.intersects(ol))
        # Testing qsi sign
        l = Line(Point(0.0, 0.0), Point(1.0, 0.0))
        ol = Line(Point(0.25, -1.0), Point(0.25, 1.0))
        self.assertTrue(l.intersects(ol))
        ol.move(dz=0.25,dy=0.0)
        self.assertTrue(l.intersects(ol))
        ol.move(dz=0.25,dy=0.0)
        self.assertTrue(l.intersects(ol))
        ol.move(dz=0.25,dy=0.0)
        self.assertTrue(l.intersects(ol))
        ol.move(dz=0.25,dy=0.0)
        self.assertFalse(l.intersects(ol))
        l = Line(Point(1.0, 0.0), Point(0.0, 0.0))
        self.assertFalse(l.intersects(ol))
        ol.move(dz = 0.0, dy = 0.5)
        self.assertFalse(l.intersects(ol))
        ol = Line(Point(1.25, 1.5), Point(1.25, -0.5))
        self.assertFalse(l.intersects(ol))

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

    def testDataMethods(self):
        v = Vertex('first_name',0.0, -2.0)
        self.assertEqual(v.data, {'name': 'first_name', 'z': 0.0, 'y': -2.0})
        v.data = {'name':'second_name', 'z': -1.2, 'y': 2.3}
        self.assertEqual(v.data, {'name':'second_name', 'z':-1.2, 'y': 2.3})

    def testNameByIndex(self):
        self.assertEqual(Vertex.__name_by_index__(0),'A')
        self.assertEqual(Vertex.__name_by_index__(1),'B')
        self.assertEqual(Vertex.__name_by_index__(2),'C')
        self.assertEqual(Vertex.__name_by_index__(24),'Y')
        self.assertEqual(Vertex.__name_by_index__(25),'Z')
        self.assertEqual(Vertex.__name_by_index__(26),'AA')
        self.assertEqual(Vertex.__name_by_index__(27),'AB')
        self.assertEqual(Vertex.__name_by_index__(28),'AC')
        self.assertEqual(Vertex.__name_by_index__(50),'AY')
        self.assertEqual(Vertex.__name_by_index__(51),'AZ')
        self.assertEqual(Vertex.__name_by_index__(52),'BA')


class TestSegment(unittest.TestCase):

    def testConstructor(self):
        lns = (('a',(0.0, 0.0), (1.0, 1.0)),
               ('3',(1.0, 1.0), (1.0, -1.0)),
               ('II',(1.0, -1.0), (-1.0, 1.0)),
               ('iv',(-1.0, 1.0), (0.0, 0.0)))
        for l in lns:
            nm, (zi, yi), (zj, yj) = l
            s = Segment(name = nm, pi = Point(zi, yi), pj = Point(zj, yj))
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
            s = Segment(name = nm, pi = Point(zi, yi), pj = Point(zj, yj))
            self.assertEqual(_str.format(nm, Point(zi, yi), Point(zj, yj)), str(s))

    def testThickness(self):
        s = Segment(name = 'A', pi = Point(0.0, 0.0), pj = Point(1.0, 1.0))
        self.assertIsNone(s.thickness)
        s.thickness = 1.0
        self.assertIsInstance(s.thickness, float)
        self.assertAlmostEqual(s.thickness, 1.0)
        s.thickness = None
        self.assertIsNone(s.thickness)

    def testVertices(self):
        s = Segment(name = 'A', pi = Point(0.0, 0.0), pj = Point(1.0, 1.0))
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
        s = Segment(name = 'A', pi = Point(0.0, 6.0), pj = Point(0.0, -6.0))
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
        s = Segment(name = 'A', pi = Point(0.0,6.0), pj = Point(0.0, -6.0))
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

    def testDataMethods(self):
        pi = Point(0.0, 0.0); pj = Point(1.0, -3.0)
        s = Segment(name='first_name',pi=pi, pj=pj)
        d = s.data
        self.assertEqual(d['name'],'first_name')
        self.assertEqual(d['pi'],pi)
        self.assertEqual(d['pj'],pj)
        self.assertIsNone(d['thickness'])
        s.thickness = 1.0

        d = s.data
        self.assertAlmostEqual(d['thickness'], 1.0)

        new_data = {}
        new_data['name'] = 'second_name'
        new_data['pi'] = Point(0.0, 0.0)
        new_data['pj'] = Point(0.0, 12.0)
        new_data['thickness'] = 0.1

        s.data = new_data

        d = s.data
        self.assertEqual(d['name'],'second_name')
        self.assertEqual(d['pi'],new_data['pi'])
        self.assertEqual(d['pj'],new_data['pj'])
        self.assertAlmostEqual(d['thickness'], 0.1)

    def testResultsMethods(self):
        pi = Point(0.0, 0.0); pj = Point(1.0, -3.0)
        s = Segment(name='first_name',pi=pi, pj=pj)
        rss = s.results
        self.assertIsNone(rss['first_vertex'])
        self.assertIsNone(rss['last_vertex'])
        self.assertEqual(rss['properties'],{})
        # Setting some 'results'
        fv = Vertex('A', 0.0, 0.0); lv = Vertex('B', 1.0, -3.0)
        s.first_vertex = fv
        s.last_vertex = lv
        d = s.results
        s.thickness = 0.1
        self.assertEqual(d['first_vertex'],fv)
        self.assertEqual(d['last_vertex'],lv)

        s.compute_properties()
        ppts = s.results['properties']
        self.assertEqual(len(ppts),8)
        new_rss = {}

        new_rss['first_vertex'] = Vertex('C', 0.0, 12.0)
        new_rss['last_vertex'] = Vertex('D', 0.0, 0.0)
        new_properties = {}
        new_properties['area'] = 1.2
        new_properties['Izl'] = 14.4
        new_properties['Iyl'] = 0.001
        new_properties['CG'] = Point(0.0, 6.0)
        new_properties['Iz'] = 14.4
        new_properties['Iy'] = 0.001
        new_properties['Izy'] = 0.0
        new_rss['properties'] = new_properties
        s.results = new_rss
        rss = s.results

        self.assertEqual(rss['first_vertex'],new_rss['first_vertex'])
        self.assertEqual(rss['last_vertex'],new_rss['last_vertex'])
        self.assertEqual(rss['properties'],new_properties)


class TestCrossSection(unittest.TestCase):

    def testConstructor(self):
        cs = Cross_section()
        self.assertIsInstance(cs, Cross_section)
        self.assertIsInstance(cs.segments, list)
        self.assertIsInstance(cs.vertices, dict)
        self.assertIsInstance(cs.results, dict)

    def testInsertAndRemoveSegment(self):
        cs = Cross_section()
        self.assertEqual(len(cs.segments), 0)
        s = Segment(Point(0.0, 0.0), Point(0.0, 6.0))
        cs.insert_segment(s)
        self.assertEqual(len(cs.segments), 1)
        self.assertIs(cs.segments[0],s)
        cs.remove_segment(s)
        self.assertEqual(len(cs.segments), 0)

    def testCreateSegment(self):
        cs = Cross_section()
        zi_ = 0.0; yi_ = 0.0
        zj_ = 0.0; yj_ = 6.0
        tks = 0.2
        cs.create_segment(zi = zi_, yi = yi_, zj = zj_, yj = yj_,
                           thickness = tks)
        s = Segment(pi = Point(zi_, yi_), pj = Point(zj_, yj_), thickness=tks)
        self.assertEqual(cs.segments[0].pi, s.pi)
        self.assertEqual(cs.segments[0].pj, s.pj)
        self.assertAlmostEqual(cs.segments[0].thickness, s.thickness)
        self.assertNotEqual(cs.segments[0],s)
        # Remove will not work because they are equal but not the same segment
        cs.remove_segment(s)
        self.assertEqual(len(cs.segments),1)

    def testGetVerticesByCoordinate(self):
        cs = Cross_section()
        cs.vertices['A'] = Vertex(name='A',z = 0.0, y = 0.0)
        cs.vertices['B'] = Vertex(name='B',z = 1.0, y = 1.0)
        cs.vertices['C'] = Vertex(name='C',z = 1.0, y = 1.0)
        cs.vertices['D'] = Vertex(name='D',z = 1.0, y = -1.0)
        vs = cs.get_vertices_by_coordinate(point = Point(0.0, 0.0),
                                           distance = 0.01)
        self.assertEqual(len(vs), 1)
        self.assertIn(cs.vertices['A'],vs)
        cs.vertices['A'].move(0.0,0.009)
        vs = cs.get_vertices_by_coordinate(point = Point(0.0, 0.0),
                                           distance = 0.01)
        self.assertEqual(len(vs), 1)
        self.assertIn(cs.vertices['A'],vs)
        vs = cs.get_vertices_by_coordinate(point = Point(0.0, -0.02),
                                           distance = 0.01)
        self.assertEqual(len(vs), 0)
        vs = cs.get_vertices_by_coordinate(point = Point(1.001, 1.001),
                                           distance = 0.01)
        self.assertEqual(len(vs), 2)
        self.assertIn(cs.vertices['B'],vs)
        self.assertIn(cs.vertices['C'],vs)
        vs = cs.get_vertices_by_coordinate(point = Point(1.01, 1.01),
                                           distance = 0.01)
        self.assertEqual(len(vs),0)
        vs = cs.get_vertices_by_coordinate(point = Point(1.001, -0.999),
                                           distance = 0.01)
        self.assertEqual(len(vs), 1)
        self.assertIn(cs.vertices['D'],vs)

    def testGetVerticesByArea(self):
        cs = Cross_section()
        cs.vertices['A'] = Vertex(name='A',z = 1.0, y = 1.0)
        cs.vertices['B'] = Vertex(name='B',z = 2.0, y = 1.0)
        cs.vertices['C'] = Vertex(name='C',z = 1.0, y = 2.0)
        cs.vertices['D'] = Vertex(name='D',z = 1.0, y = -1.0)
        cs.vertices['E'] = Vertex(name='E',z = -1.0, y = -2.0)
        cs.vertices['F'] = Vertex(name='F',z = -1.0, y = -1.0)
        cs.vertices['G'] = Vertex(name='G',z = -2.0, y = -1.0)
        cs.vertices['H'] = Vertex(name='H',z = -1.0, y = 1.0)
        vs = cs.get_vertices_by_area(box=((1.0,-1.0),(-1.0,1.0)),tolerance=0.01)
        self.assertEqual(len(vs), 4)
        self.assertIn(cs.vertices['A'],vs)
        self.assertIn(cs.vertices['D'],vs)
        self.assertIn(cs.vertices['F'],vs)
        self.assertIn(cs.vertices['H'],vs)
        vs = cs.get_vertices_by_area(box=(Point(1.0,-1.0),(-1.0,1.0)),tolerance=0.01)
        self.assertIn(cs.vertices['A'],vs)
        self.assertIn(cs.vertices['D'],vs)
        self.assertIn(cs.vertices['F'],vs)
        self.assertIn(cs.vertices['H'],vs)
        vs = cs.get_vertices_by_area(box=((1.0,-1.0),Point(-1.0,1.0)),tolerance=0.01)
        self.assertIn(cs.vertices['A'],vs)
        self.assertIn(cs.vertices['D'],vs)
        self.assertIn(cs.vertices['F'],vs)
        self.assertIn(cs.vertices['H'],vs)
        vs = cs.get_vertices_by_area(box=(Point(1.0,-1.0),Point(-1.0,1.0)),tolerance=0.01)
        self.assertIn(cs.vertices['A'],vs)
        self.assertIn(cs.vertices['D'],vs)
        self.assertIn(cs.vertices['F'],vs)
        self.assertIn(cs.vertices['H'],vs)
        vs = cs.get_vertices_by_area(box=((0.0, 0.0),Point(2.0,2.0)),tolerance=0.01)
        self.assertEqual(len(vs), 3)
        self.assertIn(cs.vertices['A'],vs)
        self.assertIn(cs.vertices['B'],vs)
        self.assertIn(cs.vertices['C'],vs)
        vs = cs.get_vertices_by_area(box=((0.0, 0.0),Point(-2.0,-2.0)),tolerance=0.01)
        self.assertEqual(len(vs), 3)
        self.assertIn(cs.vertices['E'],vs)
        self.assertIn(cs.vertices['F'],vs)
        self.assertIn(cs.vertices['G'],vs)
        vs = cs.get_vertices_by_area(box=((1.0, -1.0),Point(0.0, 0.0)),tolerance=0.01)
        self.assertEqual(len(vs), 1)
        self.assertIn(cs.vertices['D'],vs)
        vs = cs.get_vertices_by_area(box=((-1.0, 1.0),Point(0.0, 0.0)),tolerance=0.01)
        self.assertEqual(len(vs), 1)
        self.assertIn(cs.vertices['H'],vs)
        vs = cs.get_vertices_by_area(box=((-10.0, 10.0),Point(10.0, -10.0)),tolerance=0.01)
        self.assertEqual(len(vs), 8)
        vs = cs.get_vertices_by_area(box=((-0.1, 0.1),Point(0.0, -0.0)),tolerance=0.01)
        self.assertEqual(len(vs), 0)

    def testGetSegmentsByCoordinate(self):
        cs = Cross_section()
        cs.create_segment(zi = 0.0, yi = 0.0, zj = 1.0, yj = 1.0, thickness = 0.1)
        cs.create_segment(zi = 1.0, yi = -1.0, zj = 1.0, yj = 1.0, thickness = 0.1)
        #cs.segments[0].name = 'A-B'
        #cs.segments[1].name = 'C-B'
        sgs = cs.get_segments_by_coordinate(Point(0.0,0.0),distance=0.01)
        self.assertEqual(len(sgs),1)
        self.assertIn(cs.segments[0],sgs)
        for v in (0.1, 0.25, 0.5, 0.75, 0.99):
            sgs = cs.get_segments_by_coordinate(Point(v,v),distance=0.01)
            self.assertEqual(len(sgs),1)
            self.assertIn(cs.segments[0],sgs)
        sgs = cs.get_segments_by_coordinate(Point(1.0,1.0),distance=0.01)
        self.assertEqual(len(sgs),2)
        self.assertIn(cs.segments[0],sgs)
        self.assertIn(cs.segments[1],sgs)
        sgs = cs.get_segments_by_coordinate(Point(1.01,1.01),distance=0.01)
        self.assertEqual(len(sgs),0)
        for v in (0.9, 0.5, 0.0, -0.5, -1.0):
            sgs = cs.get_segments_by_coordinate(Point(1.0,v),distance=0.01)
            self.assertEqual(len(sgs),1)
            self.assertIn(cs.segments[1],sgs)
        sgs = cs.get_segments_by_coordinate(Point(1.0,-1.02),distance=0.01)
        self.assertEqual(len(sgs),0)

    def testGetSegmentsByArea(self):
        cs = Cross_section()
        cs.create_segment(zi = 0.0, yi = 0.0, zj = 1.0, yj = 1.0, thickness = 0.1)
        cs.create_segment(zi = 1.0, yi = -1.0, zj = 1.0, yj = 1.0, thickness = 0.1)
        # Box selection - Only the first one
        sgs = cs.get_segments_by_area(box=((1.0, 1.0),(0.0, 0.0)))
        self.assertEqual(len(sgs),1)
        self.assertIn(cs.segments[0],sgs)
        # box selection - No segment
        sgs = cs.get_segments_by_area(box=((1.0, 1.0),(0.01, 0.0)))
        self.assertEqual(len(sgs),0)
        # Box selection - Only the second one
        sgs = cs.get_segments_by_area(box=((1.01, -1.01),(0.99, 1.01)))
        self.assertEqual(len(sgs),1)
        self.assertIn(cs.segments[1],sgs)
        # Cross selection - Only the first one
        sgs = cs.get_segments_by_area(box=((0.0, 0.0),(0.1, 0.1)))
        self.assertEqual(len(sgs),1)
        self.assertIn(cs.segments[0],sgs)
        sgs = cs.get_segments_by_area(box=((0.1, 0.1),(0.11, 1.11)))
        self.assertEqual(len(sgs),1)
        self.assertIn(cs.segments[0],sgs)
        # Only the second one
        sgs = cs.get_segments_by_area(box=((0.99, -0.1),(01.01, 0.1)))
        self.assertEqual(len(sgs),1)
        self.assertIn(cs.segments[1],sgs)
        # Select both of them
        sgs = cs.get_segments_by_area(box=((0.0, 0.0),(1.1, 1.1)))
        self.assertEqual(len(sgs),2)
        self.assertIn(cs.segments[0],sgs)
        self.assertIn(cs.segments[1],sgs)

    def testReadFromCsv(self):
        # Trying to read from a file that does not exist
        cs = Cross_section()
        msg = cs.read_from_csv('iofiles/nofile.csv')
        self.assertEqual(msg,'File not found! (iofiles/nofile.csv)')
        # Reading 'segments.csv' file
        msg = cs.read_from_csv('iofiles/segments.csv')
        self.assertEqual(msg,'3 segment(s) read from iofiles/segments.csv!')
        self.assertEqual(len(cs.segments),3)
        self.assertEqual(cs.segments[0].pi, Point(0.0, 0.0))
        self.assertEqual(cs.segments[0].pj, Point(1.0, 1.0))
        self.assertAlmostEqual(cs.segments[0].thickness, 0.1)
        self.assertEqual(cs.segments[1].pi, Point(1.0, 1.0))
        self.assertEqual(cs.segments[1].pj, Point(1.0, -1.0))
        self.assertAlmostEqual(cs.segments[1].thickness, 0.15)
        self.assertEqual(cs.segments[2].pi, Point(1.0, -1.0))
        self.assertEqual(cs.segments[2].pj, Point(2.0, -2.0))
        self.assertAlmostEqual(cs.segments[2].thickness, 0.12)
        # Reading 'segments-header.csv' file
        cs = Cross_section()
        msg = cs.read_from_csv('iofiles/segments-header.csv')
        self.assertEqual(msg,'3 segment(s) read from iofiles/segments-header.csv!')
        self.assertEqual(len(cs.segments),3)
        self.assertEqual(cs.segments[0].pi, Point(0.0, 0.0))
        self.assertEqual(cs.segments[0].pj, Point(1.0, 1.0))
        self.assertAlmostEqual(cs.segments[0].thickness, 0.1)
        self.assertEqual(cs.segments[1].pi, Point(1.0, 1.0))
        self.assertEqual(cs.segments[1].pj, Point(1.0, -1.0))
        self.assertAlmostEqual(cs.segments[1].thickness, 0.15)
        self.assertEqual(cs.segments[2].pi, Point(1.0, -1.0))
        self.assertEqual(cs.segments[2].pj, Point(2.0, -2.0))
        self.assertAlmostEqual(cs.segments[2].thickness, 0.12)
        # Reading 'segments-header.csv' file
        cs = Cross_section()
        msg = cs.read_from_csv('iofiles/segments-dot-and-comma.csv')
        self.assertEqual(msg,'3 segment(s) read from iofiles/segments-dot-and-comma.csv!')
        self.assertEqual(len(cs.segments),3)
        self.assertEqual(cs.segments[0].pi, Point(0.0, 0.0))
        self.assertEqual(cs.segments[0].pj, Point(1.0, 1.0))
        self.assertAlmostEqual(cs.segments[0].thickness, 0.1)
        self.assertEqual(cs.segments[1].pi, Point(1.0, 1.0))
        self.assertEqual(cs.segments[1].pj, Point(1.0, -1.0))
        self.assertAlmostEqual(cs.segments[1].thickness, 0.15)
        self.assertEqual(cs.segments[2].pi, Point(1.0, -1.0))
        self.assertEqual(cs.segments[2].pj, Point(2.0, -2.0))
        self.assertAlmostEqual(cs.segments[2].thickness, 0.12)
        # Reading 'segments-with-error.csv'
        cs = Cross_section()
        msg = cs.read_from_csv('iofiles/segments-with-error.csv')
        self.assertEqual(msg,'2 segment(s) read from iofiles/segments-with-error.csv!')
        self.assertEqual(len(cs.segments),2)
        self.assertEqual(cs.segments[0].pi, Point(0.0, 0.0))
        self.assertEqual(cs.segments[0].pj, Point(1.0, 1.0))
        self.assertAlmostEqual(cs.segments[0].thickness, 0.1)
        self.assertEqual(cs.segments[1].pi, Point(1.0, -1.0))
        self.assertEqual(cs.segments[1].pj, Point(2.0, -2.0))
        self.assertAlmostEqual(cs.segments[1].thickness, 0.12)

    def testDefineVertices(self):
        # First (simple) example
        cs = Cross_section()
        msg = cs.read_from_csv('iofiles/segments.csv')
        cs.define_vertices()
        self.assertEqual(len(cs.segments),3)
        self.assertEqual(len(cs.vertices),4)
        self.assertEqual(cs.vertices['A'],cs.segments[0].pi)
        self.assertEqual(cs.vertices['B'],cs.segments[0].pj)
        self.assertEqual(cs.vertices['B'],cs.segments[1].pi)
        self.assertEqual(cs.vertices['C'],cs.segments[1].pj)
        self.assertEqual(cs.vertices['C'],cs.segments[2].pi)
        self.assertEqual(cs.vertices['D'],cs.segments[2].pj)
        self.assertEqual(cs.vertices['A'],cs.segments[0].first_vertex)
        self.assertEqual(cs.vertices['B'],cs.segments[0].last_vertex)
        self.assertEqual(cs.vertices['B'],cs.segments[1].first_vertex)
        self.assertEqual(cs.vertices['C'],cs.segments[1].last_vertex)
        self.assertEqual(cs.vertices['C'],cs.segments[2].first_vertex)
        self.assertEqual(cs.vertices['D'],cs.segments[2].last_vertex)
        self.assertEqual(cs.segments[0].name,'A-B')
        self.assertEqual(cs.segments[1].name,'B-C')
        self.assertEqual(cs.segments[2].name,'C-D')
        # Second example
        cs = Cross_section()
        msg = cs.read_from_csv('iofiles/segment-branchs.csv')
        cs.define_vertices()
        self.assertEqual(len(cs.segments),7)
        self.assertEqual(len(cs.vertices),8)
        self.assertEqual(cs.vertices['A'],cs.segments[0].pi)
        self.assertEqual(cs.vertices['B'],cs.segments[0].pj)
        self.assertEqual(cs.vertices['C'],cs.segments[1].pi)
        self.assertEqual(cs.vertices['D'],cs.segments[1].pj)
        self.assertEqual(cs.vertices['C'],cs.segments[2].pi)
        self.assertEqual(cs.vertices['B'],cs.segments[2].pj)
        self.assertEqual(cs.vertices['E'],cs.segments[3].pi)
        self.assertEqual(cs.vertices['B'],cs.segments[3].pj)
        self.assertEqual(cs.vertices['A'],cs.segments[4].pi)
        self.assertEqual(cs.vertices['F'],cs.segments[4].pj)
        self.assertEqual(cs.vertices['G'],cs.segments[5].pi)
        self.assertEqual(cs.vertices['B'],cs.segments[5].pj)
        self.assertEqual(cs.vertices['H'],cs.segments[6].pi)
        self.assertEqual(cs.vertices['A'],cs.segments[6].pj)
        self.assertEqual(cs.vertices['A'],cs.segments[0].first_vertex)
        self.assertEqual(cs.vertices['B'],cs.segments[0].last_vertex)
        self.assertEqual(cs.vertices['C'],cs.segments[1].first_vertex)
        self.assertEqual(cs.vertices['D'],cs.segments[1].last_vertex)
        self.assertEqual(cs.vertices['B'],cs.segments[2].first_vertex)
        self.assertEqual(cs.vertices['C'],cs.segments[2].last_vertex)
        self.assertEqual(cs.vertices['B'],cs.segments[3].first_vertex)
        self.assertEqual(cs.vertices['E'],cs.segments[3].last_vertex)
        self.assertEqual(cs.vertices['A'],cs.segments[4].first_vertex)
        self.assertEqual(cs.vertices['F'],cs.segments[4].last_vertex)
        self.assertEqual(cs.vertices['B'],cs.segments[5].first_vertex)
        self.assertEqual(cs.vertices['G'],cs.segments[5].last_vertex)
        self.assertEqual(cs.vertices['A'],cs.segments[6].first_vertex)
        self.assertEqual(cs.vertices['H'],cs.segments[6].last_vertex)
        self.assertEqual(cs.segments[0].name,'A-B')
        self.assertEqual(cs.segments[1].name,'C-D')
        self.assertEqual(cs.segments[2].name,'B-C')
        self.assertEqual(cs.segments[3].name,'B-E')
        self.assertEqual(cs.segments[4].name,'A-F')
        self.assertEqual(cs.segments[5].name,'B-G')
        self.assertEqual(cs.segments[6].name,'A-H')
    
    def testComputeProperties(self):
        # First (simple) example
        cs = Cross_section()
        cs.create_segment(zi = 0.0, yi = -6.0, zj = 0.0, yj = 6.0, thickness= 1.0)
        cs.define_vertices()
        self.assertEqual(len(cs.vertices),2)
        self.assertEqual(cs.vertices['A'],cs.segments[0].pi)
        self.assertEqual(cs.vertices['B'],cs.segments[0].pj)
        cs.compute_properties()
        p = cs.results
        self.assertAlmostEqual(p['A'],12.0)
        self.assertAlmostEqual(p['Iz'],144.0)
        self.assertAlmostEqual(p['Iy'],1.0)
        self.assertAlmostEqual(p['Izy'],0.0)
        self.assertAlmostEqual(p['I1'],144.0)
        self.assertAlmostEqual(p['I2'],1.0)
        self.assertAlmostEqual(p['It'],4.0)
        zcg, ycg = p['CG']
        self.assertAlmostEqual(zcg, 0.0)
        self.assertAlmostEqual(ycg, 0.0)
        self.assertAlmostEqual(p['PAA'], 0.0)
        angle = 0.0
        for _ in range(9):
            angle += 10.0
            cs.segments[0].rotate(Point(0.0, 0.0), 10.0, deg=True)
            cs.define_vertices()
            cs.compute_properties()
            p = cs.results
            self.assertAlmostEqual(p['A'],12.0)
            # self.assertAlmostEqual(p['Iz'],144.0)
            # self.assertAlmostEqual(p['Iy'],1.0)
            # self.assertAlmostEqual(p['Izy'],0.0)
            self.assertAlmostEqual(p['I1'],144.0)
            self.assertAlmostEqual(p['I2'],1.0)
            self.assertAlmostEqual(p['It'],4.0)        
            zcg, ycg = p['CG']
            self.assertAlmostEqual(zcg, 0.0)
            self.assertAlmostEqual(ycg, 0.0)
            self.assertAlmostEqual(p['PAA'], angle)
        # Cross section I
        cs = Cross_section()
        cs.create_segment(zi = 0.0, yi = 0.0, zj = 1.0, yj = 0.0, thickness = 0.2)
        cs.create_segment(zi = 1.0, yi = 0.0, zj = 2.0, yj = 0.0, thickness = 0.2)
        cs.create_segment(zi = 1.0, yi = 0.0, zj = 1.0, yj = 4.0, thickness = 0.1)
        cs.create_segment(zi = 0.0, yi = 4.0, zj = 1.0, yj = 4.0, thickness = 0.2)
        cs.create_segment(zi = 1.0, yi = 4.0, zj = 2.0, yj = 4.0, thickness = 0.2)
        cs.define_vertices()
        self.assertEqual(len(cs.vertices),6)
        self.assertEqual(cs.vertices['A'],cs.segments[0].pi)
        self.assertEqual(cs.vertices['B'],cs.segments[0].pj)
        self.assertEqual(cs.vertices['B'],cs.segments[1].pi)
        self.assertEqual(cs.vertices['C'],cs.segments[1].pj)
        self.assertEqual(cs.vertices['B'],cs.segments[2].pi)
        self.assertEqual(cs.vertices['D'],cs.segments[2].pj)
        self.assertEqual(cs.vertices['E'],cs.segments[3].pi)
        self.assertEqual(cs.vertices['D'],cs.segments[3].pj)
        self.assertEqual(cs.vertices['D'],cs.segments[4].pi)
        self.assertEqual(cs.vertices['F'],cs.segments[4].pj)
        zcg_, ycg_ = 1.0, 2.0
        for dz, dy in ((0.0, 0.0), (-2.0, 0.0), (0.0, -4.0), (0.0, 4.0)):
            for segment in cs.segments:
                segment.move(dz,dy)
            zcg_ += dz; ycg_ += dy
            cs.define_vertices()
            cs.compute_properties()
            p = cs.results
            self.assertAlmostEqual(p['A'],1.2)
            self.assertAlmostEqual(p['Iz'],3.736)
            self.assertAlmostEqual(p['Iy'],0.267)
            self.assertAlmostEqual(p['Izy'],0.0)
            self.assertAlmostEqual(p['I1'],3.736)
            self.assertAlmostEqual(p['I2'],0.267)
            self.assertAlmostEqual(p['It'],0.012)
            zcg, ycg = p['CG']
            self.assertAlmostEqual(zcg, zcg_)
            self.assertAlmostEqual(ycg, ycg_)
            self.assertAlmostEqual(p['PAA'], 0.0)
        ang = 0.0
        dang = 10.0
        for _ in range(9):
            ang += dang
            for segment in cs.segments:
                segment.rotate(Point(zcg_, ycg_), dang)
            cs.compute_properties()
            p = cs.results
            self.assertAlmostEqual(p['A'],1.2)
            # self.assertAlmostEqual(p['Iz'],3.736)
            # self.assertAlmostEqual(p['Iy'],0.267)
            # self.assertAlmostEqual(p['Izy'],0.0)
            self.assertAlmostEqual(p['I1'],3.736)
            self.assertAlmostEqual(p['I2'],0.267)
            self.assertAlmostEqual(p['It'],0.012)
            zcg, ycg = p['CG']
            self.assertAlmostEqual(zcg, zcg_)
            self.assertAlmostEqual(ycg, ycg_)
            self.assertAlmostEqual(p['PAA'], ang)
        cs.compute_properties()
        cs.compute_sectorial_area_and_torsion_center()

    def testComputeAllResults(self):
        # U Section (book, page 67)
        cs = Cross_section()
        cs.create_segment(zi = 0.0, yi = 0.0, zj = 10.0, yj = 0.0, thickness = 1.0)
        cs.create_segment(zi = 0.0, yi = 10.0, zj = 00.0, yj = 0.0, thickness = 1.0)
        cs.create_segment(zi = 10.0, yi = 0.0, zj = 10.0, yj = 10.0, thickness = 1.0)
        cs.define_vertices()
        self.assertEqual(len(cs.vertices),4)
        self.assertEqual(cs.vertices['A'],cs.segments[0].pi)
        self.assertEqual(cs.vertices['B'],cs.segments[0].pj)
        self.assertEqual(cs.vertices['C'],cs.segments[1].pi)
        self.assertEqual(cs.vertices['D'],cs.segments[2].pj)
        cs.compute_properties()
        cs.compute_sectorial_area_and_torsion_center()
        res = cs.results
        self.assertAlmostEqual(res['A'],30.0)
        self.assertAlmostEqual(res['Iz'],334.16666666)
        self.assertAlmostEqual(res['Iy'],585.0)
        self.assertAlmostEqual(res['Izy'],0.0)
        self.assertAlmostEqual(res['I1'],585.0)
        self.assertAlmostEqual(res['I2'],334.16666666)
        self.assertAlmostEqual(res['It'],10.0)
        zcg, ycg = res['CG']
        self.assertAlmostEqual(zcg, 5.0)
        self.assertAlmostEqual(ycg, 3.333333333)
        self.assertAlmostEqual(res['PAA'], 90.0)
        zct, yct = res['CT']
        self.assertAlmostEqual(zct, 5.0)
        self.assertAlmostEqual(yct, -4.26400759)
        ws = res['ws']
        self.assertAlmostEqual(ws['A'], 21.32003798)
        self.assertAlmostEqual(ws['B'], -21.32003798)
        self.assertAlmostEqual(ws['C'], -28.679962013)
        self.assertAlmostEqual(ws['D'], 28.679962013)
        self.assertAlmostEqual(res['Iw'], 5952.65580756)
        # book example, page 70
        cs = Cross_section()
        cs.create_segment(zi = -10.0, yi = -10.0, zj = 0.0, yj = -10.0, thickness = 1.0)
        cs.create_segment(zi = 0.0, yi = -10.0, zj = 0.0, yj = 10.0, thickness = 1.0)
        cs.create_segment(zi = 0.0, yi = 10.0, zj = 10.0, yj = 10.0, thickness = 1.0)
        cs.define_vertices()
        self.assertEqual(len(cs.vertices),4)
        self.assertEqual(cs.vertices['A'],cs.segments[0].pi)
        self.assertEqual(cs.vertices['B'],cs.segments[0].pj)
        self.assertEqual(cs.vertices['C'],cs.segments[1].pj)
        self.assertEqual(cs.vertices['D'],cs.segments[2].pj)
        cs.compute_properties()
        cs.compute_sectorial_area_and_torsion_center()
        res = cs.results
        self.assertAlmostEqual(res['A'],40.0)
        self.assertAlmostEqual(res['Iz'],2668.3333333)
        self.assertAlmostEqual(res['Iy'],668.3333333)
        self.assertAlmostEqual(res['Izy'],1000.0)
        self.assertAlmostEqual(res['I1'],3082.5468957)
        self.assertAlmostEqual(res['I2'],254.11977096)
        self.assertAlmostEqual(res['It'],13.3333333)
        zcg, ycg = res['CG']
        self.assertAlmostEqual(zcg, 0.0)
        self.assertAlmostEqual(ycg, 0.0)
        self.assertAlmostEqual(res['PAA'], 22.5)
        zct, yct = res['CT']
        self.assertAlmostEqual(zct, 0.0)
        self.assertAlmostEqual(yct, 0.0)
        ws = res['ws']
        self.assertAlmostEqual(ws['A'], -75.0)
        self.assertAlmostEqual(ws['B'], 25.0)
        self.assertAlmostEqual(ws['C'], 25.0)
        self.assertAlmostEqual(ws['D'], -75.0)
        self.assertAlmostEqual(res['Iw'], 41666.66666666)

if __name__ == "__main__":
    unittest.main()
