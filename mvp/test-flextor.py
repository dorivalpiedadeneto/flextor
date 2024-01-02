import unittest

from flextor import Point
from flextor import Line

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


            




    

if __name__ == "__main__":
    unittest.main()