# Minimal class (no setters/getters, no typing, no tests, minimal comments...)

from math import sin, cos, radians, sqrt, atan2, degrees
from os import path

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

    @property
    def data(self):
        return {'z': self.z, 'y': self.y}
    
    @data.setter
    def data(self, value):
        self.z = value['z']
        self.y = value['y']

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

    @property
    def data(self):
        return {'pi': self.pi.data, 'pj': self.pj.data}

    @data.setter
    def data(self, value):
        self.pi.data = value['pi']
        self.pj.data = value['pj']

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
        # length from middle pointo to pj
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

    def intersection(self, other_line, tolerance = 1.0e-6):
        '''
        Returns a tuple with the qsis of the intersection of the two lines
        (qsi, qsi_), qsi -> from self, qsi_ -> other line, or None if
        lines are collinear or parallel.
        '''
        zi, yi = self.pi.coord()
        zj, yj = self.pj.coord()
        zi_, yi_ = other_line.pi.coord()
        zj_, yj_ = other_line.pj.coord()
        A11 = 0.5 * (zj - zi); A12 = 0.5 * (zi_ - zj_)
        A21 = 0.5 * (yj - yi); A22 = 0.5 * (yi_ - yj_)
        det = A11 * A22 - A12 * A21
        if abs(det) < tolerance:
            return None
        b1 = 0.5 * (zi_ + zj_) - 0.5 * (zi + zj)
        b2 = 0.5 * (yi_ + yj_) - 0.5 * (yi + yj)
        qsi = (A22 * b1 - A12 * b2) / det
        qsi_ = (-A21 * b1 + A11 * b2) / det
        return (qsi, qsi_)

    def intersects(self, other_line, tolerance=1.0e-6):
        res = self.intersection(other_line, tolerance)
        if not res:
            return False
        qsi, qsi_ = res
        inf = -1.0 - tolerance
        sup = 1.0 + tolerance
        if (inf < qsi < sup) and (inf < qsi_ < sup):
            return True
        else:
            return False

class Vertex(Point):
    _name_letters = 'ABCDEFGHIJKLMNOPQRTSUVWXYZ'
    _nl = len(_name_letters)

    def __init__(self, name, z, y):
        self.name = name        # planning to use string to make reference
        super().__init__(z, y)

    def __repr__(self):
        return "vertex {} at z = {:.3f}, y = {:.3f}".format(self.name, self.z, self.y)
    
    def __str__(self):
        return "vertex {} at z = {:.3f}, y = {:.3f}".format(self.name, self.z, self.y)

    @property
    def data(self):
        return {'name': self.name, 'z': self.z, 'y': self.y}

    @data.setter
    def data(self, value):
        self.name = value['name']
        self.z = value['z']
        self.y = value['y']

    @classmethod
    def __name_by_index__(cls, index):
        if index // Vertex._nl:
            name = Vertex._name_letters[(index // Vertex._nl)-1] + \
            Vertex._name_letters[index % Vertex._nl]
        else:
            name = Vertex._name_letters[index % Vertex._nl]
        return name

class Segment(Line):

    def __init__(self, pi, pj, name = '', thickness = None):
        self.name = name
        super().__init__(pi, pj)
        self.thickness = thickness
        self.first_vertex = None
        self.last_vertex = None
        self.properties = {}

    @property
    def data(self):
        return {'name': self.name, 'pi': self.pi, 'pj': self.pj,
                'thickness': self.thickness}

    @data.setter
    def data(self, value):
        self.name = value['name']
        self.pi = value['pi']
        self.pj = value['pj']
        self.thickness = value['thickness']

    @property
    def results(self):
        return {'first_vertex': self.first_vertex, 
                'last_vertex': self.last_vertex,
                'properties': self.properties}        

    @results.setter
    def results(self, value):
        self.first_vertex = value['first_vertex']
        self.last_vertex = value['last_vertex']
        self.properties = value['properties']

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

    def compute_properties(self):
        if self.thickness:
            t = self.thickness
            l = self.length()
            zi, yi = self.pi.coord()
            zj, yj = self.pj.coord()
            zcg = 0.5 * (zi + zj)
            ycg = 0.5 * (yi + yj)
            A = t * l
            Izl = t * l ** 3 / 12.0     # local z axes
            Iyl = t ** 3 * l /12.0      # local y axes
            It = l * t ** 3 / 3.0
            s, c = self.unit_tangent()
            # Remember that Izyl = 0
            Iz = c**2 * Izl + s**2 * Iyl # - 2.0 * Izyl * c * s 
            Iy = s**2 * Izl + c**2 * Iyl # + 2.0 * Izyl * c * s 
            Izy = (Izl - Iyl) * c * s # + Izyl * (c**2 - s**2)
            self.properties['area'] = A
            self.properties['Izl'] = Izl
            self.properties['Iyl'] = Iyl
            self.properties['CG'] = Point(zcg, ycg)
            self.properties['Iz'] = Iz
            self.properties['Iy'] = Iy
            self.properties['Izy'] = Izy
            self.properties['It'] = It

    def sectorial_area_contribution(self, center):
        # Positive in clockwise direction
        if self.first_vertex and self.last_vertex:
            zc, yc = center.coord()
            zi, yi = self.first_vertex.coord()
            zj, yj = self.last_vertex.coord()
            w = zc * (yj - yi) + zj * (yi - yc) + zi * (yc - yj)
        else:
            w = None
        return w

class Cross_section(object):

    def __init__(self):
        self.segments = []      # list of segment instances
        self.vertices = {}      # dict of vertices instances
        self.results = {}       # results to be stored
        # Results to be stored
        # CG - Geometric Center
        # CT - Torsion Center
        # PPA - Principal Axes Angle
        # Iz - Moment of Inertia - z axis
        # Iy - Moment of Inertia - y axis
        # Iyz - Product of Inertia - z and y axes
        # I1 - First Principal Moment of Inertia
        # I2 - Second Principal Moment of Inertia
        # It - Polar Moment of Inertia
        # Iw - Sectorial Moment of Inertia
        # ws - List os sectorial areas (for each segment)

    def insert_segment(self, segment):
        self.segments.append(segment)

    def remove_segment(self, segment):
        if segment in self.segments:
            self.segments.remove(segment)

    def create_segment(self, zi, yi, zj, yj, thickness = None):
        s = Segment(Point(zi, yi), Point(zj, yj))
        if thickness:
            s.thickness = thickness
        self.segments.append(s)

    def get_vertices_by_coordinate(self, point, distance):
        vertices = []
        for vertex in self.vertices.values():
            if point.distance(vertex) <= distance:
                vertices.append(vertex)
        return vertices

    def get_vertices_by_area(self, box, tolerance):
        vertices = []
        if isinstance(box[0], Point):
            zi, yi = box[0].coord()
        else: # implied it is a list or tuple (with two values)
            zi, yi = box[0]
        if isinstance(box[1], Point):
            zj, yj = box[1].coord()
        else: # implied it is a list or tuple (with two values)
            zj, yj = box[1]
        zmin = min(zi, zj) - tolerance
        zmax = max(zi, zj) + tolerance
        ymin = min(yi, yj) - tolerance
        ymax = max(yi, yj) + tolerance
        for vertex in self.vertices.values():
            z, y = vertex.coord()
            if (zmin <= z <= zmax) and (ymin <= y <= ymax):
                vertices.append(vertex)
        return vertices

    def get_segments_by_coordinate(self, point, distance):
        segments = []
        for segment in self.segments:
            if segment.distance(point) <= distance:
                segments.append(segment)
        return segments

    def get_segments_by_area(self, box, tolerance=1.0e-6):
        segments = []
        if isinstance(box[0], Point):
            zi, yi = box[0].coord()
        else: # implied it is a list or tuple (with two values)
            zi, yi = box[0]
        if isinstance(box[1], Point):
            zj, yj = box[1].coord()
        else: # implied it is a list or tuple (with two values)
            zj, yj = box[1]
        zmin = min(zi, zj) - tolerance
        zmax = max(zi, zj) + tolerance
        ymin = min(yi, yj) - tolerance
        ymax = max(yi, yj) + tolerance
        # right to left - crossing selection - selects any object that either
        # crosses the boundary or is inside it
        # left to right - box selection - select on objects that are completely
        # within the box
        if zi < zj: # Cross selection
            for segment in self.segments:
                zsi, ysi = segment.pi.coord()
                zsj, ysj = segment.pj.coord()
                # If is totally inside, is selected (as in box selection)
                if (zmin <= zsi <= zmax) and (ymin <= ysi <= ymax) and \
                   (zmin <= zsj <= zmax) and (ymin <= ysj <= ymax):
                    segments.append(segment)
                else:
                    lines = []
                    lines.append(Line(Point(zmin, ymin), Point(zmax, ymin)))
                    lines.append(Line(Point(zmin, ymax), Point(zmax, ymax)))
                    lines.append(Line(Point(zmin, ymin), Point(zmin, ymax)))
                    lines.append(Line(Point(zmax, ymin), Point(zmax, ymax)))
                    for line in lines:
                        if segment.intersects(line):
                            segments.append(segment)
                            break


        else: # Box selection
            for segment in self.segments:
                zsi, ysi = segment.pi.coord()
                zsj, ysj = segment.pj.coord()
                if (zmin <= zsi <= zmax) and (ymin <= ysi <= ymax) and \
                   (zmin <= zsj <= zmax) and (ymin <= ysj <= ymax):
                    segments.append(segment)
        return segments

    def read_from_csv(self, filename):
        # eol -> '\r\n': windows, '\n': linux
        # sep -> ';' or ',' (spaces or tabs not valid)
        # dec -> ',' or '.'
        if path.isfile(filename):
            try:
                f = open(filename, 'r')
                data = f.read()
                f.close()
            except:
                return 'Error reading {} file!'.format(filename)
        else:
            return 'File not found! ({})'.format(filename)
        if '\r\n' in data:
            eol = '\r\n'
        else:
            eol = '\n'
        if ';' in data:
            sep = ';'
            if ',' in data:
                data = data.replace(',','.')
        else:
            sep = ','
        count = 0
        for line in data.split(eol):
            values = line.split(sep)
            if len(values) == 5:
                try:
                    zi, yi, zj, yj, tks = values
                    zi = float(zi)
                    yi = float(yi)
                    zj = float(zj)
                    yj = float(yj)
                    tks = float(tks)
                    self.create_segment(zi, yi, zj, yj, tks)
                    count += 1
                except:
                    pass
        return '{} segment(s) read from {}!'.format(count, filename)

    def define_vertices(self):       
        vid = 0     # Vertex index
        # Reset old data
        self.vertices = {}
        for segment in self.segments:
            segment.name = ''
            segment.first_vertex = None
            segment.last_vertex = None
        # Loop all segments
        for segment in self.segments:
            pi = segment.pi
            pj = segment.pj
            found_pi = False
            found_pj = False
            for vertex in self.vertices.values():
                if pi == vertex:
                    found_pi = True
                if pj == vertex:
                    found_pj = True
            if not found_pi:
                vname = Vertex.__name_by_index__(vid)
                vid += 1
                self.vertices[vname] = Vertex(name = vname, z = pi.z, y = pi.y)
            if not found_pj:
                vname = Vertex.__name_by_index__(vid)
                vid += 1
                self.vertices[vname] = Vertex(name = vname, z = pj.z, y = pj.y)
        for vertex in self.vertices.values():
            for segment in self.segments:
                if not segment.first_vertex:
                    if segment.pi == vertex:
                        segment.first_vertex = vertex
                    elif segment.pj == vertex:
                        segment.first_vertex = vertex
                elif not segment.last_vertex:
                    if segment.pi == vertex:
                        segment.last_vertex = vertex
                    elif segment.pj == vertex:
                        segment.last_vertex = vertex
        for segment in self.segments:
            segment.name = segment.first_vertex.name + '-' + \
                           segment.last_vertex.name

    def compute_properties(self):
        # A, Iz, Iy, Izy, I1, I2, It e principal axis angle
        A = 0.0; Az = 0.0; Ay = 0.0
        for segment in self.segments:
            segment.compute_properties()
            prop = segment.properties
            area = prop['area']
            zc, yc = prop['CG'].coord()
            A += area
            Az += area * zc
            Ay += area * yc
        zcg = Az / A
        ycg = Ay / A
        Iz = 0.0; Iy = 0.0; Izy = 0.0; It = 0.0
        for segment in self.segments:
            prop = segment.properties
            area = prop['area']
            zc, yc = prop['CG'].coord()
            Iz += prop['Iz']
            Iz += area * (yc - ycg) ** 2
            Iy += prop['Iy']
            Iy += area * (zc - zcg) ** 2
            Izy += prop['Izy']
            Izy += area * (yc - ycg) * (zc -zcg)
            It += prop['It']
        # I1, I2 and principal axis angle
        I1 = 0.5 * (Iz + Iy) + sqrt((0.5 * (Iz - Iy)) ** 2 + Izy ** 2)
        I2 = 0.5 * (Iz + Iy) - sqrt((0.5 * (Iz - Iy)) ** 2 + Izy ** 2)
        paa = 0.5 * degrees(atan2(2.0*Izy,Iz - Iy))
        self.results['A'] = A
        self.results['Iz'] = Iz
        self.results['Iy'] = Iy
        self.results['Izy'] = Izy
        self.results['It'] = It
        self.results['I1'] = I1
        self.results['I2'] = I2
        self.results['CG'] = (zcg, ycg)
        self.results['PAA'] = paa

    def compute_sectorial_area_and_torsion_center(self):
        # Provisory center and provisory sectorial area
        if 'CG' in self.results.keys():
            zp, yp = self.results['CG']
        else:
            return # Cannot continue
        PC = Point(zp, yp)
        if 'PAA' in self.results.keys():
            ang = radians(self.results['PAA'])
        else:
            return # Cannot continue
        # Sectorial areas (provisory)
        ws = {'A': 0.0}
        remaining_vertices = list(self.vertices.keys())
        remaining_vertices.remove('A')
        while len(remaining_vertices) > 0:
            last_size = len(remaining_vertices)
            for segment in self.segments:
                fvname = segment.first_vertex.name
                lvname = segment.last_vertex.name
                if fvname in ws.keys():
                    if not lvname in ws.keys():
                        zi, yi = segment.first_vertex.coord()
                        zj, yj = segment.last_vertex.coord()
                        w = segment.sectorial_area_contribution(PC)
                        ws[lvname] = ws[fvname] + w
                        remaining_vertices.remove(lvname)
            if len(remaining_vertices) == last_size:
                #Infinite loop while computing sectorial areas!
                return
        # At this point, the provisory sectorial areas are already computed
        vertices = {}   # coordinates in the principal axis coordinates
        for vk in self.vertices.keys():
            vertex = self.vertices[vk]
            zv, yv = vertex.coord()
            z_ = zv - zp
            y_ = yv - yp
            z = z_ * cos(ang) - y_ * sin(ang)
            y = y_ * cos(ang) + z_ * sin(ang)
            name = vertex.name
            vertices[name] = (z, y)
        # Compute the torsion center coordinate (CT)
        zd_ = 0.0; yd_ = 0.0
        for segment in self.segments:
            t = segment.thickness
            l = segment.length()
            vi = segment.first_vertex.name
            vj = segment.last_vertex.name
            zi, yi = vertices[vi]
            zj, yj = vertices[vj]
            wi = ws[vi]; wj = ws[vj]
            yd_ += t * l * (wi * (2.0 * zi + zj) + wj * (2.0 * zj + zi)) / 6.0
            zd_ -= t * l * (wi * (2.0 * yi + yj) + wj * (2.0 * yj + yi)) / 6.0
        I1 = self.results['I1']
        I2 = self.results['I2']
        yd_ /= I2
        zd_ /= I1
        # Back to glbal axis (z and y)
        yd = yd_ * cos(-ang) + zd_ * sin(-ang)
        zd = zd_ * cos(-ang) - yd_ * sin(-ang)
        yd += yp
        zd += zp
        self.results['CT'] = (zd, yd)
        # Computing principal sectorial area
        CT = Point(zd, yd)
        ws = {'A': 0.0}
        remaining_vertices = list(self.vertices.keys())
        remaining_vertices.remove('A')
        while len(remaining_vertices) > 0:
            last_size = len(remaining_vertices)
            for segment in self.segments:
                fvname = segment.first_vertex.name
                lvname = segment.last_vertex.name
                if fvname in ws.keys():
                    if not lvname in ws.keys():
                        zi, yi = segment.first_vertex.coord()
                        zj, yj = segment.last_vertex.coord()
                        w = segment.sectorial_area_contribution(CT)
                        ws[lvname] = ws[fvname] + w
                        remaining_vertices.remove(lvname)
            if len(remaining_vertices) == last_size:
                #Infinite loop while computing sectorial areas!
                return
        # Computing wc
        wc = 0.0
        for segment in self.segments:
            t = segment.thickness
            l = segment.length()
            vi = segment.first_vertex.name
            vj = segment.last_vertex.name
            wi = ws[vi]; wj = ws[vj]
            wc += t * l * (wi + wj) / 2.0
        wc /= self.results['A']
        # Finally computing the principal sectorial area values
        for k in ws:
            ws[k] -= wc
        self.results['ws'] = ws
        # Computing IW
        Iw = 0.0
        for segment in self.segments:
            t = segment.thickness
            l = segment.length()
            vi = segment.first_vertex.name
            vj = segment.last_vertex.name
            wi = ws[vi]; wj = ws[vj]
            Iw += l * t * (wi ** 2 + wi * wj + wj ** 2) / 3.0
        self.results['Iw'] = Iw

if __name__ == "__main__":
    pass