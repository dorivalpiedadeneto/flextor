# Remarks from previous tests

2023-12-28

- Always remember: code according to the MVP philosophy
- Simple (very simple) class design (MVP)
- Only straight lines (at least in the first release version)
- Class Point - index (int or None), y and z (float)
- Class Line - index (int or None), start and end points (pi and pj)
- Class Segment - derived from line, includes additional data: t (thickness), vertices and other properties (e.g. Iy, Iz, Ixy, I1, I2, ...)
- Class Vertex - derived from point, also stores the coordinate for the principal axes
- Class Cross_section: includes data as segments, vertexes, CG (point), D (point), P (point), rotation of principal axes, ...