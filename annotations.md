# Remarks from previous tests

2023-12-28

- Always remember: code according to the MVP philosophy
- Simple (very simple) class design (MVP)
- Only straight lines (at least in the first release version)
- Class Point - z and y (float)
- Class Line - start and end points (pi and pj)
- Class Segment - derived from line, includes additional data: ~~index~~ (name), t (thickness), vertices and other properties (e.g. Iy, Iz, Ixy, I1, I2, ...)
- Class Vertex - derived from point, stores index ~~and the coordinate for the principal axes~~
- Class Cross_section: includes data as segments, vertexes, CG (point), D (point), P (point), rotation of principal axes, ...

2024-01-03

 - For later: methods to idenfify lines' intersection; to be used to trim or extend methods
 - ~~Maybe a class for representing the coordinate system (besides the global one)~~ (maybe this is too complicated; just store the principal coordinate system) information in the cross section instance, as well as the methods to compute the coordinates for the principal directions)
 - Using another coordinate system it is not necessary to store the coordinate for the principal axes

2024-01-04

 - Vertices are going to be created after all segments are included; the software will name them. So it makes sense to use string as their names (letters are better than number for plotting the results in this context). One vertex may be shared by two or more segments.
 - Segments are going to be created by including lines in the model. The segments name will also be created by the software.
 - Input by file using json: probably will indicate coordinates and thickness (no name - they will be given by the code);
 - Input by GUI: the coordinates may be given by click or coordinates input, with thickness value;
 - Input by dxf: user needs to inform a thickness for each layer (import by layer), or for all lines in file.
 - To select a line in GUI, user can use click or a selection box; after selected, the could be deleted, or change.
 - As there is no reference to the line id, maybe they can be stored in lists instead of dictionaries in the cross section instance,
 - Considering these cases, there is no need to inform names in the input file.
 - Output file: json; maybe it is a good idea to export the segments names and vertices, as well as the results. The software version either.
 - Next methods to be implemented: read from input files, export to output files.