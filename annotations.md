# Remarks from previous tests

2023-12-28

- Always remember: code according to the MVP philosophy
- Simple (very simple) class design (MVP)
- Only straight lines (at least in the first release version)
- Class Point - z and y (float)
- Class Line - start and end points (pi and pj)
- Class Segment - derived from line, includes additional data: ~~index~~ (name), t (thickness), vertices and other properties (e.g. Iy, Iz, Ixy, I1, I2, ...)
- Class Vertex - derived from point, stores index ~~and the coordinate for the principal axes~~
- Class Cross_section: includes data as segments, vertices, CG (point), D (point), P (point), rotation of principal axes, ...

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
  
2024-01-11

 - Segments have some parameters that may by None or 'nothing' (like an empty dictonary for instance). It would be a good idea that in these cases the json file can just omit these (and while reading or exporting these parameters are just omitted)?
   - Just name and points:
    ```{name: 'S01', 'pi': {'z': 0.0, 'y': 0.0}, 'pj': {'z': 0.0, 'y': 1.0}} ```
   - All of them (nothing omitted):
    ```{name: 'S01', 'pi': {'z': 0.0, 'y': 0.0}, 'pj': {'z': 0.0, 'y': 1.0}, 'thickness': null, 'first_vertex': null, 'last_vertex': null, 'properties':{}}```
   Another issue: How to store a vertex instance as json (it is only a reference to an existing vertex in the vertex list). Maybe, the json should only have the vertex name (an string). In this case, it appears that the best way to deal with this issues it defining an static method to deal with this. But the vertex is something full of the segment scope. Probably, the cross section class must be where this is dealt with.
 - Probably the above is not the best solution. Maybe, the data method must just deal with the input data (name, that can be an empty string - in  most cases it will be, in fact - the points and thickness, that also can be null). The vertices and properties are compute by the software, so they must be treated as a result, that may be available or not. This appears to be the best approach.
    ```{name: 'S01', 'pi': {'z': 0.0, 'y': 0.0}, 'pj': {'z': 0.0, 'y': 1.0}, 'thickness': null} ```
 - Then, results method (get and set) are also implemented in the segment class.
 - Results are going to be stored in a different dictionary in the input/output files, so there must be a way to easily find the segment in which the result (if available is going to be stored); maybe still keep it as a list, but use the elements' name to make reference to them (it is more confortable to keep it as a list);
 - The segment's name will be defined by the software (the user will not define it); probably it will be given considering the vertices names (A-B,A-D,F-C,...);
 - If no results are available, segments name will be set as '' (empty string)
 - So, lets' change the cross section class to store data as defined above.
 - Maybe results are not a good name for the vertices methods; maybe it is better to use the name 'properties' instead of 'results' (results is a better name to be used in the cross section class)
 - Implemented methods to insert and remove segment to cross section; also to create by giving coordinates.
 - In the sequence: create method to get segment by coordinate or by a selection box. In fact, to have the same behaviour as CADs for selection box (from left to right or right to left), it is necessary to implement another two methods in line's class (also it would be a good opportunity to implement a method in point class too, to select using selection boxes).
  
2024-01-12

 - At least for the current implementation, the cross section is the one that will manage the other classes. I believe that it is better to let it deal it in the cross section class.
 - For segment: right to left - crossing selection - selects any object that either crosses the boundary or is inside it
 - For segment: left to right - box selection - select on objects that are completely within the box
 - For inputing line (segments): probably the easiest way is to use a csv file! (just an idea)

2024-01-15

 - Cross selection demands finding line's intersection; so, write this methods first!
 - The derived formulation is wrong (when determinant is negative); needs reviewing of intersection method in line.

2024-01-22

 - Implement method to identify the vertices
 - First, define a class method in vertex to define a name (with uppercase letters) to the vertices (test it also)
  
2024-01-25

 - Implemented sectorial_area_contribution in vertex class (but not tested yet)
 - Started implementing compute_sectorial_area_and_torsion_center in cross section class (not finished yet)
 - Next step: Computing torsion center position