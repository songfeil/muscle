# Geometry Processing – Muscle Mesh Project

## Compilation

Starting in this directory, issue:

    mkdir build
    cd build
    cmake ..
    make 

## Execution

Once built, you can execute the program from inside the `build/` using 

    ./muscle
    
## Documentation

An overview of each header file, documenting each function.

### Muscle and Tendon Generation

...

### User Input

Much of the user interface functionality was directly coded into the main.cpp file, however some repeatedly used or more complex chunks of code were wrapped into their own functions and files.

#### ui_helpers.h
A variety of functions that are used for user-input.

- intersection_with_xy_plane()
    - User input currently defaults to the xy plane for simplification in 3D-- this function returns an intersection point in 3D using raycasting.
- add_face_to_patch()
    - Used to select faces on the bone mesh, and group them into "patches" if they share an edge (at the time the face is selected).
- verts_within_x_range()
    - Used for the experimental inflation/deflation mode, highlights vertices within the x-range of some point (the cursor)

#### mesh_editing.h
- xflate_mesh()
    - Either inflate or deflate a whole given mesh, by moving all vertices in the direction of the normal.
- xflate_verts()
    - Either inflate or deflate only a subset of vertices on a mesh.
- xflate_verts_in_xrange()
    - Either inflate or deflate verts within the x-range of a given point.
- smooth_mesh_with_fixed()
    - Smooth mesh but pin fixed verts back into place (work-around of non-volume-preserving smoothing)

#### generate_bone.h
- start_flinstone_bone()
    - Given a point, generate new "flinstone bone" mesh with one end pinned to that point.
- transform_flinstone_bone()
    - Given a start point and an end point, transform a given "flinstone bone" mesh lengthwise along the vector between the points
    

#### generate_bone.h
