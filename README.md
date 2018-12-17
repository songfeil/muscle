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

#### generate_muscle.h
- generate_muslce
    - Given control points, generate Catmull Rom Spline curve, sample points for Poisson Surface Reconstruction and smooth the generated muscle mesh
- attach_muscle
    - Given the whole muscle mesh and selected bone faces, perform single-face attachment from muscle to each desired attachment sites on bone using ARAP deformation
- attach_muscle_multiface
    - Given the whole muscle mesh and selected bone faces, perform multi-face attachment from muscle to each desired attachment sites on bone using ARAP deformation
- attach_tendon
    - Given the whole muscle mesh and selected bone faces, generate tendon mesh connecting the muscle and the desired attachment sites on bone.
    
#### poisson_surface_reconstruction.h
- poisson_surface_reconstruction
    - Takes input sample points P and input normals N and gives a watertight mesh using a simplified version of [Kazhdan et. al 2006]

#### PCA_elipse.h
- PCA_param
    - Given a patch mesh, fit an ellipse using PCA. Output the parameters defining an ellipse including short/long axis length, short/long axis direction, center of the ellipse and normal of the ellipse plane.

#### volume_along_curve.h
- volume_along_curve
    - Given points along a curve, genereate point clouds and normals from circles around the curve for the purpose of Poisson Surface Reconstruction.
- ellipse_along_curve
    - Given points and other interpolated ellipse parameters along a curve, genereate point clouds and normals from ellipses around the curve for the purpose of Poisson Surface Reconstruction.

#### triangle_hunt.h
- triangle_hunts
    - Given a surface mesh and a query point, find out the closest face on the mesh from the query point.

- triangle_hunt_lst
    - Given a surface mesh and a query point, return a sorted list of the distances of faces on the mesh from the query point.

### User Input

Much of the user interface functionality was directly coded into the main.cpp file, however some repeatedly used or more complex chunks of code were wrapped into their own functions and files.

#### ui_helpers.h
A variety of functions that are used for user-input.

- intersection_with_xy_plane
    - User input currently defaults to the xy plane for simplification in 3D-- this function returns an intersection point in 3D using raycasting.
- add_face_to_patch
    - Used to select faces on the bone mesh, and group them into "patches" if they share an edge (at the time the face is selected).
- verts_within_x_range
    - Used for the experimental inflation/deflation mode, highlights vertices within the x-range of some point (the cursor)

#### mesh_editing.h
- xflate_mesh
    - Either inflate or deflate a whole given mesh, by moving all vertices in the direction of the normal.
- xflate_verts
    - Either inflate or deflate only a subset of vertices on a mesh.
- xflate_verts_in_xrange
    - Either inflate or deflate verts within the x-range of a given point.
- smooth_mesh_with_fixed
    - Smooth mesh but pin fixed verts back into place (work-around of non-volume-preserving smoothing)

#### generate_bone.h
- start_flinstone_bone
    - Given a point, generate new "flinstone bone" mesh with one end pinned to that point.
- transform_flinstone_bone
    - Given a start point and an end point, transform a given "flinstone bone" mesh lengthwise along the vector between the points
