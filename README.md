# Geometry Processing – Muscle Mesh Project

> **To get started:** Fork this repository then issue
> 
>     git clone --recursive http://github.com/[username]/geometry-processing-introduction.git
>

## Prerequisite installation

On all platforms, we will assume you have installed cmake and a modern c++
compiler on Mac OS X[¹](#¹macusers), Linux[²](#²linuxusers), or
Windows[³](#³windowsusers).

We also assume that you have cloned this repository using the `--recursive`
flag (if not then issue `git submodule update --init --recursive`). 

## Layout

    README.md
    CMakeLists.txt
    main.cpp
    include/
      function1.h
      function2.h
      ...
    src/
      function1.cpp
      function2.cpp
      ...
    shared/
      libigl/
        include/
          igl/
            ...
      ...

## Compilation

Starting in this directory, issue:

    mkdir build
    cd build
    cmake ..
    make 

## Execution

Once built, you can execute the program from inside the `build/` using 

    ./introduction
    
## Documentation

An overview of each header file, documenting each function.

### Muscle and Tendon Generation

- bezier.h
- poisson_surface_reconstruction.h
- (all the other code

### User Input

Much of the user interface functionality was directly coded into the main.cpp file, however some repeatedly used or more complex chunks of code were wrapped into their own functions and files.

#### ui_helpers.h
A variety of functions that are used for user-input.

- intersection_with_xy_plane()
    - User input currently defaults to the xy plane for simplification in 3D-- this function returns an intersection point in 3D using raycasting.
- add_face_to_patch()
    - Used to select faces on the bone mesh, and group them into "patches" if they share an edge.
- verts_within_x_range()
    - Used for the experimental inflation/deflation mode, highlights vertices within the x-range of some point (the cursor)

#### mesh_editing.h

#### generate_bone.h
