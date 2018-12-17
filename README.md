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

Once built, you can execute the assignment from inside the `build/` using 

    ./introduction
