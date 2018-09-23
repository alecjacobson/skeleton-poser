# skeleton poser

> **Get started with:**
>
>     git clone --recursive http://github.com/alecjacobson/skeleton-poser

This small application loads in a 3D model as a mesh and a [linear blend
skinning](http://skinning.org) skeleton as a [.tgf
file](http://libigl.github.io/libigl/file-formats/tgf/). The application will
attempt to compute [bounded biharmonic
weights](http://libigl.github.io/libigl/tutorial/#bounded-biharmonic-weights),
then an interactive OpenGL window opens where the user can select and manipulate
bones with a rotation widget.

## Compilation

    mkdir build
    cd build
    cmake ../
    make

## Execution 

    ./skeleton_poser mesh.obj skeleton.tgf [weights.dmat]
