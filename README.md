# TriMe++

TriMe++: a multi-threaded software library for 2D geometry meshing using the Delaunay triangulation. 

By Jiayin Lu (1,2), Chris H. Rycroft (2)

(1) Harvard University

(2) University of Wisconsin-Madison

Overview
========================================================================
TriMe++ is a multi-threaded software library designed for generating adaptive 2D meshes for intricate 
geometric shapes using Delaunay triangulation. TriMe++ implements three iterative meshing algorithms: 
DistMesh (1), Centroidal Voronoi Diagram (2) meshing, and a hybrid method of the two. Multi-threaded parallel 
computing is used throughout the meshing procedure via OpenMP (3). 

(1)
 - P.-O. Persson, G. Strang, A simple mesh generator in matlab, SIAM Review 46 (2) (June 2004) 329–345, 
   http://persson.berkeley.edu/distmesh/.

 - P.-O. Persson, Mesh generation for implicit geometries, Ph.D. thesis, Massachusetts Institute of 
   Technology (2005).

(2)
 - Q. Du, V. Faber, M. Gunzburger, Centroidal Voronoi tessellations: Applications and algorithms, 
   SIAM Review 41 (4) (1999) 637–676.

 - Q. Du, D. Wang, Tetrahedral mesh generation and optimization based on centroidal voronoi tessellations, 
   International Journal on Numerical Methods in Engineering 56 (9) (2003) 1355–1373.

(3)
 - L. Dagum, R. Menon, OpenMP: an industry standard API for shared-memory pro- gramming, IEEE Computational 
   Science and Engineering 5 (1) (1998) 46–55. doi: 10.1109/99.660313. https://www.openmp.org.


Method
=============================================================
TriMe++ is built around several C++ classes that follows a clearly structured meshing pipeline. 
During initialization, the shape_2d class reads the geometry input and generates a signed distance field 
using a grid-based data structure to represent the shape. The sizing_2d class subsequently produces 
adaptive element sizing and density fields for the mesh. It utilizes an adaptive quad-tree data structure, 
enabling efficient refinement of sizing and density values in areas with complex geometries. 
In the meshing procedure, the parallel_meshing_2d class iteratively improves point positions in the mesh. 
In each meshing iteration, we employ the multi-threaded Voro++ to facilitate parallel computation 
in generating the Delaunay triangulation of the points. Users have the option to select from three meshing 
algorithms, the DistMesh algorithm in the mesh_alg_2d_dm class, the Centroidal Voronoi Diagram meshing 
algorithm in the mesh_alg_2d_cvd class, and a hybrid method of the two in the mesh_alg_2d_hybrid class. 
Throughout this meshing workflow, we utilize OpenMP for multi-threaded parallel computations.




Compilation - Linux / Mac OS / Windows with Cygwin
==================================================
The code is written in ANSI C++, and compiles on many system architectures. The
package contains the C++ source code and example files. On Linux, Mac OS, and 
Windows (using Cygwin), the compilation and installed can be carried out using GNU Make.

To begin, the user should review the file "config.mk" in the top level
directory, to make sure that the compilation and installation settings are
appropriate for their system. Typing "make" will then compile the static
library and examples.


Related programs
================
We use the multi-threaded version of Voro++ in the meshing algorithms. 

- Chris H. Rycroft, "Voro++: A three-dimensional Voronoi cell library in C++",
  Chaos 19, 041111 (2009).

- Lu, Jiayin, Lazar, Emanuel, and Rycroft, Chris. (2023). "An extension to Voro++ for 
  multithreaded computation of Voronoi cells". Computer Physics Communications. 291. 108832 (2023). 


Contents
========
examples - examples making use of the library

src - source code files


Usage
=====
TriMe++ is released as a free software.

I am very interested to hear from users of the software, so if you find this
useful, please email me at jiayinlu19960224work@gmail.com. Also, if you plan to publish an
academic paper using this software, please consider citing the following
publications:

- Lu, Jiayin, Rycroft, Chris. (2023). TriMe++: Multi-threaded geometry meshing using the Delaunay Triangulation

We detail in the reference paper the methods, designs and structure of the meshing pipeline
and algorithms. We also describe the implementation of multi-threaded parallelizaiton throughout 
the meshing procedure. An example code is given in the paper. We also test the performance of the 
software in the paper. We compare the three meshing methods, and show that the hybrid method has 
the advantages of the other two, while avoiding the disadvantages. We demonstrate with examples that 
the meshing achieves significant speed up in time using parallel computing. We also show that TriMe++ 
is able to handle extremely complicated geometry shapes and generates adaptive meshes of high quality.

Further improvements of TriMe++ will be made in the future. 



Acknowledgement
=====
This research was supported by a grant from the United States--Israel Binational
Science Foundation (BSF), Jerusalem, Israel through grant number 2018/170.
C.~H.~Rycroft was partially supported by the Applied Mathematics
Program of the U.S. DOE Office of Advanced Scientific Computing Research under
contract number DE-AC02-05CH11231.