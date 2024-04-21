# <span style="font-variant:small-caps;">TriMe++</span>

<span style="font-variant:small-caps;">TriMe++</span>: a **multi-threaded** software library for **fast** generation of **large-scale**, **adaptive** 2D triangle meshes for **intricate 
geometric shapes** using Delaunay triangulation. 

By Jiayin Lu [^a][^b], Chris H. Rycroft [^a]

![header plot](/docs/github_header_plot.jpg)

[^a]: University of Wisconsin-Madison

[^b]: Harvard University

Please consider [citing our paper](#citation-and-paper).

I am very interested to hear from users of the software, so if you find this
useful or have any questions, please email me at jiayinlu19960224work@gmail.com.


Table of contents
================================================
1. [Overview](#overview)

      1.1. [Directory structure](#directory-structure)

      1.2. [Method](#method)

      1.3. [Related programs](#related-programs)

      1.4. [Citation and paper](#citation-and-paper)

1. [Quick start](#quick-start)

      3.1. [Compilation](#compilation---linux--mac-os--windows-with-cygwin)
      
      3.2. [Run the examples](#running-the-examples)

      3.3. [Understand the output files](#understanding-the-output-files)

      3.4. [Customization](#customization)

1. [Performance](#performance)

1. [Code updates](#updates)

1. [Acknowledgement](#acknowledgement)

1. [References](#references)


Overview
========================================================================
<span style="font-variant:small-caps;">TriMe++</span> is especially suited for the fast generation of **large-scale**, **adaptive** and **high-quality** triangle mesh for complicated shapes in 2D. 

<span style="font-variant:small-caps;">TriMe++</span> implements three iterative meshing algorithms: 
DistMesh [[1]](#references), Centroidal Voronoi Diagram [[2]](#references) meshing, and a hybrid method of the two. 

In each meshing iteration, the multi-threaded version of <span style="font-variant:small-caps;">Voro++</span> [[3]](#references) is used to obtain fast and efficient generation of the Delaunay triangulation. 

Multi-threaded parallel 
computing is implemented throughout the meshing procedure via OpenMP [[4]](#references). 


Directory structure
-------

<code>TriMe</code>: Root-directory
><code>src_voro</code>: Sub-directory with source code files of the multi-threaded <span style="font-variant:small-caps;">Voro++</span> program
>
><code>src_TriMe</code>: Sub-directory with source code files of the <span style="font-variant:small-caps;">TriMe++</span> program
>
><code>Examples</code>: Sub-directory with example files of meshing on geometry shapes making use of the library
>
><code>config.mk</code>: File specifying <span style="font-variant:small-caps;">C++</span> compilation and installation settings


Method
------------
<span style="font-variant:small-caps;">TriMe++</span> is built around several <span style="font-variant:small-caps;">C++</span> classes that follows a clearly structured meshing pipeline. Throughout the meshing workflow, we utilize OpenMP for multi-threaded parallel computations.

During initialization, the <code>shape_2d</code> class reads the geometry input and generates a signed distance field 
using a grid-based data structure to represent the shape. The <code>sizing_2d</code> class subsequently produces 
adaptive element sizing and density fields for the mesh. It utilizes an adaptive quad-tree data structure, 
enabling efficient refinement of sizing and density values in areas with complex geometries. 

In the meshing procedure, the <code>parallel_meshing_2d</code> class iteratively improves point positions in the mesh. 
In each meshing iteration, we employ the multi-threaded <span style="font-variant:small-caps;">Voro++</span> to facilitate fast and efficient parallel computation 
in generating the Delaunay triangulation of the points. Users have the option to select from three meshing 
algorithms, the DistMesh algorithm in the <code>mesh_alg_2d_dm</code> class, the Centroidal Voronoi Diagram meshing 
algorithm in the <code>mesh_alg_2d_cvd</code> class, and a hybrid method of the two in the <code>mesh_alg_2d_hybrid</code> class. 

Related programs
---------------
We use the multi-threaded version of <span style="font-variant:small-caps;">Voro++</span> [[3]](#references) in the meshing algorithms. 

Citation and paper
------
<span style="font-variant:small-caps;">TriMe++</span> is released as a free software.

If you plan to publish an academic paper using this software, please consider citing the following
publications:

> Lu, Jiayin, Rycroft, Chris. (2023). <span style="font-variant:small-caps;">TriMe++</span>: Multi-threaded triangular meshing in two dimensions
>[[ArXiv]](https://arxiv.org/abs/2309.13824)
>
> **Abstract**
>
>Geometry meshing gives discrete representation of continuous geometry shapes. It has important applications in computer graphics and in scientific computing, especially for numerical simlations using the finite element method. In some large-scale simulations, high-resolution geometry meshes using hundreds of thousands of points may be needed. In this paper, we present <span style="font-variant:small-caps;">TriMe++</span>, a multi-threaded software library designed for generating adaptive 2D meshes for intricate geometric shapes using Delaunay triangulation. Multi-threaded parallel computing is implemented throughout the meshing procedure, making it especially suitable for fast generation of large-scale meshes. Three iterative meshing algorithms are provided to users: DistMesh, Centroidal Voronoi Diagram meshing, and a hybrid method of the two. We test the performance of <span style="font-variant:small-caps;">TriMe++</span>. We compare the three meshing methods, and show that the hybrid method has the advantages of the other two, while avoiding the disadvantages. We demonstrate that the software library achieves significant speed up in time with parallel computing when generating a large-scale mesh of $10^6$ points. We also show that <span style="font-variant:small-caps;">TriMe++</span> is able to handle extremely complicated geometry shapes and generates adaptive meshes of high quality.

Further improvements of <span style="font-variant:small-caps;">TriMe++</span> will be made and documented in the [Code updates](#code-updates) section. 



Quick start
==================================================


Compilation - Linux / Mac OS / Windows with Cygwin
---------------
The code is written in ANSI <span style="font-variant:small-caps;">C++</span>, and compiles on many system architectures. The
package contains the <span style="font-variant:small-caps;">C++</span> source code and example files. On Linux, Mac OS, and 
Windows (using Cygwin), the compilation and installed can be carried out using GNU Make.

To begin, the user should review the file "config.mk" in the top level
directory, to make sure that the compilation and installation settings are
appropriate for their system. 

Typing "make" in the Example directory will then compile the static
library and examples.

Run the examples
---------------


Understand the output files
---------------


Customization
---------------


Performance
================================================


Code updates
================================================
Oct 19th, 2023: 

1. Added in fixed point input for meshing with fixed points.

2. Implemented half-edge structure for final mesh clean up and outputting boundary vertices in
    counter-clockwise (CCW) order.

3. Triangle vertex ID's are outputted in counter-clockwise (CCW) order.


Acknowledgement
=====
This research was supported by a grant from the United States--Israel Binational
Science Foundation (BSF), Jerusalem, Israel through grant number 2018/170.
C.~H.~Rycroft was partially supported by the Applied Mathematics
Program of the U.S. DOE Office of Advanced Scientific Computing Research under
contract number DE-AC02-05CH11231.


References
=============================================
[1] DistMesh 
 - P.-O. Persson, G. Strang, A simple mesh generator in matlab, SIAM Review 46 (2) (June 2004) 329–345. [[website]](http://persson.berkeley.edu/distmesh/)

 - P.-O. Persson, Mesh generation for implicit geometries, Ph.D. thesis, Massachusetts Institute of 
   Technology (2005).

[2] Centroidal Voronoi Diagram meshing
 - Q. Du, V. Faber, M. Gunzburger, Centroidal Voronoi tessellations: Applications and algorithms, 
   SIAM Review 41 (4) (1999) 637–676.

 - Q. Du, D. Wang, Tetrahedral mesh generation and optimization based on centroidal voronoi tessellations, 
   International Journal on Numerical Methods in Engineering 56 (9) (2003) 1355–1373.

[3] Multi-threaded <span style="font-variant:small-caps;">Voro++</span>
- Chris H. Rycroft, "<span style="font-variant:small-caps;">Voro++</span>: A three-dimensional Voronoi cell library in C++",
  Chaos 19, 041111 (2009). [[website]](https://math.lbl.gov/voro++/)

- Lu, Jiayin, Lazar, Emanuel, and Rycroft, Chris. (2023). "An extension to <span style="font-variant:small-caps;">Voro++</span> for 
  multithreaded computation of Voronoi cells". Computer Physics Communications. 291. 108832 (2023). [[website]](https://github.com/chr1shr/voro/tree/dev)

[4] OpenMP
 - L. Dagum, R. Menon, OpenMP: an industry standard API for shared-memory pro- gramming, IEEE Computational 
   Science and Engineering 5 (1) (1998) 46–55. doi: 10.1109/99.660313. [[website]](https://www.openmp.org)


