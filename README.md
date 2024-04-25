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

      3.1. [Compilation](#compilation): Linux / Mac OS / Windows with Cygwin
      
      3.2. [Run a basic example](#run-a-basic-example): built-in shape primitives, mesh adaptivity control, parallel threads control

      3.3. [Understand the output files](#understand-the-output-files): vertices, triangles, edges, boundaries,  mesh quality measures

      3.4. [Customization](#customization): user-defined sizing field, shape boolean operations, custom shape from user-defined signed distance function, custom shape from contour line segments, custom shape from binary image

1. [Performance](#performance)

1. [Code updates](#code-updates)

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
> **Code**
>
><code>config.mk</code>: File specifying <span style="font-variant:small-caps;">C++</span> compilation and installation settings
>
><code>src_voro</code>: Sub-directory with source code files of the multi-threaded <span style="font-variant:small-caps;">Voro++</span> program
>
><code>src_TriMe</code>: Sub-directory with source code files of the <span style="font-variant:small-caps;">TriMe++</span> program
>
><code>Examples</code>: Sub-directory with example files of meshing on geometry shapes making use of the library
>
> **Documentation**
>
><code>README.md</code>: Documentation and quick start guide
>
><code>docs</code>: Sub-directory with figures for the README.md file


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


Compilation
---------------
The code is written in ANSI <span style="font-variant:small-caps;">C++</span>, and compiles on many system architectures. The
package contains the <span style="font-variant:small-caps;">C++</span> source code and example files. On Linux, Mac OS, and 
Windows (using Cygwin), the compilation and installed can be carried out using GNU Make.

To begin, the user should review the file "config.mk" in the top level
directory, to make sure that the compilation and installation settings are
appropriate for their system. 

Then type <code>make</code> in the command line in the Example directory will compile the static
library and examples. And we now see executable example files in the directory. 

Run a basic example
---------------

In the <code>/Examples</code> directory, let's first look at the most basic example file, <code>primitive_shape_meshing.cc</code>. This file does meshing on built-in primitive shapes provided by the library: rectangle, circles, triangles. 

To run its meshing, we simply need to run the corresponding executable file by typing <code>./primitive_shape_meshing</code> in the command line. 

Let's now look at the code in detail in <code>primitive_shape_meshing.cc</code>. 

### Set up
On the top, use <code>#include "parallel_meshing_2d.hh"</code> to import the <span style="font-variant:small-caps;">TriMe++</span> code. 

In the main code, these few lines at the beginning of the code specify:
- <code>num_t</code>: The number of parallel threads to use.
- <code>method_ind</code>: The meshing algorithm to use. $0$ DistMesh; $1$ CVD; $2$ Hybrid.
- <code>Ntotal</code>: The total number of meshing points.
- <code>K</code>: Adaptivity of mesh. $0$ Uniform; $>1$ Adaptive; Larger <code>K</code> increases the adaptivity of mesh.
- <code>output_interval</code>: File output frequency. $0$ No output; $-1$ Only output at termination; any positive integer, e.g. $10$, represents outputting every $10$ triangulations during the meshig procedure. 

Here, we are using $4$ parallel threads, the hybrid meshing method, $5000$ meshing vertices, with an adaptivity $K=0.1$, and output the meshing data at termination. 
```c++
//Specify number of parallel threads
int num_t=4;
//Method index, 0 Distmesh; 1 CVD; 2 Hybrid
int method_ind=2; 
//Number of vertices in the mesh
int Ntotal=5000; 
//0:uniform sizing field
double K=0.1; 
//File output frequency. 0 no output; -1 output at termination; 10, every 10 triangulations output
int output_interval=-1; 
```

The next few lines create a directory <code>/Examples/case_name_base</code> to store the output files. 
```c++
//File output name prefix
char case_name_base[256];
sprintf(case_name_base,"triangle_mesh_N_%d_K_%g",Ntotal,K);
//Create directory to store file outputs
mkdir(case_name_base,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
```
### Multi-threaded <span style="font-variant:small-caps;">Voro++</span>
```c++
//------------------1.Create container----------------------
int cnx=sqrt(Ntotal/3.3); int cny=cnx;
container_2d con(0.0,1.0,0.0,1.0,cnx,cny,false,false,16,num_t);
```

### Built-in 2D primitive shapes

A few built-in primitive shapes are provided in the library. They are rectangle, circle and triangle: 
><code>shape_2d_rectangle(container_2d &con, int num_t, double x0, double x1, double y0, double y1)</code>: A rectangular shape defined on the domain $[x0,x1]\times[y0,y1]$.
>
><code>shape_2d_circle(container_2d &con, int num_t, double r, double x0, double y0)</code>: A circular shape with radius $r$ centered at $(x0,y0)$. 
>
><code>shape_2d_triangle(container_2d &con, int num_t, double x0, double y0, double x1, double y1, double x2, double y2)</code>: A triangle defined by three vertices, $(x0,y0)$, $(x1,y1)$ and $(x2,y2)$.

In the syntax of the primitive shapes creation above, we need to pass in the <span style="font-variant:small-caps;">Voro++</span> <code>container_2d</code> object created, and pass in the number of parallel threads to use for calculating the Voronoi diagram in the Delaunay triangulation. 

In the following code, we create a triangle shape defined by the three vertices $(0.1,0.1)$,$(0.3,0.8)$ and $(0.8,0.5)$. The commented code below also provides examples for creating a circle and a rectangle primitive. 
```c++
//-------------------2.Create shape------------------------
//A circle with radius 0.3 centered at (0.5,0.5)
//shape_2d_circle shp(con,num_t,0.3,0.5,0.5);

//A rectangle on domain [0.1,0.9]x[0.3,0.7]
//shape_2d_rectangle shp(con,num_t,0.1, 0.9, 0.3, 0.7);

//A triangle defined by the three vertices (0.1,0.1),(0.3,0.8) and (0.8,0.5)
shape_2d_triangle shp(con,num_t,0.1,0.1,0.3,0.8,0.8,0.5);
```
The plot here shows three example meshes using the different primitive shapes and adaptivity paramemters. 
![primitives meshes plot](/docs/primitives_meshes_plot.jpg)


### Meshing
Next, we create the adaptive sizing field for meshing, based on the shape defined and the adaptivity parameter <code>K</code>. 
```c++
//-------------------3.Create adaptive sizing field------------------
sizing_2d_automatic size_field(&shp,K);
```

We then create the parallel_meshing_2d object for the meshing procedure later. It takes in as parameters: 
- the <span style="font-variant:small-caps;">Voro++</span> <code>container_2d</code> object for Delaunay triangulation, 
- the <code>shape_2d</code> object that defines the geometry to mesh, 
- the adaptive sizing field <code>size_field</code> that defines the relative triangle element sizes in the mesh, 
- the number of parallel threads to use in the meshing algorithm <code>num_t</code>, 
- the frequency of file outputs <code>output_interval</code>,
- the output file names prefix <code>case_name_base</code>.
```c++
//------------------4.Create parallel_meshing_2d object----------------------
//Create the pm2d object
parallel_meshing_2d pm2d(&con, &shp, &size_field, num_t, output_interval, case_name_base);
```

We then initialize meshing points. 
```c++
//Optional: set seed for rand() for generating randome points in pt_init
srand(10);
//Initialize meshing points
pm2d.pt_init(Ntotal);
```

Afterwards, we start the iterative meshing procedure, depending on the meshing algorithm chosen. 
For a chosen meshing algorithm, the code first creates an object <code>mesh_method</code> that implements the meshing algorithm. It reads in the <code>pm2d</code> object to obtain all the information of meshing specifications, such as the shape, adaptivity and number of threads to use. 
The <code>pm2d</code> then calls <code>pm2d.meshing(&mesh_method)</code> to do the iterative meshing using that chosen method. 
```c++
//======================Parallel meshing=====================
if(method_ind==0){ //DistMesh
      mesh_alg_2d_dm mesh_method(&pm2d);
      pm2d.meshing(&mesh_method); 
}
else if(method_ind==1){ //CVD
      mesh_alg_2d_cvd mesh_method(&pm2d);
      pm2d.meshing(&mesh_method); 
}
else if(method_ind==2){ //Hybrid
      mesh_alg_2d_hybrid mesh_method(&pm2d);
      pm2d.meshing(&mesh_method); 
}
```


Understand the output files
---------------
Let's now look at the output files from running <code>./primitive_shape_meshing</code>. They are in the directory <code>/Example/triangle_mesh_N_5000_K_0.1</code>. 

### Triangulation data

In <span style="font-variant:small-caps;">TriMe++</span>, each point is associated with an ID <code>pid</code> and x and y coordinates <code>(px, py)</code>. 
> point: <code>pid, (px, py)</code>

A triangle is associated with an ID <code>tid</code>, the three point IDs of its three vertices, <code>(vid0, vid1, vid2)</code>, and two quality measures, aspect ratio <code>$\alpha$</code> and edge ratio <code>$\beta$</code>.

> triangle: <code>tid, (vid0, vid1, vid2), ($\alpha$, $\beta$)</code>

The aspect ratio <code>$\alpha$</code> and edge ratio <code>$\beta$</code> are calculated by 

$$\alpha=\frac{R_{\text{circum}}}{2\cdot R_{\text{in}}}, \quad \beta=\frac{l_{\max}}{l_{\min}}$$

where $R_{\text{circum}}$ and $R_{\text{in}}$ are the circumradius and the inradius of the triangle, respectively. $l_{\max}$ and $l_{\min}$ are the longest and shortest edge lengths of the triangle, respectively. Both $\alpha$ and $\beta$ have minimum values of $1$, corresponding to an equilateral triangle. $\alpha$ has large values for badly shaped triangles, including both needle-type and flat-type triangles shown in plot (a) and (b) in the figure below. $\beta$ has large values for needle-type triangles.

<p align="center">
<img src="/docs/needle_flat_triangles.png" width="300" />
</p>

### Output files

The filename prefix is <code>fp="triangle_mesh_N_5000_K_0.1"</code>. The next description phrase <code>ti</code> describe the meshing iteration outputed: <code>fp_ti_=fp_1_</code> is the initial triangulation; <code>fp_10_</code> is the triangulation of the $10^{th}$ iteration. If <code>output_interval=-1</code> is set to only output at termination, then the description phrase is <code>fp_final_</code>. 

The rest of the description phrases in the filenames desribes the data being outputted. Suppose we have $N$ points and $M$ triangles generated. 

**Points**

> <code>fp_ti_xy_id.txt</code>: Each row format is <code>[x y]</code>. The particle ID and coordinates. The particle ID <code>[0,1,...,N-1]</code> is implicity implied by the line number. The first line corresponds to the coordinates of point $0$, and the tenth line corresponds to the coordinates of point $9$.

**Edges**

> <code>fp_ti_tria_bar_ids.txt</code>>: Each row format is <code>[vid0 vid1]</code>, the point IDs corresponding to the end points of each edge in the triangulation. The triangulation edges are unique and non-overlapping in the output. 
> 
> <code>fp_ti_tria_bar_coords.txt</code>: The coordinates of the end points for each triangulation edge. Each edge consists of two rows, each row is the coordiate of one of the end poins: 
>> <code>[x0 y0]</code>
>>
>> <code>[x1 y1]</code>
>
> An empty line then follows, separating the current edge output and the next edge output. 

**Triangles**

> <code>fp_ti_tria_vertex_ids.txt</code>: Each row format is <code>[tid vid0 vid1 vid2]</code>, the triangle ID followed by the three particle IDs of its three vertices. 
>
> <code>fp_ti_tria_vertex_coords.txt</code>: Each row format is <code>[tid x0 y0 x1 y1 x2 y2]</code>, the triangle ID followed by the three particle coordinates of its three vertices. 

**Mesh quality**

> <code>fp_ti_tria_quality_stat_ar.txt</code>: This file consists of a single line, listing the aspect ratios $\alpha_i$ of each triangle $i$, separated by space. The format is: <code>[$\alpha_0$ $\alpha_1$ $\alpha_2$ ... $\alpha_M$]</code>.
>
> <code>fp_ti_tria_quality_stat_er.txt</code>: This file consists of a single line, listing the edge ratios $\beta_i$ of each triangle $i$, separated by space. The format is: <code>[$β0$ $\beta_1$ $\beta_2$ ... $\beta_M$]</code>.
>
> <code>fp_tria_quality_stat_overall.txt</code>: This file outputs the overall mesh quality at the corresponding outputted triangulation iterations. Each row consists of $12$ numbers, they are (in order): 
> - (0) the triangulation iteration number,
> - (1) the number of triangles,
> - (2) the maximum $\alpha$, (3) the maximum $\beta$,
> - (4) the median $\alpha$, (5) the median $\beta$,
> - (6) the mean $\alpha$, (7) the mean $\beta$,
> - (8) the $\frac{1}{2}$-mean of $\alpha$, (9) the $\frac{1}{2}$-mean of $\beta$, 
> - (10) the standard deviation of $\alpha$, (11) the standard deviation of $\beta$. 

The $\frac{1}{2}$-mean mentioned above is calculated as, 
$$M_{\frac{1}{2}}(x_1,\ldots,x_n)=\left( \frac{1}{n}\sum_{i=1}^n x_i^{\frac{1}{2}} \right)^2,$$
which is less sensitive to large outliers than the arithmetic mean, and thus more suitable as an indicator for overall mesh quality.

**Shape boundaries**

Suppose the shape have $B$ boundaries, $b_1, b_2, ...b_B$. We use an integer <code>bi</code> to denote the $b_i^{\text{th}}$ boundary. The boundary numbering $b_i$ is based on descending order of the number of boundary points on each boundary. 

>fp_ti_bdry_vertices_ids_CCW_bi.txt</code>: The boundary point IDs in counter-clockwise order of the $b_i^{\text{th}}$ boundary. Each row is a point ID. Suppose the boundary has $n$ points, then there are $n$ rows. That is, the beginning point is not repeated at the end. 
>
>fp_ti_bdry_vertices_coords_CCW_bi.txt</code>: The boundary point coordinates in counter-clockwise order of the $b_i^{\text{th}}$ boundary. The points have the same ordering and corresponds to the point IDs in <code>fp_ti_bdry_vertices_ids_CCW_bi.txt</code>. Each row format is <code>[x y]</code>, the coordinates of the boundary point. Again, the number of rows is $n$, as the beginning point is not repeated at the end. 

Customization
---------------
### user-defined sizing field

### shape boolean operations

### custom shape from user-defined signed distance function

### custom shape from contour line segments

### custom shape from binary image





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


