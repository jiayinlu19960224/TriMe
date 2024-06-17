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

      3.4. [Customization](#customization): 
            
      - [Changing number of parallel threads](#changing-number-of-parallel-threads)
      - [User-defined sizing field](#user-defined-sizing-field)
      - [Shape boolean operations](#shape-boolean-operations)
      - [Custom shape from user-defined signed distance function](#custom-shape-from-user-defined-signed-distance-function-sdf)
      - [Custom shape from contour line segments](#custom-shape-from-contour-line-segments)
      - [Setting fixed points in the mesh](#setting-fixed-points-in-the-mesh)

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
  >Best efficient performance is reached when the number of threads $\leq$ the number of physical cores in the computer.
  
  >When the number of threads is more than that, we can usually still get improvement in computation time, just less efficient. 
  
  >It is recommended that the number of threads is never more than twice the number of physical cores. 
- <code>method_ind</code>: The meshing algorithm to use. $0$ DistMesh; $1$ CVD; $2$ Hybrid.
- <code>Ntotal</code>: The total number of meshing points specified. At meshing termination, we employed a procedure to [delete problematic points from the mesh](#may-19th-2024-delete-problematic-vertices-at-termination). Therefore, the final mesh number of points is $\leq$ ``Ntotal``. 
- <code>K</code>: Adaptivity of mesh. $0$ Uniform; $>1$ Adaptive; Larger <code>K</code> increases the adaptivity of mesh.
- <code>output_interval</code>: File output frequency. 
  >$0$ No output; 
  >
  >$-1$ Only output at termination; 
  >
  >$-2$ Output initial and final triangulation data
  >
  >Any positive integer, e.g. $10$, represents outputting every $10$ triangulations during the meshig procedure. 

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
Next, we generate automatically the adaptive sizing field for meshing, based on the shape defined and the adaptivity parameter <code>K</code>. 
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

A triangle is associated with an ID <code>tid</code>, the three point IDs of its three vertices, <code>(vid0, vid1, vid2)</code>, and two quality measures, aspect ratio <code>α</code> and edge ratio <code>β</code>.

> triangle: <code>tid, (vid0, vid1, vid2), (α, β)</code>

The aspect ratio $\alpha$ and edge ratio $\beta$ are calculated by 

$$\alpha=\frac{R_{\text{circum}}}{2\cdot R_{\text{in}}}, \quad \beta=\frac{l_{\max}}{l_{\min}}$$

where $R_{\text{circum}}$ and $R_{\text{in}}$ are the circumradius and the inradius of the triangle, respectively. $l_{\max}$ and $l_{\min}$ are the longest and shortest edge lengths of the triangle, respectively. Both $\alpha$ and $\beta$ have minimum values of $1$, corresponding to an equilateral triangle. $\alpha$ has large values for badly shaped triangles, including both needle-type and flat-type triangles shown in plot (a) and (b) in the figure below. $\beta$ has large values for needle-type triangles.

<p align="center">
<img src="/docs/needle_flat_triangles.png" width="300" />
</p>

### Output files

The filename prefix is <code>fp="triangle_mesh_N_5000_K_0.1"</code>. The next description phrase <code>ti</code> describe the meshing iteration outputted: <code>fp_ti_=fp_1_</code> is the initial triangulation; <code>fp_10_</code> is the triangulation of the $10^{th}$ iteration. For the last triangulation at termination, the phrase is set to <code>fp_final_</code>.

The rest of the description phrases in the filenames desribes the data being outputted. Suppose we have $N$ points and $M$ triangles generated. 

**Points**

> <code>fp_ti_xy_id.txt</code>: Each row format is <code>[x y]</code>. The particle ID and coordinates. The particle ID <code>[0,1,...,N-1]</code> is implicity implied by the line number. The first line corresponds to the coordinates of point $0$, and the tenth line corresponds to the coordinates of point $9$.

**Edges**

> <code>fp_ti_tria_bar_ids.txt</code>: Each row format is <code>[vid0 vid1]</code>, the point IDs corresponding to the end points of each edge in the triangulation. The triangulation edges are unique and non-overlapping in the output. 
> 
> <code>fp_ti_tria_bar_coords.txt</code>: The coordinates of the end points for each triangulation edge. Each edge consists of two rows, each row is the coordiate of one of the end poins: 
>> <code>[x0 y0]</code>
>>
>> <code>[x1 y1]</code>
>
> An empty line then follows, separating the current edge output and the next edge output. 

**Triangles**

> <code>fp_ti_tria_vertex_ids.txt</code>: Each row format is <code>[tid vid0 vid1 vid2]</code>, the triangle ID followed by the three particle IDs of its three vertices. At termination, the triangle vertices IDs are outputted in counter-clockwise (CCW) order. 
>
> <code>fp_ti_tria_vertex_coords.txt</code>: Each row format is <code>[tid x0 y0 x1 y1 x2 y2]</code>, the triangle ID followed by the three particle coordinates of its three vertices. At termination, the triangle vertices coordinates are outputted in CCW order. 

**Mesh quality**

> <code>fp_ti_tria_quality_stat_ar.txt</code>: This file consists of a single line, listing the aspect ratios $\alpha_i$ of each triangle $i$, separated by space. The format is: <code>[α0 α1 α2 ... αM]</code>.
>
> <code>fp_ti_tria_quality_stat_er.txt</code>: This file consists of a single line, listing the edge ratios $\beta_i$ of each triangle $i$, separated by space. The format is: <code>[β0 β1 β2 ... βM]</code>.
>
> <code>fp_tria_quality_stat_overall.txt</code>: This file outputs the overall mesh quality at the corresponding outputted triangulation iterations. Each row consists of $12$ numbers, they are (in order): 
> - (0) the triangulation iteration number,
> - (1) the number of vertices, 
> - (2) the number of triangles,
> - (3) the maximum $\alpha$, (4) the maximum $\beta$,
> - (5) the median $\alpha$, (6) the median $\beta$,
> - (7) the mean $\alpha$, (8) the mean $\beta$,
> - (9) the $\frac{1}{2}$-mean of $\alpha$, (10) the $\frac{1}{2}$-mean of $\beta$, 
> - (11) the standard deviation of $\alpha$, (12) the standard deviation of $\beta$. 

The $\frac{1}{2}$-mean mentioned above is calculated as, 
$$M_{\frac{1}{2}}(x_1,\ldots,x_n)=\left( \frac{1}{n}\sum_{i=1}^n x_i^{\frac{1}{2}} \right)^2,$$
which is less sensitive to large outliers than the arithmetic mean, and thus more suitable as an indicator for overall mesh quality.

**Shape boundaries**

We can output the vertices and edges made up of the shape boundaries in counter-clockwise order at termination. Note that this is only outputted at termination.

Suppose at termination, the shape have $B$ boundaries, $b_1, b_2, ...b_B$. We use an integer <code>bi</code> to denote the $b_i^{\text{th}}$ boundary. 

><code>fp_final_bdry_vertices_ids_CCW_bi.txt</code>: The boundary point IDs in counter-clockwise order of the $b_i^{\text{th}}$ boundary. Each row is a point ID. Suppose the boundary has $n$ points, then there are $n$ rows. That is, the beginning point is not repeated at the end. 
>
><code>fp_final_bdry_vertices_coords_CCW_bi.txt</code>: The boundary point coordinates in counter-clockwise order of the $b_i^{\text{th}}$ boundary. The points have the same ordering and corresponds to the point IDs in <code>fp_final_bdry_vertices_ids_CCW_bi.txt</code>. Each row format is <code>[x y]</code>, the coordinates of the boundary point. Again, the number of rows is $n$, as the beginning point is not repeated at the end. 

Customization
---------------
### Changing number of parallel threads

Sometimes we may get a better performance by setting different numbers of parallel threads for setting up and for the actual meshing. To do so, we can follow the implementation provided in ``primitive_shape_meshing_change_num_threads.cc``. The code structure is similar to the basic example. Here are the few differences we need to make: 

Firstly, we specify two different numbers of parallel threads to use, for setting up and for the actual meshing. 
```c++
    //Specify number of parallel threads for setting up and the actual meshing
    int num_t_setup=4;
    int num_t_mesh=8;
```

In creating the container, shape, sizing field, and the parallel_meshing_2d class object, we use ``num_t_setup``. 

Then, in the parallel meshing part, we call ``mesh_method.change_number_thread(num_t_mesh);`` before the actual meshing, to properly switch to using ``num_t_mesh`` threads. For example, for the DistMesh meshing method,

```c++
      mesh_alg_2d_dm mesh_method(&pm2d);
      //Call this function to properly switch to using num_t_mesh threads for meshing
      mesh_method.change_number_thread(num_t_mesh);
      pm2d.meshing(&mesh_method); 
```


### User-defined sizing field
In the basic example above, the (triangle) element sizing field was calculated automatically, based on the mesh gradation parameter $K$, using the <code>sizing_2d_automatic</code> sizing field object. A possible customization is that users may define their own sizing field function. 

A sizing field function gives the ***relative*** sizes of the triangles. For example, a sizing field of $1+x$ on a square domain $[0,1]\times[0,1]$ means the triangles close to the right edge $x=1$ are about two times larger than the triangles close to the left edge $x=0$. 

 An example implementation is provided in the file ``primitive_shape_meshing_sizing_function.cc``. Notice that we first need to create our own sizing field class ``sizing_2d_func``, which is derived from the base class ``sizing_2d``. In ``sizing_2d_func``, we can define our own sizing field function in ``double getSizingVal(double x, double y)``. 

 ```c++
 /**
 * @brief A class representing a user-defined sizing function in 2D.
 *
 * This class derives from the base class `sizing_2d` and provides user-defined implementations for
 * computing sizing and density values at a given point.
 */
class sizing_2d_func : public sizing_2d {
    public:

      /**
       * @brief Constructor for the `sizing_2d_func` class.
       *
       * @param shp_ A pointer to the underlying shape_2d object.
       */
      sizing_2d_func(shape_2d *shp_) : sizing_2d(shp_){}

      /**
       * @brief Destructor for the `sizing_2d_func` class.
       */
      ~sizing_2d_func(){};

      /**
       * @brief Computes the sizing value at a given point (x, y).
       *
       * This function implements the user-defined sizing function.
       *
       * @param x The x-coordinate of the point.
       * @param y The y-coordinate of the point.
       * @return The computed sizing value at the point (x, y).
       */
      double getSizingVal(double x, double y){
        return 1.0+3*x+3*y;
      };

};
 ```
 
 The execution code in ``int main()`` are the same as previous, except that we use ``sizing_2d_func size_field(&shp);`` for the user-defined sizing field, instead of the automatic sizing field ``sizing_2d_automatic size_field(&shp,K);``. 
 

 The example code generates the following mesh, where the triangles close to the origin are the smallest, and gradually grows larger towards $(1,1)$. 

<p align="center">
<img src="/docs/rectangle_user_defined_sizing_field.png" width="300" />
</p>



### Shape boolean operations
We can perform boolean operations on shapes. Suppose $d_{A}$ and $d_{B}$ are the signed distance fields for two shapes, $A$ and $B$, then the signed distance field from their boolean operations can be computed as follows [[1]](#references):

> Union: $\quad d_{A\cup B}(x,y)=\min (d_{A}(x,y), d_{B}(x,y))$
>
> Difference: $\quad d_{A\backslash B}(x,y)=\max (d_{A}(x,y), -d_{B}(x,y))$
>
> Intersection: $\quad d_{A\cap B}(x,y)=\max (d_{A}(x,y), d_{B}(x,y))$

These operations are implemented in the code as ``shape_2d`` derived class objects that we can create, and they are ``shape_2d_union``, ``shape_2d_difference`` and ``shape_2d_intersection``, respectively. 

An example implementation is provided in the file ``primitive_shape_meshing_boolean_sdf.cc``. The final mesh looks as follows: 

<p align="center">
<img src="/docs/boolean_mesh_plot.png" width="400" />
</p>

To create this shape, we first create four circles with some overlapping regions. And we do intersections for pairs of circles, to obtain the four leafy shapes.  

```c++
//A. Intersection
shape_2d_circle cdl(con,num_t,0.15,0.4,0.4);
shape_2d_circle cdr(con,num_t,0.15,0.6,0.4);
shape_2d_circle cul(con,num_t,0.15,0.4,0.6);
shape_2d_circle cur(con,num_t,0.15,0.6,0.6);

shape_2d_intersection leafd(con,num_t,&cdl,&cdr);
shape_2d_intersection leafl(con,num_t,&cdl,&cul);
shape_2d_intersection leafu(con,num_t,&cul,&cur);
shape_2d_intersection leafr(con,num_t,&cdr,&cur);
```

Next, we create a circle in the middle, and perform union of all the shapes created. 
```c++
//B. Union
shape_2d_circle cm(con,num_t,0.1,0.5,0.5);

shape_2d_union m0(con,num_t,&leafd,&leafl);
shape_2d_union m1(con,num_t,&m0,&leafu);
shape_2d_union m2(con,num_t,&m1,&leafr);
shape_2d_union m3(con,num_t,&m2,&cm);
```

Lastly, we create a large square, and take the difference of the square and the previous union shape. 
```c++
//C. Difference
shape_2d_rectangle rec(con,num_t,0.1,0.9,0.1,0.9);

shape_2d_difference shp_final(con,num_t,&rec,&m3);
```


### Custom shape from user-defined signed distance function (SDF)
We can also define our own SDF. The example code is provided in ``primitive_shape_meshing_own_sdf.cc``. The code generates mesh for a six-pointed hexagram star shape with a hollow circle in the middle -- like a flower: 

<p align="center">
<img src="/docs/shape_user_defined_SDF_mesh.png" width="400" />
</p>

Let's see how to obtain our hexagram star shape. Firstly, we need to create our own SDF class, which is a derived class from the base class ``shape_2d``. 
```c++
/**
 * @brief A class representing a user-defined signed distance function in 2D. 
 * Here, we construct a SDF for a hexagram star shape. 
 * The code of the SDF is adapted from https://iquilezles.org/articles/distfunctions2d/. 
 *
 * This class derives from the base class `shape_2d` and provides user-defined implementations for
 * computing signed distance value at a given point.
 */
class shape_2d_mySDF : public shape_2d{...}
```

For the class constructor, it needs to take in the parameters ``container_2d &con_`` and ``int num_t_``, since they are needed for the initialization of the ``shape_2d`` class. The rest of the input parameters are user-created and based on the shape one defines. Here, we require input parameters ``(centerX_, centerY_)`` where the hexagram centers, and ``radius`` defining the radius of the circumscribed circle of the hexagram. In the constructor, ``get_geometryGrid();`` needs to be called, to initialize and create an underlying geometry grid required for the meshing code. The geometry grid categorize each grid cell as inside the shape, outside the shape, or on the shape boundary, significantly speed up the meshing procedure. 

```c++
    public: 
        /**
         * Constructor for shape_2d_mySDF, a shape given by user-defined signed distance field.
         *
         * @param con_ The container_2d object.
         * @param num_t_ Number of parallel threads.
         * @param centerX_ Center of the hexagram shape, x-coordinate.
         * @param centerY_ Center of the hexagram shape, y-coordinate.
         * @param radius_ Radius of the circumscribed circle.
         */
         shape_2d_mySDF(container_2d &con_, int num_t_,double centerX_, double centerY_, double radius_)
            :shape_2d(con_,num_t_), centerX(centerX_), centerY(centerY_), radius(radius_)
             {
                get_geometryGrid();
             }

        /**
          * Destructor for shape_2d_mySDF.
          * Frees any resources allocated by the object.
          */
         ~shape_2d_mySDF(){};

        double centerX; //Center of the hexagram shape, x-coordinate.
        double centerY; //Center of the hexagram shape, y-coordinate.
        double radius; //Radius of the circumscribed circle.
```

Next, we can define our own SDF in the function ``double sdf(double x, double y)``. The function format must be such that it takes in the coordinates of a point, $(x,y)$, and returns the signed distance to the shape. The SDF calculation of the hexagram shape below is adapted from the code by [[Inigo Quilez ]](https://iquilezles.org/articles/distfunctions2d/). 
```c++
         /**
          * Signed distance function calculation of a hexagram star shape. 
          * Code adapted from https://iquilezles.org/articles/distfunctions2d/. 
          * 
          * @param x The x-coordinate.
          * @param y The y-coordinate.
          * @return The signed distance of the point (x, y) to the shape.
          *         Negative values indicate inside, positive values indicate outside.
          */
         double sdf(double x, double y){
            double kx=-0.5;
            double ky=0.86602540378;
            double kz=0.57735026919;
            double kw=1.73205080757;
            //Translate (x,y) to be centered around the origin, since the hexagram SDF code is centered at the origin. 
            x-=centerX; 
            y-=centerY; 
            double px=abs(x);
            double py=abs(y);
            double fac1=min((kx*px+ky*py),0.0)*2.0;
            double fac2=min((ky*px+kx*py),0.0)*2.0;
            px-=fac1*kx;
            py-=fac1*ky;
            px-=fac2*ky;
            py-=fac2*kx;
            double fac3=min(max(px,radius*kz),radius*kw);
            px-=fac3;
            py-=radius;
            double d=sqrt(px*px+py*py);
            if(py<0.0){
                d=-d;
            }
            return d;

         };
```

Next, in ``int main() `` implementation, we can create our ``shape_2d_mySDF`` shape. Furthermore, we create a circular shape by using the built-in primitive ``shape_2d_circle``. And we perform boolean difference of the hexagram and the circle to obtain our final flower shape. 

```c++
    //User-defined SDF shape
    shape_2d_mySDF myShp(con,num_t,0.5,0.5,0.2);
    //Circular shape
    shape_2d_circle cir(con,num_t,0.1,0.5,0.5);
    //Boolean difference of the two shapes
    shape_2d_difference shp(con,num_t,&myShp,&cir);
```

### Custom shape from contour line segments
We can also define our shape using the contour line segments of the shape. The line segments are stored in a variable ``std::vector<std::vector<double>> boundaries``. 

* Number of shape boundaries : A shape can have multiple shape boundaries. 
> If the shape has one boundary, the line segments of that boundary are stored in a vector, e.g. ``boundary1``, and then further stored in the variable ``boundaries``. That is, ``boundaries[0]`` is that boundary vector.

> If the shape has multiple boundaries, each boundary has a vector that stores its line segments, e.g. ``boundary1``, ``boundary2``, etc. These boundary vectors are further stored in the variable ``boundaries``. That is, ``boundaries[0]`` is the vector of the first boundary, ``boundaries[1]`` is the vector of the second boundary, and so on. 

* Orientation: 
> A boundary enclosing the shape interior is defined as a sequence of points in ***clockwise*** order, forming a closed loop with the same end point ``[x0,y0]``. For example, ``boundary1`` can be defined as: ``boundary1 = [x0,y0, x1,y1, x2,y2, x3,y3, x0,y0]``, where the ***interior*** of the loop is the ***interior*** of the shape.

> If the shape has holes, the hole boundary is inputted as a sequence of points in  ***counter-clockwise*** order, again forming a close loop. In this case, the ***interior*** of the loop is the ***exterior*** of the shape.

#### ***Example 1: Custom poker***

An example is provided in ``custom_poker_meshing.cc``. Here is the mesh it produces:

<p align="center">
<img src="/docs/custom_shape_contour_line_seg.png" width="400" />
</p>

Here, we have one boundary. First, we have ``boundaryx`` and ``boundaryy`` storing the $x$ and $y$ coordinates of the shape boundary. Notice that the coordinates form a closed loop with the same end point. In a ``for`` loop, we store the coordinates in the boundary vector ``boundary1``, in the format ``boundary1 = [x0,y0, x1,y1, ..., x0,y0]``. We also keep a count of the number of boundary line segments ``bdry_ct``. 

Then, we store ``boundary1`` into the viariable ``boundaries``. 

```c++
//Define a shape from custom contour line segments
    double boundaryx[117]={ 0.585 , ... , 0.585};
    double boundaryy[117]={ 0.0284, ... , 0.0284};
    std::vector<std::vector<double>> boundaries;
    std::vector<double> boundary1;
    int bdry_ct=0;
    for(int i=116;i>=0;i--){
        boundaryx[i]=boundaryx[i];
        boundaryy[i]=boundaryy[i];
        boundary1.push_back(boundaryx[i]);
        boundary1.push_back(boundaryy[i]);
        bdry_ct++;
    }
    boundaries.push_back(boundary1);
```

In this code, we also distinguish a few different parallel number of threads of use. Note that this is optional - for simplicity, one can as well just use the same number of parallel threads throughout the code. 

``num_t_meshing`` is the number of threads to use during the actually meshing procedure. ``num_t_setup`` is the number of threads to use when we are setting up the shape, container, and sizing field, etc, before the actual meshing. ``num_t_setup_cshp`` is the number of threads to use specifically for the custom shape creation. If there are only a small number of boundary line segments ``bdry_ct``, then we may set ``num_t_setup_cshp=1``, to avoid parallel overhead. But when there are a large number of boundary line segments, then we may set ``num_t_setup_cshp=num_t_setup``. Here, an arbitraty threshold of $200$ is used. 

```c++
      int num_t_meshing=2*physical_core;
      int num_t_setup=physical_core;
      //for custom_shape_2d: if the boundary contour line segments count is small, we use serial code nt=1 to avoid parallel overhead
      int num_t_setup_cshp=1; 
      if(bdry_ct>200){
            num_t_setup_cshp=num_t_setup;
      }
```

In the shape creation, we can create an ``shape_2d_contour_lines`` for our custom shape defined by the shape boundaries. Let's first take a look at parts of the constructor of the class below. We see that it is a derived class from the base class ``shape_2d``, therefore, it needs to take as input parameters ``container_2d &con_`` and ``int num_t_``, which are used to initiate the base class ``shape_2d``. Furthermore, an input parameter ``int num_t_cshp_`` is needed to create the custom shape. The variable ``std::vector<std::vector<double>> boundaries_`` is also required, which has the shape boundaries information. Lastly, a boolean variable ``bool normalize_model_`` is required, which indicates whether to normalize the model. If true, the model will be translated to fit and centered into the container domain, and the model will be scaled. Suppose the container domain size is ``[ax, bx] x [ay, by]``, with side length ``Lx = bx-ax``, ``Ly = by-ay``. The model will be scaled so that its longest dimension fits in ``[ax+0.1*Lx, ax+0.9*Lx] x [ay+0.1*Ly, ay+0.9*Ly]``. 

```c++
/**
* Constructor for shape_2d_contour_lines, a custom shape defined by user-input contour line segments.
*
* @param con_ The container_2d object.
* @param num_t_ Number of parallel threads.
* @param num_t_cshp_ Number of parallel threads for the custom shape.
* @param boundaries_ A vector of vectors representing the contour line segments of the shape.
* @param normalize_model_ Flag indicating whether to normalize the model.
*/
shape_2d_contour_lines(container_2d &con_, int num_t_, int num_t_cshp_, std::vector<std::vector<double>> boundaries_, bool normalize_model_):shape_2d(con_,num_t_),
cshp(custom_shape_2d(boundaries_, con_.ax, con_.bx, con_.ay, con_.by, normalize_model_,num_t_cshp_))
{...}
```

Based on the constructor description, we can create our custom shape as follows. We first feed in the ``boundaries`` to the ``shape_2d_contour_lines`` object, and require it to normalize the model. We also make a hollow circule inside the shape. 

```c++
      //final test shape, complicated poker shape
      bool normalize_model=true;
      shape_2d_contour_lines shp1(con,num_t_setup,num_t_setup_cshp,boundaries,normalize_model);
      shape_2d_circle shp2(con,num_t_setup,0.1,0.5,0.7);
      shape_2d_difference shp(con,num_t_setup,&shp1,&shp2);
```

The following code are similar. Except that, since we want to use a different number of parallel threads for the meshing, we use the function ``mesh_method.change_number_thread(num_t_meshing)`` to properly set the number of threads in the code before doing the actual meshing. 

```c++
if(method_ind==0){
         mesh_alg_2d_dm mesh_method(&pm2d);
         mesh_method.change_number_thread(num_t_meshing);
         pm2d.meshing(&mesh_method); 
    }
    else if(method_ind==1){
         mesh_alg_2d_cvd mesh_method(&pm2d);
         mesh_method.change_number_thread(num_t_meshing);
         pm2d.meshing(&mesh_method); 
    }
    else if(method_ind==2){
         mesh_alg_2d_hybrid mesh_method(&pm2d);
         mesh_method.change_number_thread(num_t_meshing);
         pm2d.meshing(&mesh_method); 
    }
```

#### ***Example 2: North America geography map***

Another example is meshing on the North America geography map:

![NA mesh plot](/docs/north_america_map_mesh_plot.jpg)

The mesh is generated using the file ``North_America_map_meshing.cc``. To obtain the shape, we first read in the contour line segments in the file ``north_america_single_bdry_xy.txt``, and store them in the variable ``boundaries``. Each **line** in ``north_america_single_bdry_xy.txt`` gives a boundary of the shape, in the format: ``[x0 y0 x1 y1 ... xk yk ... x0 y0]``

After we define the shape from the contour line segments, we follow the similar procedure as before to obtain the mesh. 

A few comments:

- The time it takes to initialize the shape is relatively long here. This is because, in creating the shape, we build an underlying geometry grid, where each grid cell is categorized as inner/boundary/outer grid. In this process, each grid cell need to computer their center signed distance, which is expensive, since we have a highly detailed geometry. The calculation will be improved in the near future. We will implement a quad-tree structure for the geometry grid, to avoid computing the signed distance for every grid cell. 

- The mesh has good quality overall. In this testing, we obtain $19851$ vertices and $35972$ triangles, with median aspect and edge ratios to be $1.03383$ and $1.2194$, respectively. We see some meshing inaccuracy along the shape boundary, up to the resolution of the local triangle sizes. 


### Setting fixed points in the mesh

We can also add a list of fixed points into the mesh. These points stays at the some positions throughout the meshing procedure. An example implementation is provided based on the custom poker shape meshing code, in ``custom_poker_add_fixed_pts_meshing.cc``. 


Once we created the ``parallel_meshing_2d pm2d`` class object, we add the fixed points in. As shown in the following, we first create a list of the fixed point coordinates, with format ``[x0 y0 x1 y1 ...]``. To add them into the mesh, 
> One option is to call ``pm2d.add_fixed_points_normailze(Nfixed, &shp1, fixed_pt_list)``, which perform the same normalization (scale + translation) on the fixted points as the ``shape_2d`` object ``shp1``. Note that only shapes built from shape contour line segments, ``shape_2d_contour_lines`` class objects, has the normalization object at initialization. Therefore, this option to normalize the fixed points is only relevant for ``shape_2d_contour_lines`` shape objects. 

> Another option is to use ``pm2d.add_fixed_points(Nfixed, fixed_pt_list)``, which simple adds the points as they are into the mesh as fixed points. 

```c++
     //Add fixed points into the mesh
     int Nfixed=10; //Number of fixed points
     //Coordinates of the fixed points stored into a list
     double *fixed_pt_list=new double[2*Nfixed];
     double fixed_pt_x[Nfixed]={0.2,0.3,0.4,0.2,0.2,0.5,0.5,0.7,0.7,0.7};
     double fixed_pt_y[Nfixed]={0.2,0.3,0.4,0.3,0.4,0.6,0.1,0.4,0.3,0.2};
     //The fixed point list has format [x0, y0, x1, y1, ...]
     for(int i=0;i<Nfixed;i++){
        fixed_pt_list[2*i]=fixed_pt_x[i];
        fixed_pt_list[2*i+1]=fixed_pt_y[i];
     }
     //Option 1: Use the same normalization (scale + translation) on the fixed points as in shp1
     pm2d.add_fixed_points_normailze(Nfixed, &shp1, fixed_pt_list);
     //Option 2: Simple add in the points as they are as fixed points
     //pm2d.add_fixed_points(Nfixed, fixed_pt_list);
```

After adding fixed points, we can initialize meshing points as follows. Notice that here, the total number of meshing points, ``Ntotal`` equals to the sum of the number of fixed points and the number of unfixed points. 
```c++
     pm2d.pt_init(Ntotal); //Ntotal=Nfix+Nmove: total number of fixed and unfixed points in the mesh
```

Performance
================================================

For detailed performance analysis, please refer to Sec. 6 and 7 of the [paper](https://arxiv.org/abs/2309.13824). 
Here, we provide some brief summaires on two cases:

- [A. Square, uniform sizing field](#a-square-uniform-sizing-field)

- [B. Square, uniform sizing field](#b-custom-poker-adaptive-sizing-field)

A few aspects we test are: 

- **Parallel performance**. We measure the wall clock time $t_p$ using a variable number of $p$ threads, and calculate the parallel efficiency using
    $$T_e(p) =\frac{t_1}{p \cdot t_p}$$
    which measures the effective slowdown from the hypothetically perfect parallel
    scaling.

- **Optimal point scheme**. The above tests are repeated using two different point schemes: 
      
      I. All ``Ntotal`` points are initialized prior to meshing.

      II. A small fraction of points are first initialized, and during the meshing iterations, additional points are inserted to refine the mesh, until ``Ntotal`` points are reached.
    We compare the timing performance and mesh quality of the point schemes.

- **Mesh quality**. We provide mesh quality statistics from a typical run using the preferred point scheme.
   
- **System size**. Using the preferred point scheme, we vary the number of meshing points ``Ntotal``, and examine the parallel performance of the code as system size varies. 

Our tests are done on a Ubuntu Linux computer with 256 GB of memory and dual Intel Xeon E5-2650L v4 processors with 14 low-power cores using a 1.7 GHz base clock speed. Since the computer has $28$ physical cores, we focus our attention to the first $28$ threads. 


## A. Square, uniform sizing field. 
Here is a schematic plot of an example mesh with ``Ntotal``$=10^3$ points.

<p align="center">
<img src="/docs/square_mesh_1000.png" width="300" />
</p>

### Recommended settings
Uniform meshing of a square represents a situation where both the geometry and the sizing field of the mesh are simple. In general, we found that for ***simple and convex shapes*** (e.g. rectangles, circles) where points can easily move around without obstruction, and when the ***element sizing field varies gradually***, there is no need to use the point addition scheme. 
>Instead, ***initializing all points to start meshing*** gives better performance. 
>>This can be done by modifying the parameter ``const double pt_init_frac=1.0`` in file ``/src_TriMe/config.hh``.

>Furthermore, we recommend to use ***CVD meshing or hybrid meshing*** in these situations. 

### Performance
Here, we show performance statistics (parallel time, mesh quality) with ``Ntotal``$=10^6$ points, of the three meshing methods, using the recommended setting of initializing all points prior to meshing:

![uniform-square-performance](/docs/Github-performance-uniform-square.png)

 Plot (a) shows that parallel computation significantly speed up the meshing time. Indeed, the serial codes run close to $300$ s for both CVD and hybrid meshing. But the $28$-thread parallel code run only $10\sim 15$ s. Here, we see DistMesh takes much longer. This can sometimes happen when a specific set of meshing points results in problematic triangles during the meshing process, causing more iterations to terminate. This situation rarely happens using CVD meshing. Since the shape and the sizing field are both simple, CVD meshing can effectively and steadily improve mesh quality and stably lead to terminate. Furthermore, as hybrid meshing switches to CVD meshing once the overall mesh achieves a good quality, it also stably leads to termination. 
 
 Plot (b) shows the corresponding parallel efficiency. As expected, CVD meshing has very high parallel efficiency close to $0.9$, since in the CVD procedure, the Voronoi cell centroids can be computed simultaneously among parallel threads. And we see that DistMesh has the lowest efficiency, as the algorithm is less efficient for parallelization. Hybrid meshing has parallel efficiency in between the other two. 

 Plot (c) shows the median aspect ratios of the meshes against the number of meshing iterations from a typical run. We see all methods produce overall very high quality meshes at termination.
 
Plot (d) shows the maximum aspect ratio, corresponding to the worst triangles in the meshes. We see the ratios are very reasonable for all three methods, meaning all triangles are of high quality and close to equilateral triangles. 


### Parallel performance for varying system size

<p align="center">
<img src="/docs/Github-performance-uniform-square-Npt.png" width="600" />
</p>

The plot here shows the parallel efficiency as the number of points increases. As expected, CVD meshing has the best efficiency, and the efficiency of DistMesh is the lowest. Furthermore, we see a gradually increasing trend of parallel efficiency for the three meshing methods from $10^4$ to $10^6$ particles. After this, parallel efficiencies remain roughly at the same level as the system size continues to increase. The higher parallel efficiency for large systems is desirable for large-scale meshing. 

## B. Custom poker, adaptive sizing field. 
This is the case illustrated in the [custom poker example](#example-1-custom-poker). For performance test, we mesh it with ``Ntotal``$=10^6$ points. The corresponding point density field is shown here: 

<p align="center">
<img src="/docs/large_system_density_field.png" width="400" />
</p>

### Recommended settings
Adaptive meshing on the custom shape represents a case where the ***geometry and the sizing field are complicated***. In general, we found that when either geometry or sizing field is complicated, it is recommended to use: 

>***Point addition scheme***, where we initialize ``pt_init_frac`` of ``Ntotal`` points to start meshing. And we add points during meshing from time to time to refine regions of the mesh that need points the most.   
>> This can be done by modifying the parameter ``const double pt_init_frac=0.2`` in file ``/src_TriMe/config.hh``.
>>
>> Reason: According to the highly non-uniform density field, points in the mesh require more movements and adjustments. However, since the geometry is complicated, points cannot move around freely without obstruction. Using the point addition scheme effectively refines regions of the mesh that need points most during the meshing iterations, and leads to much better shaped triangles conforming to the density field. 

>Furthermore, we recommend to use ***hybrid meshing*** in these situations. 


### Performance
Here, we show performance statistics (parallel time, mesh quality) with ``Ntotal``$=10^6$ points, of the three meshing methods, using the recommended point addition scheme:

![adaptive-poker-performance](/docs/Github-performance-adaptive-poker.png)

Plot (a) shows parallelization significantly speed up the computation time. Furthermore, we see that CVD meshing is the most expensive in serial code, while DistMesh and hybrid meshing takes similar time. This is reasonable, since the Voronoi cell centroid calculation is more expensive than the force balancing procedure in DistMesh. And because hybrid meshing uses mostly the cheaper DistMesh iterations, its computation time is similar to that of DistMesh, and hybrid also tends to terminate stably like CVD. Using $28$-thread parallel code, the computation time of the three methods are similar.

Plot (b) shows the corresponding parallel efficiency. As expected, CVD meshing has high parallel efficiency, DistMesh the lowest, and hybrid in between. In this case, hybrid meshing has only a few CVD iterations at the end, so its efficiency is similar to that of DistMesh. 

Plot (c) shows the median aspect ratio. The peaks in the plots correspond to when points are added into the mesh according to the point addition scheme. All three meshing methods achieve very good median ratios at termination.

Plot (d) shows the maximum aspect ratios of the worst triangles. DistMesh typically has highest maximum aspect ratios, corresponding to flatter shaped triangles in the mesh. In comparison, CVD meshing has better shaped triangles and lower maximum aspect ratio. Since hybrid meshing uses CVD iterations as the end as refinement steps, it also has lower maximum aspect ratio similar to that of CVD meshing.

### Parallel performance for varying system size
As shown in the plot, similar to the uniform meshing of sqaure case, we observe similar trend of a gradual increase in parallel efficiency as the system size increases from $10^4$ to $10^6$ particles. Then the efficiency remains at roughly the same level as the system size continues to increase. Furthermore, CVD meshing has the highest efficiency, DistMesh the lowest, and hybrid meshing in between. 

<p align="center">
<img src="/docs/Github-performance-adaptive-poker-Npt.png" width="600" />
</p>



Code updates
================================================

### May 19th, 2024: Delete problematic vertices at termination
1. At termination, we check for any problematic vertices, defined as: 
      - Single vertice not incident to any triangles in the mesh. 
      - Vertices incident to any single triangle not connected to any other triangles. 
      - Vertices connected to more than two boundary edges. 

      In a while loop, we search and delete these vertices, until no more problematic vertices are left. 

      Therefore, the final mesh may have $\leq$ ``Ntotal`` number of vertices that we specified. This choice is made for better mesh qualities. 

### Oct 19th, 2023: Half-edge data structure; Fixed points input

1. Added in fixed point input for meshing with fixed points.

2. Implemented half-edge structure for final mesh clean up and outputting boundary vertices in counter-clockwise (CCW) order.

3. Triangle vertex ID's are outputted in counter-clockwise (CCW) order at termination.


Acknowledgement
=====
This research was supported by a grant from the United States-Israel Binational
Science Foundation (BSF), Jerusalem, Israel through grant number 2018/170.
C. H. Rycroft was partially supported by the Applied Mathematics
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


