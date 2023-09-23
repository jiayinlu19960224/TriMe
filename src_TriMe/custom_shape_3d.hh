#ifndef CUSTOM_SHAPE_3D_HH
#define CUSTOM_SHAPE_3D_HH

#include "basic_calculation.hh"


namespace voro {
    
    class custom_shape_3d : public basic_calculation_3d{
        //***Read in .stl mesh and get triangles, normals
        //***closest point on mesh, closestx, closesty, closestz
        //***signed distance to the mesh
      public:
        
       //domain container grid information
        const int nx, ny, nz;
        const double ax, bx, ay, by, az, bz;
        const double dx, dy, dz;

       //mesh information
        //read in .stl file and get the mesh/shape info
        custom_shape_3d(std::string file_name, double ax_, double bx_, double ay_, double by_, double az_, double bz_, int nx_, int ny_, int nz_);
        std::vector<std::vector<std::vector<double> > > triangles; // list of triangles, indexed by triangle index
        std::vector<std::vector<double> > normals; // list of normal for each triangle face, index by triangle index
        std::vector<std::vector<double> > edge_normals;//list of normal for each triangle edge, index by edge index
        std::vector<std::vector<double> > ver_normals;//list of normal for each triangle vertice, index by vertex index
        std::vector<std::vector<double> > vertices;// list of vertices
        std::vector<std::vector<int> > elements; // list of three vertex indexes of each triangle, index by triangle index
        std::vector<std::vector<int> > edges; //list of three triangle edges indexes for each triangle, index by triangle index
        std::vector<std::vector<int> > ver_tria; //list of triangles indexes sharing a common vertex, indexed by vertex index
        std::vector<std::vector<int> > edge_tria; //list of triangles indexes sharing a common edge, indexed by edge index
        std::vector<std::vector<int> > cell_tria; //list of triangles in a grid cell; index by grid cell (i,j,k) of the domain with nx, ny, nz cells
        int num_edges; //number of different edges

      //mesh information
        //Read in a .stl file and store the triangles and the triangle normals to "triangles", and "normals"
        void ReadSTL(std::string file_name);
        void normalize_model(std::vector<std::vector<std::vector<double> > >& _triangles);
        //get unique vertex info, vertex coordinates index by vertex indices
        void get_vertex();
        //get triangle elements (list of three vertex index for each triangle)
        void get_elements();
        //get ver_tria list: list of triangles sharing a vertex
        void get_ver_tria();
        //get edge_tria list: list of triangles sharing an edge
        void get_edge_tria();
        //get cell_tria: list of triangles in each grid cell, stored in cell_tria, index by cell i, j, k: cell_tria[nx*ny*k+nx*j+i]
        void get_cell_tria();
        //index all edges for each triangle.
        //generated _edges[i][j] represents the index of the jth edge in ith triangle.
        void CreateEdges();
        //find index of a given vertice (x, y, z)
        int find_index(double x, double y, double z);
        //compare two vertex
        bool vertexCmp(double x0,double y0,double z0,double x1,double y1,double z1);
        bool tria_cube_test3(double ax, double ay, double az, 
            double hx, double hy, double hz, 
            double v0x, double v0y, double v0z, double v1x, double v1y, double v1z, double v2x, double v2y, double v2z);

      //shape distance
        //return the signed distance of a point X(x,y,z) to the shape boundary
        double f_b_pts(double x, double y, double z) const;
        
    };

}


#endif