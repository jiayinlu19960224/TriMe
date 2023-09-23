#ifndef MESH_ALG_2D_CVD_HH
#define MESH_ALG_2D_CVD_HH


#include "mesh_alg_2d.hh"
//#include "parallel_meshing_2d.hh"

namespace voro {

    /**
     * @brief The `mesh_alg_2d_cvd` class is a derived class from `mesh_alg_2d` that implements functions specific to CVD (Centroidal Voronoi Diagram) meshing.
     * It inherits the base class `mesh_alg_2d` and adds CVD-related variables and functions.
     */
    class mesh_alg_2d_cvd :  virtual public mesh_alg_2d {
      public:

        /**
         * @brief Constructs a new `mesh_alg_2d_cvd` object.
         * @param pm2d_ Pointer to the `parallel_meshing_2d` object.
         */
        mesh_alg_2d_cvd(parallel_meshing_2d *pm2d_);

        /**
         * @brief Destroys the `mesh_alg_2d_cvd` object.
         */
        ~mesh_alg_2d_cvd();
        
        //CVD related variables and functions

        /**
         * @brief Calculates the centroid of a Voronoi cell using quadrature on its subdivision triangles.
         * 
         * This function computes the centroid of a Voronoi cell by performing quadrature on its subdivision triangles.
         * It takes the number of triangles in the subdivision, the X and Y coordinates of the triangle vertices,
         * and computes the centroid by evaluating the weighted sum of the triangle centroids.
         * 
         * @param tria_make_ct Number of triangles in the Voronoi cell's subdivision.
         * @param tria_vertices_x X-coordinates of the triangle vertices.
         * @param tria_vertices_y Y-coordinates of the triangle vertices.
         * @param centroid_x_ Computed X-coordinate of the centroid.
         * @param centroid_y_ Computed Y-coordinate of the centroid.
         */
        void calculate_voro_centroid_from_triangles(int tria_make_ct,
            std::vector<double> &tria_vertices_x, std::vector<double> &tria_vertices_y,
            double &centroid_x_, double &centroid_y_);


        /**
         * @brief Calculates the centroid and area of a Voronoi cell using quadrature on its subdivision triangles.
         * 
         * This function computes the centroid and area of a Voronoi cell by performing quadrature on its subdivision triangles.
         * It takes the number of triangles in the subdivision, the X and Y coordinates of the triangle vertices,
         * and computes the centroid by evaluating the weighted sum of the triangle centroids. Additionally, it computes
         * the area of the Voronoi cell using the signed area of each triangle in the subdivision.
         * 
         * @param tria_make_ct Number of triangles in the Voronoi cell's subdivision.
         * @param tria_vertices_x X-coordinates of the triangle vertices.
         * @param tria_vertices_y Y-coordinates of the triangle vertices.
         * @param centroid_x_ Computed X-coordinate of the centroid.
         * @param centroid_y_ Computed Y-coordinate of the centroid.
         * @param voro_area_ Computed area of the Voronoi cell.
         */
        void calculate_voro_centroid_from_triangles(int tria_make_ct,
            std::vector<double> &tria_vertices_x, std::vector<double> &tria_vertices_y,
            double &centroid_x_, double &centroid_y_,double &voro_area_);


        /**
         * @brief Computes the adaptive Voronoi centroid for a particle.
         * 
         * This function calculates the adaptive Voronoi centroid for a particle with the given ID (`pid`). The resulting
         * centroid coordinates are stored in the `cx` and `cy` variables passed by reference. The function uses the provided
         * Voronoi vertex count (`voro_vertex_ct`) and the coordinates of the Voronoi vertices (`voro_vertex_xy`) to perform
         * the computation.
         * 
         * The adaptive Voronoi centroid computation involves subdividing the Voronoi cell into triangles and using quadrature
         * on these triangles to approximate the centroid. The subdivision process continues until the computed centroid
         * satisfies a specified error tolerance.
         * 
         * @param pid The ID of the particle.
         * @param cx Reference to the variable storing the X-coordinate of the computed centroid.
         * @param cy Reference to the variable storing the Y-coordinate of the computed centroid.
         * @param voro_vertex_ct The count of Voronoi vertices for the particle.
         * @param voro_vertex_xy The coordinates of the Voronoi vertices.
         */
        void adaptive_compute_voro_centroid(int pid,double &cx,double &cy,
            int voro_vertex_ct, std::vector<double> &voro_vertex_xy);

        double *prev_pt_mvmt_chrtrt_fac; /**< Factor used in Voronoi centroid error tolerance calculation; Ratio of previous point movement over characteristic edge length at the point. */
        //int *voro_tria_subdivide_ct; //number to triangles to sub-divide a Voronoi cell into, to calcualte its centroid with Trapezoid rule
        //void update_voro_tria_subdivide_ct(int i, int voro_vertex_ct);
       // void compute_voro_centroid(int pid,double &cx,double &cy,
       //       int voro_vertex_ct, std::vector<double> &voro_vertex_xy);

        //virtual functions in base class
        /**
         * Updates the position of a point with the CVD algorithm.
         *
         * @param i The index of the point.
         * @param c The Voronoi cell associated with the point.
         * @param x The new x-coordinate of the point.
         * @param y The new y-coordinate of the point.
         * @param inner_pt_over_stop_Continue_mvmt_thres_ct_ Reference to the counter for inner points over the stop Continue movement threshold.
         * @param use_random_point If using CVD meshing: Decide whether to move point to the centroid of the needle type triangle it connects to instead, to eliminate the needle triangle.
         * @param random_point_x The x-coordinate of the needle type triangle centroid the point should move to.
         * @param random_point_y The y-coordinate of the needle type triangle centroid the point should move to.
         */
        void update_pt_position_cvd(int i, voronoicell_neighbor_2d &c,double x, double y, int &inner_pt_over_stop_Continue_mvmt_thres_ct_, 
            bool use_random_point, double random_point_x, double random_point_y) override;
        
        /**
         * Updates the position of points and determines reTria for CVD.
         * The positions are updated in the xy_id array.
         * reTria is set to true for CVD.
         */
        void update_pt_position_and_determine_reTria() override; 
        
        /**
         * Computes the Voronoi cell and extracts information (Voronoi vertices). 
         * This function is specific to the CVD algorithm.
         */
        void voro_compute_retria_extract_info() override; 

    };

    /**
     * @brief The cvd_tria_info class represents a subdivision triangle of a Voronoi cell and computes related information.
     */
    class cvd_tria_info{
      public:

        /**
         * @brief Constructs a cvd_tria_info object with the given coordinates and associated parallel_meshing_2d object.
         *
         * @param x1_ x-coordinate of the first vertex of the triangle.
         * @param y1_ y-coordinate of the first vertex of the triangle.
         * @param x2_ x-coordinate of the second vertex of the triangle.
         * @param y2_ y-coordinate of the second vertex of the triangle.
         * @param x3_ x-coordinate of the third vertex of the triangle.
         * @param y3_ y-coordinate of the third vertex of the triangle.
         * @param pm2d_ Pointer to the associated parallel_meshing_2d object.
         */
        cvd_tria_info(double x1_,double y1_,double x2_,double y2_,double x3_,double y3_,parallel_meshing_2d *pm2d_);
        parallel_meshing_2d *pm2d; /**< Pointer to the associated parallel_meshing_2d object. */
        double x1(){return _x1;}; /**< Get the x-coordinate of the first vertex. */
        double y1(){return _y1;}; /**< Get the y-coordinate of the first vertex. */
        double x2(){return _x2;}; /**< Get the x-coordinate of the second vertex. */
        double y2(){return _y2;}; /**< Get the y-coordinate of the second vertex. */
        double x3(){return _x3;}; /**< Get the x-coordinate of the third vertex. */
        double y3(){return _y3;}; /**< Get the y-coordinate of the third vertex. */
            
        double area(){return _area;}; /**< Get the area of the triangle. */
        double integrand_rho(){return _integrand_rho;}; /**< Get the density and triangle area factor */
        double integrand_xrho(){return _integrand_xrho;}; /**< Get the density and triangle area weighted x-coordinate of the triangle centroid. */
        double integrand_yrho(){return _integrand_yrho;}; /**< Get the density and triangle area weighted y-coordinate of the triangle centroid. */

        /**
         * @brief Computes the density at the given coordinates using the associated parallel_meshing_2d object.
         *
         * @param x x-coordinate.
         * @param y y-coordinate.
         * @return The density value at the given coordinates.
         */
        double density(double x, double y);
        
        double _x1; /**< x-coordinate of the first vertex. */
        double _y1; /**< y-coordinate of the first vertex. */
        double _x2; /**< x-coordinate of the second vertex. */
        double _y2; /**< y-coordinate of the second vertex. */
        double _x3; /**< x-coordinate of the third vertex. */
        double _y3; /**< y-coordinate of the third vertex. */


        double _area; /**< Area of the triangle. */
        double _integrand_rho; /**< Density and triangle area factor */
        double _integrand_xrho; /**< Density and triangle area weighted x-coordinate of the triangle centroid. */
        double _integrand_yrho; /**< Density and triangle area weighted y-coordinate of the triangle centroid. */
        
        /** 
         * @brief Computes the information of a subdividing triangle of a Voronoi cell. 
         * It calculates various properties of the triangle, including its area, density, 
         * density and triangle area weighted centroid, and related integrands.
         */
        void compute_info();
    };
}



#endif



