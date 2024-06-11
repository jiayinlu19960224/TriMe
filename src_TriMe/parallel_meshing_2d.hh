/**
 * @file parallel_meshing_2d.h
 * @brief Header file for the parallel_meshing_2d class.
 */

#ifndef PARALLEL_MESHING_2D_HH
#define PARALLEL_MESHING_2D_HH

#include "shape_2d.hh"
#include "adf_2d.hh"
#include "sizing_2d.hh"
#include "mesh_alg_2d.hh"
#include "mesh_alg_2d_dm.hh"
#include "mesh_alg_2d_cvd.hh"
#include "mesh_alg_2d_hybrid.hh"

#include <algorithm>
#include <limits>
#include <vector>
#include "voro++.hh"
#include <fstream>


namespace voro {

    /**
     * @brief A class for parallel 2D meshing calculations.
     *
     * The `parallel_meshing_2d` class extends the `basic_calculation_2d` class and provides functionality for parallel meshing operations. It includes various member variables and functions for setting up the meshing domain, handling geometry grids, initializing points, performing meshing algorithms, and timing measurements.
     *
     * The class supports the following features:
     * - Setting up the meshing domain and geometry grids
     * - Handling shape, container, and sizing field objects
     * - Computing signed distance functions (SDF), sizing values, and density values
     * - Handling point categories and characteristic element size lengthscale (h)
     * - Constructing adaptive arrays for boundary grids
     * - Error diffusion for point initialization
     * - Printing initial points to a file
     * - Obtaining geometric tolerances for boundary grids using adaptive adf quad cells
     * - Meshing algorithms for generating the mesh
     * - Timing measurements for generating points in a grid
     *
     * The class is designed to be used in parallel meshing calculations in a 2D domain.
     */
    class parallel_meshing_2d : public basic_calculation_2d {
    public:
        /**
         * @brief Constructor for the parallel_meshing_2d class.
         * @param con_ Pointer to the container_2d class.
         * @param shp_ Pointer to the shape_2d class.
         * @param size_field_ Pointer to the sizing_2d class.
         * @param num_t_ Number of threads to use in parallel meshing.
         * @param output_interval_ Output interval for meshing.
         * @param file_name_prefix_ Output file name prefix.
         */
        parallel_meshing_2d(container_2d* con_, shape_2d* shp_, sizing_2d* size_field_, int num_t_, int output_interval_,const char *file_name_prefix_);
        
        /**
         * @brief Destructor for the parallel_meshing_2d class.
         */
        ~parallel_meshing_2d();

         /**
         * @brief Check the point category.
         */
        void check_pt_ctgr();
        //------------------Skeleton functions and variables---------------------

        // A. setting up
        const char* file_name_prefix; /**< Output file name prefix */
        int num_t; /**< Number of threads to use in parallel meshing */


        // Meshing domain variables
        double ax; /**< X-coordinate of the domain lower-left corner */
        double bx; /**< X-coordinate of the domain upper-right corner */
        double ay; /**< Y-coordinate of the domain lower-left corner */
        double by; /**< Y-coordinate of the domain upper-right corner */
        double lx; /**< Length of the domain in the x-direction */
        double ly; /**< Length of the domain in the y-direction */

        
        // Meshing geometry grid variables
        int gnx;    /**< Number of grid points in the x-direction */
        int gny;    /**< Number of grid points in the y-direction */
        int gnxy;   /**< Total number of grid points */
        double gdx; /**< Grid spacing in the x-direction */
        double gdy; /**< Grid spacing in the y-direction */
        double diag_gdxy; /**< Diagonal length of a grid cell */
        double inv_gdx; /**< Inverse of the grid spacing in the x-direction */
        double inv_gdy; /**< Inverse of the grid spacing in the y-direction */
        double gdxy; /**< Product of gdx and gdy: gdxy = gdx * gdy */
        int geo_bgrid_ct; /**< Number of geometry boundary grids */
        int geo_igrid_ct; /**< Number of geometry inner grids */

        shape_2d* shp; /**< Pointer to the shape_2d class providing geometry information */
        container_2d* con; /**< Pointer to the container_2d class */
        sizing_2d* size_field; /**< Pointer to the sizing_2d class providing element sizing information */

        /**
         * @brief Get the signed distance function (SDF) value at a given point.
         * @param x X-coordinate of the point.
         * @param y Y-coordinate of the point.
         * @return SDF value at the point (x, y).
         */
        double sdf(double x, double y) { return shp->sdf(x, y); }

        /**
         * @brief Get the sizing value at a given point.
         * @param x X-coordinate of the point.
         * @param y Y-coordinate of the point.
         * @return Sizing value at the point (x, y).
         */
        double sizing(double x, double y) { return size_field->getSizingVal(x, y); }

        /**
         * @brief Get the density value at a given point.
         * @param x X-coordinate of the point.
         * @param y Y-coordinate of the point.
         * @return Density value at the point (x, y).
         */
        double density(double x, double y) { return size_field->getDensityVal(x, y); }

        /**
         * @brief Get the geometry grid value at a given point.
         *
         * The value represents the grid category:
         * - Inside cells: -1, -2, -3
         * - Boundary cells: 1, 2, 3
         * - Outside cells: gnxy+1, gnxy+2, ...
         *
         * @param x X-coordinate of the point.
         * @param y Y-coordinate of the point.
         * @return Geometry grid value at the point (x, y).
         */
        int geo_grid(double x, double y) { return shp->pt_geo_grid_val(x, y); }

        /**
         * @brief Get the geometry grid value at a given grid index.
         * 
         * * The value represents the grid category:
         * - Inside cells: -1, -2, -3
         * - Boundary cells: 1, 2, 3
         * - Outside cells: gnxy+1, gnxy+2, ...
         * 
         * @param ij Grid index.
         * @return Geometry grid value at the grid index.
         */
        int geo_grid(int ij) { return shp->geo_grid[ij]; }

        /**
         * @brief Get the boundary grid index (ij) from a given index (ind).
         * @param ind Boundary grid index.
         * @return Boundary grid index (ij).
         */
        int geo_bgrid_ij(int ind) { return shp->geo_bgrid_ij[ind]; }

        /**
         * @brief Get the inner grid index (ij) from a given index (ind).
         * @param ind Inner grid index.
         * @return Inner grid index (ij).
         */
        int geo_igrid_ij(int ind) { return shp->geo_igrid_ij[ind]; }

        /**
         * @brief Get the geometric tolerance for boundary grids at a given point.
         * @param x X-coordinate of the point.
         * @param y Y-coordinate of the point.
         * @return Geometric tolerance for the boundary grids at the point (x, y).
         */
        double get_bgrid_geps(double x, double y);


        
        /**
         * @brief Calculate the spatial step size (delta x, delta y) for use in finite difference for calculating spatial derivatives in the projection step.
         *
         * If the point (x, y) is in the boundary grid, this function calculates the spatial step size (deps) using the boundary grid information.
         * If the point is not in the boundary grid, it returns the previously calculated valid spatial step size (deps_prime).
         *
         * @param x X-coordinate of the point.
         * @param y Y-coordinate of the point.
         * @param deps_prime Previously calculated valid spatial step size.
         * @return Spatial step size (delta x, delta y) for the point in the boundary grid, or the previously calculated valid step size.
         */
        double deps_func(double x, double y, double deps_prime);
        
        /**
         * @brief Calculate the geometric tolerance (geps) for a point in the boundary grid.
         *
         * If the point (x, y) is in the boundary grid, this function calculates the geometric tolerance (geps) using the boundary grid information.
         * If the point is not in the boundary grid, it returns the previously calculated valid geometric tolerance (geps_prime).
         *
         * @param x X-coordinate of the point.
         * @param y Y-coordinate of the point.
         * @param geps_prime Previously calculated valid geometric tolerance.
         * @return Geometric tolerance for the point in the boundary grid, or the previously calculated valid tolerance.
         */
        double geps_func(double x, double y, double geps_prime);
        
        /**
         * @brief Calculate the signed distance function (sdf) at a given point (x, y).
         *
         * This function is virtual and defined in the basic_calculation class.
         * It returns the signed distance function (sdf) value at the specified point.
         *
         * @param x X-coordinate of the point.
         * @param y Y-coordinate of the point.
         * @return Signed distance function value at the point (x, y).
         */
        double sdf_func(double x, double y){return sdf_adf(x,y);};

        /**
         * @brief Check if a point is in the outer grid.
         *
         * This function is virtual and defined in the basic_calculation class.
         * It checks if the specified point (x, y) is in the outer grid.
         *
         * @param x X-coordinate of the point.
         * @param y Y-coordinate of the point.
         * @return True if the point is in the outer grid, false otherwise.
         */
        bool pt_in_outer_grid(double x,double y){return size_field->pt_in_outer_grid(x,y);};
        


        // Point initialization variables
        double* xy_id; /**< If undefined: x=ax+bx, y=ay+by; The first Nfixed points are the fixed points. */
        double* xy_id_fixed_pt_init; /**< Temporary variable to store the inputted fixed points */
        int* pt_ctgr; /**< Point category: inner -1, boundary 1 (abs(sdf)<=geps), undefined 0 */
        int inner_pt_ct; /**< Number of inner points */
        double chrtrt_len_h_avg; /**< Average characteristic element size lengthscale h (adaptive) */
        double* chrtrt_len_h; /**< Characteristic element size lengthscale h (adaptive); -1 undefined */
        double* bgrid_geps; /**< Geometric tolerance for boundary grids (adaptive) */
        double* bgrid_deps; /**< Projection delta_x, delta_y for boundary grids (adaptive) */
        adf_2d** bgrid_adf; /**< ADF quad cell for boundary grids (adaptive) */
        adf_stat_2d** bgrid_adf_stat; /**< Statistics for boundary ADFs: error, depth, cell count */
        int Ntotal; /**< Total number of meshing points: fixed + moving points */
        int Ncurrent; /**< Current number of points */
        int Nremain; /**< Remaining number of points */
        int Nfixed; /**< Number of fixed points> */
        bool bdry_adf_construction; /**< Flag indicating whether boundary ADF construction is enabled */


        /**
         * @brief Calculate the normalized density arrays and the characteristic lengthscale h for all inner and boundary grids.
         * @param bgrid_ingeo_ct Number of boundary grids inside the geometry.
         * @param bgrid_outgeo_ct Number of boundary grids outside the geometry.
         * @param bgrid_ingeo_ind Indices of boundary grids inside the geometry.
         * @param bgrid_outgeo_ind Indices of boundary grids outside the geometry.
         * @param density_igrid Normalized density array for inner grids.
         * @param density_bgrid_ingeo Normalized density array for boundary grids inside the geometry.
         * @param Ncurrent_rho_igrid Current density array for inner grids.
         * @param Ncurrent_rho_bgrid_ingeo Current density array for boundary grids inside the geometry.
         */
        void get_normalized_density_chrtrt_h(int bgrid_ingeo_ct,int bgrid_outgeo_ct,
            std::vector<int> &bgrid_ingeo_ind,std::vector<int> &bgrid_outgeo_ind, 
            double *&density_igrid, double *&density_bgrid_ingeo,
            double *&Ncurrent_rho_igrid,double *&Ncurrent_rho_bgrid_ingeo);
        /**
         * @brief Get adaptive arrays for boundary grids.
         * @param scale_adf_error_tol Scale factor for the ADF error tolerance.
         */
        void get_bdry_grid_adaptive_arrays(double scale_adf_error_tol);

        /**
         * @brief Collect valid neighbors for error diffusion in point initialization.
         * @param nei_ij Neighbor grid index.
         * @param nei_ct Neighbor count.
         * @param inner_nei_gi Inner grid neighbor indices.
         * @param ingeo_bdry_nei_gii Boundary grid neighbor indices inside the geometry.
         * @param bgrid_ingeo_outgeo_ctgr Boundary grid category (inside/outside) array.
         */
        void error_diffusion_collect_nei(int nei_ij, int &nei_ct,
            std::vector<int> &inner_nei_gi,std::vector<int> &ingeo_bdry_nei_gii,
            int *&bgrid_ingeo_outgeo_ctgr);

        /**
         * @brief Generate points in a grid.
         * @param generate_pt Flag indicating whether to generate points.
         * @param Npt Number of points to generate.
         * @param i Grid row index.
         * @param j Grid column index.
         * @param ij Grid index.
         * @param gi Geometry index (1: inside, 2: boundary, 0: outside).
         * @param grid_ctgr Grid category.
         * @param pt_counter Point counter.
         * @param xy_id_temp Temporary array for point coordinates.
         * @param pt_ctgr_temp Temporary array for point categories.
         * @param seed Random seed.
         * @param inner_pt_ct_temp Temporary inner point count.
         */
        void generate_pt_in_grid(bool generate_pt,int Npt,
            int i, int j, int ij, int gi, int grid_ctgr,int &pt_counter,
            double *&xy_id_temp,int *&pt_ctgr_temp,
            unsigned int &seed, int &inner_pt_ct_temp);

        /**
         * @brief Perform error diffusion in point initialization.
         * @param generate_pt Flag indicating whether to generate points.
         * @param diff Diffusion factor.
         * @param ij Grid index.
         * @param bgrid_ingeo_outgeo_ctgr Boundary grid category (inside/outside) array.
         * @param Ncurrent_rho_igrid Current density array for inner grids.
         * @param Ncurrent_rho_bgrid_ingeo Current density array for boundary grids inside the geometry.
         */
        void error_diffusion(bool generate_pt, double diff,
            int ij, int *&bgrid_ingeo_outgeo_ctgr,
            double *&Ncurrent_rho_igrid,double *&Ncurrent_rho_bgrid_ingeo);

        /** 
         * Add fixed points into the mesh. Rescale the points using the same scaling as the scaling of shp.
         * 
         * @param Nfixed The number of fixed points
         * @param shp A custom_shape_2d model, whose scaling parameters we will use for normalizing the fixed points
         * @param fixed_pt_list the x and y coordinations of the fixed points
         */
        void add_fixed_points_normailze(int Nfixed_, shape_2d* shp_, double* fixed_pt_list);
        /** 
         * Add fixed points into the mesh.
         * 
         * @param Nfixed The number of fixed points
         * @param fixed_pt_list the x and y coordinations of the fixed points
         */
        void add_fixed_points(int Nfixed_, double* fixed_pt_list);

        /**
         * Adding fixed points to the mesh. Used in pt_init() function.
         */
        void add_fixed_points_procedure();

        /**
         * @brief Initialize points by generating random points.
         * @param Ntotal_ Total number of points.
         */
        void pt_init(int Ntotal_);

        /**
         * @brief Initialize points by reading in a list of points.
         * @param Ntotal_ Total number of points.
         */
        void pt_init(double *pt_list, int Ninit_, int Ntotal_);

        /**
         * @brief Print initial points to a file.
         * @param case_name Name of the case.
         */
        void print_init_pts_to_file(const char *case_name);

        /**
         * @brief Get the ADF value on the boundary grid using ADF quad cells.
         * @param x X-coordinate of the point.
         * @param y Y-coordinate of the point.
         * @return ADF value on the boundary grid.
         */
        double geo_bgrid_adf(double x, double y); //use adf quad cells on the geometry bdry grids, to obtain sdf fast of a point on a bdry grid.
        
        /**
         * @brief Get the signed distance function (SDF) value using ADF if available; otherwise, use shp->sdf.
         * @param x X-coordinate of the point.
         * @param y Y-coordinate of the point.
         * @return SDF value at the point (x, y).
         */
        double sdf_adf(double x, double y); //use bdry_adf if available; otherwise, use shp->sdf
        
        /**
         * @brief Get the characteristic lengthscale h of a point.
         * @param x X-coordinate of the point.
         * @param y Y-coordinate of the point.
         * @return Characteristic lengthscale h at the point (x, y).
         */
        double get_chrtrt_len_h(double x, double y); //chrtrt_len_h of a point (x,y)


        // Meshing variables
        mesh_alg_2d* mesh_alg; /**< Meshing algorithm to choose from */
        int output_interval; /**< Output interval: 0 no output, -1 last final output, -2 initial and final outputs, n>0 -output every n triangulations */


        /**
         * @brief Perform meshing using the specified meshing algorithm.
         *
         * The corresponding meshing class will be created upon choosing the meshing method.
         * The meshing process involves three steps:
         * 1. Initialization: meshing_init()
         * 2. Iteration: meshing_iter()
         * 3. Finalization: meshing_final()
         *
         * @param mesh_alg Pointer to the meshing algorithm to use.
         */
        void meshing(mesh_alg_2d *mesh_alg_){ 
            mesh_alg=mesh_alg_; 
            mesh_alg->meshing_init(); 
            mesh_alg->meshing_iter(); 
            mesh_alg->meshing_final(); 
        } 


        //Timing variables
        double t_generate_pt_in_grid; /**< Time taken for point generation in grids */

        
    };
    


}

#endif
