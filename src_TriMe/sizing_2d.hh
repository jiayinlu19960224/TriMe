#ifndef SIZING_2D_HH
#define SIZING_2D_HH

#include "voro++.hh"
#include "basic_calculation.hh"
#include "shape_2d.hh"
#include "config.hh"

namespace voro {

    /**
   * @brief The sizing_2d class represents a element sizing calculation of a 2D shape.
   * 
   * This class is derived from the basic_calculation_2d class and provides methods for performing sizing calculations
   * based on the underlying geometry grid.
   */
    class sizing_2d : public basic_calculation_2d {
      public:
        /**
       * @brief Constructs a sizing_2d object with the given shape_2d.
       * 
       * @param shp_ Pointer to a shape_2d object providing the underlying geometry grid and sdf function.
       */
        sizing_2d(shape_2d *shp_);

        ~sizing_2d(){}; /**< Destructor for the sizing_2d class. */

        int num_t; /**< The number of threads for parallel execution. */

        // A. get the underlying geometry grid info
        shape_2d* shp; /**< Pointer to a shape_2d object that provides the sdf(x,y) function and underlying geometry grid. */

        double ax, bx, ay, by; /**< The x and y coordinates of the domain boundaries. */
        int gnx, gny, gnxy; /**< The number of grid cells for the geometry grid. */
        double gdx, gdy, diag_gdxy, inv_gdx, inv_gdy; /**< Grid spacing and related values. */
        int geo_bgrid_ct; /**< The number of geometry boundary grids. */
        int geo_igrid_ct; /**< The number of geometry inner grids. */
        int geo_ogrid_ct; /**< The number of geometry outer grids. */

         /**
         * @brief Calculates and returns the signed distance field (SDF) value for the given coordinates (x, y).
         * 
         * @param x The x-coordinate.
         * @param y The y-coordinate.
         * @return The signed distance field value.
         */
        double sdf(double x, double y){return shp->sdf(x,y);};

        /**
         * @brief Convenience function that calls sdf(x, y).
         * 
         * This function is provided for compatibility purposes.
         * 
         * @param x The x-coordinate.
         * @param y The y-coordinate.
         * @return The signed distance field value.
         */
        double sdf_func(double x, double y){return sdf(x,y);}

        /**
         * @brief Checks if the given point (x, y) is inside the outer grid.
         * 
         * @param x The x-coordinate of the point.
         * @param y The y-coordinate of the point.
         * @return True if the point is inside the outer grid, false otherwise.
         */
        bool pt_in_outer_grid(double x,double y);

        /**
         * @brief Returns the geometry grid value for the given cell index ij.
         * 
         * @param ij The cell index.
         * @return The geometry grid value.
         */
        int geo_grid(int ij){return shp->geo_grid[ij];}; 

        /**
         * @brief Returns the ij value of the boundary grid at the given index.
         * 
         * @param ind The index of the boundary grid.
         * @return The ij value of the boundary grid.
         */
        int geo_bgrid_ij(int ind){return shp->geo_bgrid_ij[ind];}; 

        /**
         * @brief Returns the ij value of the inner grid at the given index.
         * 
         * @param ind The index of the inner grid.
         * @return The ij value of the inner grid.
         */
        int geo_igrid_ij(int ind){return shp->geo_igrid_ij[ind];}; 

        /**
         * @brief Returns the ij value of the outer grid at the given index.
         * 
         * @param ind The index of the outer grid.
         * @return The ij value of the outer grid.
         */
        int geo_ogrid_ij(int ind){return shp->geo_ogrid_ij[ind];}; //outer grid ij's array


        /**
         * @brief Calculates and returns the sizing value for the given coordinates (x, y).
         * 
         * This method should be overridden in derived classes to provide the actual sizing calculation.
         * 
         * @param x The x-coordinate.
         * @param y The y-coordinate.
         * @return The sizing value.
         */
        virtual double getSizingVal(double x, double y){return 0;};  
        
        /**
         * @brief Calculates and returns the density value for the given coordinates (x, y).
         * 
         * A default calculation is implemented here. 
         * This method should be overridden in the derived sizing_2d_automatic classes to provide the actual density calculation in that class.
         * 
         * @param x The x-coordinate.
         * @param y The y-coordinate.
         * @return The density value.
         */
        virtual double getDensityVal(double x, double y){
          double sizing_val=getSizingVal(x,y);
          if(density_sizing_exp==2){
            return 1.0/(sizing_val*sizing_val);
          }
          else{
            return 1.0/pow(sizing_val,density_sizing_exp);  //rho=1/(h**2)
          }
          
        };
        

        /**
         * @brief Only called in CVD meshing for automatic sizing field for outer grids.
         */
        virtual void get_ogrid_sizing_density_grid(){}; 
        

        /**
         * @brief Prints the sizing and density fields to files.
         * 
         * @param case_name The name of the case or directory.
         */
        void print_fields_to_file(const char *case_name); //sizing, density

    };


    /**
   * @brief Class for automatic element sizing and density computation in 2D.
   *
   * This class extends the sizing_2d class and adds functionality for automatic sizing and density field computation.
   */
    class sizing_2d_automatic : public sizing_2d {
      public:
        /**
         * @brief Constructor for sizing_2d_automatic class.
         * @param shp_ Pointer to the shape_2d object.
         * @param K_ Degradation control parameter.
         */
        sizing_2d_automatic(shape_2d *shp_, double K_) : sizing_2d(shp_), K(K_),

        igrid_sizing(new double[geo_igrid_ct]),igrid_density(new double[geo_igrid_ct]),
        bgrid_sizing(new double[geo_bgrid_ct]),bgrid_density(new double[geo_bgrid_ct]),
        ogrid_sizing(new double[1]),ogrid_density(new double[1]),ogrid_sizing_density(false),

        bdry_pt_ct(0),medial_pt_ct(0),
        bdry_pt_xy(new double[1]),bdry_pt_lfs(new double[1]),medial_pt_xy(new double[1])
        {
          get_density_grid();
        }

        /**
         * @brief Destructor for sizing_2d_automatic class.
         */
        ~sizing_2d_automatic();

        

        //automatic sizing and density field computation
        double K;                          /**< Degradation control parameter. 0 is uniform. */
        double* igrid_sizing;              /**< Sizing field for geo_inner_grids. */
        double* igrid_density;             /**< Density field for geo_inner_grids. */
        double* bgrid_sizing;              /**< Sizing field for geo_bdry_grids. */
        double* bgrid_density;             /**< Density field for geo_bdry_grids. */
        double* ogrid_sizing;              /**< Sizing field for geo_outer_grids. Only created for CVD meshing method. */
        double* ogrid_density;             /**< Density field for geo_outer_grids. Only created for CVD meshing method. */
        bool ogrid_sizing_density;         /**< True if CVD meshing is used. */

        int bdry_pt_ct;                    /**< Count of boundary points. */
        int medial_pt_ct;                  /**< Count of medial axis points. */
        double* bdry_pt_xy;                /**< Array of boundary point coordinates (size: 2 * bdry_pt_ct). */
        double* bdry_pt_lfs;               /**< Array of boundary point local feature sizes (size: bdry_pt_ct). */
        double* medial_pt_xy;              /**< Array of medial axis point coordinates. */
        double deps_temp;                  /**< Temporary variable. Spatial stepsize to use in finite difference in point projection. */
        double geps_temp;                  /**< Temporary variable: Geometry boundary tolerance. */

        /**
         * @brief Computes the value of deps (spatial stepsize to use in finite difference in point projection) for a given point (x, y).
         *
         * This function returns the adaptive value of deps if it is available (point lies in boundary grid). 
         * Otherwise it returns the default value deps_prime.
         *
         * @param x The x-coordinate of the point.
         * @param y The y-coordinate of the point.
         * @param deps_prime The default value of deps if adaptive value is not available.
         * @return The computed value of deps for the point (x, y).
         */
        double deps_func(double x, double y, double deps_prime){return deps_temp;};

        /**
         * @brief Computes the value of geps (geometry boundary tolerance) for a given point (x, y).
         *
         * This function returns the adaptive value of geps if it is available (i.e., the point lies in a boundary grid).
         * Otherwise, it returns the default value geps_prime.
         *
         * @param x The x-coordinate of the point.
         * @param y The y-coordinate of the point.
         * @param geps_prime The default value of geps if the adaptive value is not available.
         * @return The computed value of geps for the point (x, y).
         */
        double geps_func(double x, double y, double geps_prime){return geps_temp;};

        /**
         * @brief Function to compute boundary points.
         */
        void get_bdry_pts();

        /**
         * @brief Function to compute medial axis points.
         */
        void get_medial_pts();

        /**
         * @brief Function to compute boundary point local feature size values.
         */
        void get_bdry_pt_lfs();

        /**
         * @brief Function to compute the sizing grid.
         */
        void get_sizing_grid();

        /**
         * @brief Function to compute the density grid.
         */
        void get_density_grid();  

        /**
         * @brief Function to compute the sizing and density grid for the outer grids.
         * Only called in CVD meshing.
         */
        void get_ogrid_sizing_density_grid(); 

        /**
         * @brief Function to get the sizing value at a point by reading the corresponding sizing grid value.
         * @param x X-coordinate of the point.
         * @param y Y-coordinate of the point.
         * @return The sizing value at the given point.
         */
        double getSizingVal(double x, double y);  

        /**
         * @brief Function to get the density value at a point by reading the corresponding density grid value.
         * @param x X-coordinate of the point.
         * @param y Y-coordinate of the point.
         * @return The density value at the given point.
         */
        double getDensityVal(double x, double y);   

        /**
         * @brief Function to print medial points and boundary points to a file.
         * @param case_name The name of the output file.
         */
        void print_pts_to_file(const char *case_name); 
    };


}


#endif

