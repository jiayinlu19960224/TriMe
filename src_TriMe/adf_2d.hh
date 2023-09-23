#ifndef ADF_2D_HH
#define ADF_2D_HH


#include "shape_2d.hh"

namespace voro {

    /**
     * @class adf_stat_2d
     * @brief Class representing statistics for an adaptive signed distance field in 2D.
     */
    class adf_stat_2d {
      public:
        /**
         * @brief Default constructor.
         * Initializes the statistics with default values.
         */
        adf_stat_2d():cell_ct(0),depth(0),err(0){};

        /**
         * @brief Destructor.
         * Cleans up any resources associated with the adf_stat_2d object.
         */
        ~adf_stat_2d(){};

        int cell_ct; /**< Number of cells in the quad tree. */
        int depth; /**< Depth of the quad tree. */
        double err; /**< Maximum error of the quad tree. */
    };
    
    //construct adaptive signed distance field based on shape contour line input
    class adf_2d {
      public:
         shape_2d *shp; /**< Pointer to the shape_2d object. */

        /**
         * @brief Constructor with specified parameters.
         * @param adf_stat Pointer to the adf_stat_2d object.
         * @param shp_ Pointer to the shape_2d object.
         * @param x0_ Lower x-coordinate of the quad cell domain.
         * @param x1_ Upper x-coordinate of the quad cell domain.
         * @param y0_ Lower y-coordinate of the quad cell domain.
         * @param y1_ Upper y-coordinate of the quad cell domain.
         * @param v0_ Value of the first corner of the quad cell.
         * @param v1_ Value of the second corner of the quad cell.
         * @param v2_ Value of the third corner of the quad cell.
         * @param v3_ Value of the fourth corner of the quad cell.
         * @param err_tol_ Error tolerance of the subdivision rule.
         * @param max_level_ Maximum level of the quad tree allowed.
         * @param level_ Level of the quad tree.
         */
        adf_2d(adf_stat_2d *adf_stat, shape_2d *shp_, double x0_, double x1_, double y0_, double y1_, double v0_, double v1_, double v2_, double v3_, double err_tol_, int max_level_, int level_);
        
        /**
         * @brief Constructor with specified parameters (without corner values).
         * @param adf_stat Pointer to the adf_stat_2d object.
         * @param shp_ Pointer to the shape_2d object.
         * @param x0_ Lower x-coordinate of the quad cell domain.
         * @param x1_ Upper x-coordinate of the quad cell domain.
         * @param y0_ Lower y-coordinate of the quad cell domain.
         * @param y1_ Upper y-coordinate of the quad cell domain.
         * @param err_tol_ Error tolerance of the subdivision rule.
         * @param max_level_ Maximum level of the quad tree allowed.
         * @param level_ Level of the quad tree.
         */
        adf_2d(adf_stat_2d *adf_stat, shape_2d *shp_, double x0_, double x1_, double y0_, double y1_, double err_tol_, int max_level_, int level_);
        
        /**
         * @brief Destructor.
         * Cleans up any resources associated with the adf_2d object.
         */
        ~adf_2d(){delete ul;delete ur;delete ll;delete lr;};


        double err_tol; /**< Error tolerance of the subdivision rule. */
        int max_level; /**< Maximum level of the quad tree allowed. */
        bool child; /**< Indicates whether the quad has children or not. */
        int level; /**< Level of the quad tree. */
        double x0,x1,y0,y1; /**< Domain of the quad cell. */
        double v0,v1,v2,v3; /**< Four corner values of the quad cell. */
        double v4,v5,v6,v7,v8; /**< The four midpoint values of the edge and the middle point value. */
        double biv4,biv5,biv6,biv7,biv8; /**< The values given by bilinear interpolation of the current cell. */
        adf_2d* ul; /**< Upper left quad child. */
        adf_2d* ur; /**< Upper right quad child. */
        adf_2d* ll; /**< Lower left quad child. */
        adf_2d* lr; /**< Lower right quad child. */

        /**
         * @brief Check if a point is inside the quad cell.
         * @param x x-coordinate of the point.
         * @param y y-coordinate of the point.
         * @return true if the point is inside the quad cell, false otherwise.
         */
        bool point_in_cell(double x, double y);

        /**
         * @brief Locate the cell in the quad tree and get the point value with bilinear interpolation of the smallest cell.
         * @param x x-coordinate of the point.
         * @param y y-coordinate of the point.
         * @return The interpolated point value.
         */
        double getVal(double x, double y); 

        /**
         * @brief Calculate the transformed epsilon value for a given x-coordinate. Used for bilinear interpolation.
         * @param x x-coordinate.
         * @return The epsilon value.
         */
        double epsilon(double x);

        /**
         * @brief Calculate the transformed eta value for a given y-coordinate. Used for bilinear interpolation.
         * @param y y-coordinate.
         * @return The eta value.
         */
        double eta(double y);

        /**
         * @brief Perform bilinear interpolation to get the point value in the cell.
         * @param x x-coordinate of the point.
         * @param y y-coordinate of the point.
         * @return The interpolated point value.
         */
        double bilinear_i(double x, double y); 

        /**
         * @brief Determine whether to subdivide the quad cell or not by comparing v4, v5, v6, v7, v8 with their bilinear interpolated values
         * @return true if subdivision is needed, false otherwise.
         */
        bool subdivide(); 

        /**
         * @brief Calculate the bilinear interpolated values biv4, biv5, biv6, biv7, biv8 for the cell.
         */
        void get_biv45678();  

        /**
         * @brief Calculate the shape signed distance values v4, v5, v6, v7, v8 for the cell.
         */
        void get_v45678(); 

        /**
         * @brief Print the quad tree data to a file (used in print_tree() function).
         * @param[out] outFile Reference to the output file stream.
         */
        void print_tree_func(FILE *&outFile); //used in print_tree(fp)

        /**
         * @brief Print the quad tree data to a text file.
         * @param[in] fp Path to the output text file.
         */
        void print_tree(const char *fp); 

        /**
         * @brief Calculate the maximum error of the quad tree.
         * @param adf_stat Pointer to the adf_stat_2d object.
         */
        void get_err(adf_stat_2d *adf_stat); 
    };
}


#endif



