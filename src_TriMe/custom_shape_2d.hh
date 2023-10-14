#ifndef CUSTOM_SHAPE_2D_HH
#define CUSTOM_SHAPE_2D_HH


#include "basic_calculation.hh"

namespace voro {
    
    /**
     * @brief Class representing a custom 2D shape within a domain container.
     *
     * This class inherits from the `basic_calculation_2d` class and represents a custom 2D shape
     * defined by a set of boundary points. It provides functionality to compute the closest point
     * on the shape, the signed distance to the shape, and other related operations.
     * 
     * It can: 
     * 1. Read in shape boundaries points information.
     *
     * Boundaries are input in the following format:
     * - `std::vector<std::vector<double>> boundaries;`
     * - `std::vector<double> boundary1;`
     * - `std::vector<double> boundary2;`
     * - ...
     * - `boundaries.push_back(boundary1); boundaries.push_back(boundary2);`
     * - ...
     *
     * Each boundary, `boundaries[i]`, represents the ith boundary of the shape. A shape can have multiple shape boundaries.
     * Each boundary is defined as a sequence of points in clockwise order, forming a closed loop with the same end point [x0,y0].
     * For example, `boundary1` can be defined as: `boundary1 = [x0,y0, x1,y1, x2,y2, x3,y3, x0,y0]`, where the interior of the loop is the interior of the shape.
     * 
     * 2. Calculate closest point on shape, closestx, closesty, closestz.
     * 
     * 3. Calculate signed distance to the shape.
     */
    class custom_shape_2d : public basic_calculation_2d{
      public:

        // Domain container grid information
        const double ax; /**< Lower x-coordinate of the container domain. */
        const double bx; /**< Upper x-coordinate of the container domain. */
        const double ay; /**< Lower y-coordinate of the container domain. */
        const double by; /**< Upper y-coordinate of the container domain. */
        const double lx; /**< Length of the container domain in the x-direction (bx - ax). */
        const double ly; /**< Length of the container domain in the y-direction (by - ay). */
        int nx; /**< Number of cells in the x-direction. */
        int ny; /**< Number of cells in the y-direction. */
        double dx; /**< Cell size in the x-direction (lx / nx). */
        double dy; /**< Cell size in the y-direction (ly / ny). */
        double inv_dx; /**< Inverse of the cell size in the x-direction (1 / dx). */
        double inv_dy; /**< Inverse of the cell size in the y-direction (1 / dy). */

        bool normalize_model; /**< If true, the model will be normalized to fit and centered into the container domain.
                               * The model's longest dimension will be in [ax+0.1*(bx-ax), ax+0.9*(bx-ax)]x[ay+0.1*(by-ay),ay+0.9*(by-ay)]. */
        double scale_min_domain_range; /**< scaling parameter, std::min(bx-ax,by-ay). */
        double scale_xmid; /**< scaling parameter, midpoint x-coordinate of the input shape. */
        double scale_ymid; /**< scaling parameter, midpoint y-coordinate of the input shape. */
        double scale_max_range; /**< scaling parameter, max(height,width) of the input shape. */

        int num_t; /**< Number of parallel threads. */


       //shape information
        /**
         * @brief Constructs the `custom_shape_2d` object with the specified shape boundaries, domain parameters, and number of parallel threads.
         *
         * This constructor initializes the `custom_shape_2d` object with the given shape boundaries, domain parameters, and the number
         * of parallel threads. The shape boundaries are specified as a vector of vectors, where each inner vector represents a boundary
         * and contains the x and y coordinates of the boundary points in clockwise order. The `ax`, `bx`, `ay`, and `by` parameters define
         * the domain container's boundaries. If `normalize_model` is set to true, the model will be normalized to fit and centered within
         * the container domain. The `num_t` parameter specifies the number of parallel threads to be used during calculations.
         *
         * @param boundaries_ The vector of shape boundaries.
         * @param ax_ The x-coordinate of the domain container's lower-left corner.
         * @param bx_ The x-coordinate of the domain container's upper-right corner.
         * @param ay_ The y-coordinate of the domain container's lower-left corner.
         * @param by_ The y-coordinate of the domain container's upper-right corner.
         * @param normalize_model_ Flag indicating whether to normalize the model.
         * @param num_t_ The number of parallel threads to be used.
         */
        custom_shape_2d(std::vector<std::vector<double>> boundaries_,double ax_, double bx_, double ay_, double by_, bool normalize_model_, int num_t_);
        
        /**
         * @brief Destroys the `custom_shape_2d` object.
         *
         * This destructor frees the memory allocated for the boundary points and line segment counts.
         */
        ~custom_shape_2d();

        double **b_pts; /**< Boundary points for each boundary. */
        int b_ct;       /**< Number of boundaries. */
        int *seg_ct;    /**< Number of line segments for each boundary. */
        std::vector<std::vector<std::pair<int, int>>> cell_lineSeg; /**< List of boundary line segments in a grid cell. 
                                                                     * Index by grid cell (i, j) of the domain with nx, ny cells. 
                                                                     * Line segments represented by (bi, sj): bi-th b_ct, and sj-th seg_ct[bi]. */


        /**
         * @brief Creates the boundary points for each boundary of the shape.
         *
         * This function takes a vector of boundaries and generates the boundary points
         * for each boundary of the shape. The boundary points are stored in the `b_pts`
         * member variable.
         *
         * @param boundaries_ A vector of vectors representing the boundaries of the shape.
         */
        void create_b_pts(std::vector<std::vector<double>> boundaries_);

        /**
         * @brief Retrieves the list of boundary line segments in each grid cell.
         *
         * This function calculates and stores the list of boundary line segments that intersect
         * each grid cell in the domain. The line segments are represented by pairs of indices:
         * the index of the boundary (`bi`) and the index of the line segment within that boundary (`sj`).
         * The list of line segments for each grid cell is stored in the `cell_lineSeg` member variable.
         */
        void get_cell_lineSeg();

       //shape distance
        /**
         * @brief Calculates the signed distance from a point (x, y) to the shape boundary.
         *
         * This function calculates the signed distance from a point (x, y) to the shape boundary.
         * The signed distance is returned as a double value.
         *
         * @param x The x-coordinate of the point.
         * @param y The y-coordinate of the point.
         * @return The signed distance from the point to the shape boundary.
         */
        double f_b_pts(double x, double y);


        //timing variable
        double t_c1;
        double t_c2;
        double t_c3;
        double t_c4;
        double t_c5;
    };

    
}


#endif