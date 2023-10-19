#ifndef MESH_ALG_2D_HH
#define MESH_ALG_2D_HH

#include "voro++.hh"
#include "basic_calculation.hh"
#include <unordered_map>

namespace voro {

    class parallel_meshing_2d;
    struct HE_edge;
    struct HE_vert;
    struct HE_face;

    /**
     * @breif Edge struct for constructing half-edge (HE) data structure.
     */
    struct HE_edge
    { 
      HE_edge(): eu(-1),ev(-1),prev(0),next(0),pair(0),vert(0),face(0){};
      int eu;
      int ev;
      HE_edge* prev; //previous half-edge around the triangle face
      HE_edge* next; //next HE around the triangle face
      HE_edge* pair; //oppositely oriented adjacent HE
      HE_vert* vert; //vertex at end of HE
      HE_face* face; //triangle face the HE borders
    };

    /**
     * @breif Vertex struct for constructing half-edge (HE) data structure.
     */
    struct HE_vert
    {
      HE_vert(): HE_vert_id(-1),edge(0){};
      int HE_vert_id; 
      HE_edge* edge; //One of the HE emanating from the vertex
    };

    /**
     * @breif Triangle face struct for constructing half-edge (HE) data structure.
     */
    struct HE_face
    {
      HE_face():HE_face_id(-1),edge(0){};
      int HE_face_id;
      HE_edge* edge; //One of the HE bordering the triangle face.
    };

    /**
     * @brief Wall object for cutting a Voronoi cell to an octagon shape based on the local characteristic element edge length.
     * 
     * This class inherits from the base wall_2d class and provides functionality to cut a cell into an octagon shape.
     * The size of the octagon depends on the local characteristic element edge length.
     */
    class wall_is_2d : public wall_2d {
     public:
            /**
             * @brief Constructor for the wall_is_2d class.
             * @param pm2d_ Pointer to the parallel_meshing_2d object.
             */
             wall_is_2d(parallel_meshing_2d *pm2d_);

             /**
             * @brief Destructor for the wall_is_2d class.
             */
             ~wall_is_2d();

             /**
             * @brief Update the triangle length criteria based on the local characteristic length.
             */
             void update_tria_length_cri();

            /**
             * @brief Update the triangle length criteria for a specific geometry grid.
             * @param ij The index of the geometry grid.
             * @param fac The factor by which to update the triangle length criteria.
             */
             void update_tria_length_cri(int ij, double fac);

            /**
             * @brief Check if a point is inside the walls.
             * @param x The x-coordinate of the point.
             * @param y The y-coordinate of the point.
             * @return True if the point is inside the walls, false otherwise.
             */
             bool point_inside(double x,double y) {return true;};

             /**
               * Base function to cut a Voronoi cell of a point (x,y).
               *
               * @tparam v_cell_2d The type of voronoicell class.
               * @param c The Voronoicell object to be cut.
               * @param x The x-coordinate of the point.
               * @param y The y-coordinate of the point.
               * @return True if the cell is successfully cut, false otherwise.
               */
             template<class v_cell_2d>
             bool cut_cell_base(v_cell_2d &c,double x,double y);
             
            /**
             * Cuts a voronoicell_2d object at the given point (x, y) with octagon walls.
             *
             * @param c The voronoicell_2d object to be cut.
             * @param x The x-coordinate of the point.
             * @param y The y-coordinate of the point.
             * @return True if the cell is successfully cut, false otherwise.
             */
             bool cut_cell(voronoicell_2d &c,double x,double y) {return cut_cell_base(c,x,y);};
             
             /**
               * Cuts a voronoicell_neighbor_2d object for the given point (x, y) with octagon walls.
               *
               * @param c The voronoicell_neighbor_2d object to be cut.
               * @param x The x-coordinate of the point.
               * @param y The y-coordinate of the point.
               * @return True if the cell is successfully cut, false otherwise.
               */
             bool cut_cell(voronoicell_neighbor_2d &c,double x,double y) {return cut_cell_base(c,x,y);};
             
             voronoicell_2d v; /**< The voronoicell_2d object. */
             parallel_meshing_2d *pm2d; /**< Pointer to the parallel_meshing_2d object. */
             double *tria_length_cri;  /**< Pointer to the array of lengths used to initialize cells to be squares with length tria_length_cri. */
             double cut_cell_base_octE_fac; /**< The octagon factor used in the cut_cell_base function. */
             int gnxy,gnx,gny; /**< Variables related to the geometry grid. */
             double gdx,gdy,inv_gdx, inv_gdy; /**< Variables related to the spatial step size and its inverse. */
             double ax,ay,bx,by; /**< Variables representing the full domain size. */
    };
    
    /**
     * @brief Construct adaptive signed distance field based on shape contour line input.
     * 
     * This class inherits from the basic_calculation_2d class and provides functionality for meshing algorithms in 2D.
     */
    class mesh_alg_2d : public basic_calculation_2d {
      public:

        /**
         * @brief Constructor for the mesh_alg_2d class.
         * @param pm2d_ Pointer to the parallel_meshing_2d object.
         */
        mesh_alg_2d(parallel_meshing_2d *pm2d_);

        /**
         * @brief Destructor for the mesh_alg_2d class.
         */
        ~mesh_alg_2d();

        parallel_meshing_2d *pm2d; /**< Pointer to the parallel_meshing_2d object. */
        wall_is_2d wis; /**< Object of the class wall_is_2d, used for bounding Voronoi cells with octagon walls. */
        int num_t; /**< Number of threads to use in parallel meshing. */

        /**
         * @brief Change the number of threads used in parallel meshing.
         * @param new_num_t The new number of threads.
         */
        void change_number_thread(int new_num_t);

        double ax, bx, ay, by, lx, ly; /**< Meshing domain size variables. */
        int gnx, gny, gnxy; /**< Meshing geometry grid variables. */
        double gdx, gdy, diag_gdxy, inv_gdx, inv_gdy; /**< Meshing geometry grid spacing variables. */
        int geo_bgrid_ct; /**< Number of geometry boundary grids. */
        int geo_igrid_ct; /**< Number of geometry inner grids. */
        //pt info
        int Ntotal; /**< Total number of meshing points. */
        int Ncurrent; /**< Current number of meshing points. */
        int Nremain; /**< Remaining number of meshing points, Nremain=Ntotal-Ncurrent. */
        int Nfixed;

        //helper functions and variables
        /**
         * @brief Get the category of a point.
         * @param id The ID of the point.
         * @return The category of the point: inner -1, boundary 1, undefined 0.
         */
        int pt_ctgr(int id); 

        /**
         * @brief Get the characteristic edge length at a point.
         * @param x The x-coordinate of the point.
         * @param y The y-coordinate of the point.
         * @return The characteristic edge length at the point.
         */
        double get_chrtrt_len_h(double x, double y);

        /**
         * @brief Get the characteristic edge length of a grid cell.
         * @param ij The grid cell index.
         * @return The characteristic edge length of the grid cell.
         */
        double chrtrt_len_h(int ij);


        /**
         * Get the geometry boundary tolerance (geps) for a given point (x, y) lie in boundary grid.
         * @param x The x-coordinate of the point.
         * @param y The y-coordinate of the point.
         * @return The geometry boundary tolerance (geps) for the point.
         */
        double get_bgrid_geps(double x, double y);
        /**
         * Get the spatial step size (deps) for finite difference calculation in the point projection step,
         * for a given point (x, y) lie in boundatry grid.
         * @param x The x-coordinate of the point.
         * @param y The y-coordinate of the point.
         * @return The spatial step size (deps) for the point.
         * @throws std::exception if the point is not in the geometry boundary grid.
         */
        double get_bgrid_deps(double x, double y);

        /**
         * @brief Returns the category of the grid cell at the given index in the underlying geometry grid.
         * @param ij The index of the grid cell.
         * @return The status of the grid cell: -1, -2, -3 for inside; 1, 2, 3 for boundary; gnxy+1, gnxy+2, ... for outside.
         */
        int geo_grid(int ij);

        /**
         * @brief Returns the category of the underlying geometry grid cell that the point (x, y) lies in.
         * @param x The x-coordinate of the point.
         * @param y The y-coordinate of the point.
         * @return The status of the grid cell: -1, -2, -3 for inside; 1, 2, 3 for boundary; gnxy+1, gnxy+2, ... for outside.
         */
        int geo_grid(double x, double y);

        /**
         * @brief Get the current x-coordinate of a point.
         * @param id The ID of the point.
         * @return The current x-coordinate of the point.
         */
        double current_x(int id); 

        /**
         * @brief Get the current y-coordinate of a point.
         * @param id The ID of the point.
         * @return The current y-coordinate of the point.
         */
        double current_y(int id); 

        /**
         * @brief Update the x-coordinate of a point with a new position.
         * @param id The ID of the point.
         * @param new_x The new x-coordinate of the point.
         */
        void update_x(int id, double new_x); 

        /**
         * @brief Update the y-coordinate of a point with a new position.
         * @param id The ID of the point.
         * @param new_y The new y-coordinate of the point.
         */
        void update_y(int id, double new_y); 

        int inner_pt_ct_temp; /**< Temporary variable to store the count of inner points. */

        //Thread private variables to collect voro and tria info needed 
        int *tria_ct_private; /**< Thread-private variable to collect the number of triangles. */
        std::vector<int> *tria_vertex_private; /**< Thread-private variable to collect the vertex indices of triangles. */
        std::vector<double> *tria_centroid_private; /**< Thread-private variable to collect the centroids of triangles. */
        std::vector<double> *tria_ccrd_h_private; /**< Thread-private variable to collect the circumradius-to-characteristic-length ratios of triangles. */

        //triangl quality: aspect ratio, edge ratio
        std::vector<double> *tria_ar_private; /**< Thread-private variable to collect the aspect ratios of triangles. */
        std::vector<double> *tria_er_private; /**< Thread-private variable to collect the edge ratios of triangles. */
        double *max_tria_ar_private; /**< Thread-private variable to store the maximum triangle aspect ratio. */

        //triangle bars/edges
        std::vector<int> *barid_private; /**< Thread-private variable to collect the IDs of unique/non-duplicate bars. */
        int *bar_ct_private; /**< Thread-private variable to store the number of unique/non-duplicate bars. */

        //triangle quality 
        double *alpha_mean_arsum_private; /**< Thread-private variable to store the sum of aspect ratios for calculating the alpha-mean. */
        double *alpha_mean_ersum_private; /**< Thread-private variable to store the sum of edge ratios for calculating the alpha-mean. */

        //meshing framework control variables
        double *xy_id_new; /**< Array to store the updated point positions. */
        int output_interval; /**< Interval for outputting results. 0: No output. -1: Last final output. 10: Output every 10 triangulations. */
        bool Continue; /**< Flag indicating whether to continue the meshing process. */
        bool reTria; /**< Flag indicating whether to re-triangulate. */
        bool addPt; /**< Flag indicating whether to add new points. */
        int all_iter_ct;  /**< Total number of iterations. */
        int tria_iter_ct; /**< Number of re-triangulations. */
        int inbtw_tria_iter_ct; /**< Number of iterations since the previous re-triangulation. */
        double *stop_Continue_mvmt_thres; /**< Array to store the meshing stop/continue inner point movement thresholds for each geometry grid. Dummy -1. If all inner points have movement <thres, Continue=false. */
        int inner_pt_over_stop_Continue_mvmt_thres_ct; /**< Counter for the number of inner points that have movement over the meshing stop/continue inner point movement threshold. */
        double *pt_mvmt_dis_thres; /**< Array to store the point movement distance thresholds for each geometry grid. Dummy -1. */

        //triangle related variables
        int tria_ct; /**< Number of triangles in the mesh. */
        int *tria_vertex; /**< Array storing the vertex ID's of mesh triangles. */
        int tria_vertex_tct; /**< Current tria_vertex array size. */
        double *tria_ccrd_h; /**< Array storing triangles' circumradius/chrtrt_h_i ratios. */
        int tria_ccrd_h_tct; /**< Current tria_ccrd_h array size. */
        double *tria_centroid; /**< Array storing triangles' centroids. */
        int tria_centroid_tct; /**< Current tria_centroid array size. */
        //triangle quality measure: edge ratio, aspect ratio, and their max, mean, standard deviation
        double *tria_ar; /**< Array storing triangles' aspect ratios. */
        int tria_ar_tct; /**< Current tria_ar array size. */
        double *tria_er; /**< Array storing triangles' edge ratios. */
        int tria_er_tct; /**< Current tria_er array size. */
        //triangle edge info
        int bar_ct; /**< Number of unique/non-duplicate bars. */
        int *barid;  /**< Array storing the IDs of unique/non-duplicate bars. */
        int barid_barct; /**< Current size of the barid array. */
        //Half-edge(HE) data structure to represent the triangulation: variables
        bool HE_exist;
        //Cantor pairing function for hash function
        inline size_t KKey(int i,int j) {
          return (i+j)*(i+j+1)/2+i;
        }
        //Cantor pairing function for hash function
        inline size_t KKey(std::pair<unsigned int, unsigned int> uv) {
          int i=uv.first; int j=uv.second;
          return KKey(i,j);
        }
        std::unordered_map<size_t, HE_edge*> Edges; /**< The half-edges of the triangulation. */
        HE_face* Faces; /**< The triangle faces in the triangulation. */
        HE_vert* Vertices; /**< The vertices in the triangulation. */
        std::unordered_map<unsigned int, std::pair<unsigned int, unsigned int>> Bdry_Edges; /**< The boundary edges. Key: u; Value: <v,0/1>. 0/1 means have been connected with another bdry edge or not. */
        std::vector<std::pair<unsigned int, unsigned int>> Bdry_Edges_Start; /**< The starting HE of the boundaries. Each boundary has one. */
        /**
         * @brief Construct the Half-edge data structure for the triangulation.
         *      The boundary HE are those with NULL face pointer.
         *      Requires tria_vertex to have been stored in the voro-tria calculation.
         */
        void construct_HE();
        void print_bdry_CCW();

        //add_pt and termination criteria
        double alpha_mean_tria_ar; /**< Alpha-mean value of triangles' aspect ratio */
        double alpha_mean_tria_er; /**< Alpha-mean value of triangles' edge ratio */
        double previous_alpha_mean_tria_ar; /**< Alpha-mean value of triangles' aspect ratio in previous iteration */
        double previous_alpha_mean_tria_er; /**< Alpha-mean value of triangles' edge ratio in previous iteration */

        //other mesh quality statistics
        double max_tria_ar_prev; /**< Maximum triangle aspect ratio in the previous iteration */
        double max_tria_ar; /**< Maximum triangle aspect ratio */
        double max_tria_er; /**< Maximum triangle edge ratio */
        double mean_tria_ar; /**< Mean triangle aspect ratio */
        double mean_tria_er; /**< Mean triangle edge ratio */
        double median_tria_ar; /**< Median triangle aspect ratio */
        double median_tria_er; /**< Median triangle edge ratio */
        double std_tria_ar; /**< Standard deviation of triangle aspect ratios */
        double std_tria_er; /**< Standard deviation of triangle edge ratios */
        
        //meshing framework sub-functions
        /**
         * Updates the stop/continue pt movement threshold for the given geometry grid index.
         *
         * @param ij The index of the geometry grid.
         */
        void update_stop_Continue_mvmt_thres(int ij);
        
        /**
         * Updates the point movement distance threshold for the given geometry grid index.
         *
         * @param ij The index of the geometry grid.
         */
        void update_pt_mvmt_dis_thres(int ij);
        
        /**
         * Checks if the movement of the inner point (px_old, py_old) to (px_new, py_new) exceeds the stop/continue movement threshold.
         *
         * @param pi The index of the inner point.
         * @param px_old The old x-coordinate of the point.
         * @param py_old The old y-coordinate of the point.
         * @param px_new The new x-coordinate of the point.
         * @param py_new The new y-coordinate of the point.
         * @return True if the movement exceeds the threshold, false otherwise.
         */
        bool inner_pt_movement_over_stop_Continue_mvmt_thres(int pi, double px_old, double py_old, double px_new, double py_new);
        
        /**
         * Checks if the movement of the inner point with the given index exceeds the stop/continue movement threshold.
         *
         * @param pi The index of the inner point.
         * @param pt_movement The movement distance of the point.
         * @param stop_Continue_mvmt_thres_ij The stop/continue movement threshold for the corresponding geometry grid.
         * @return True if the movement exceeds the threshold, false otherwise.
         */
        bool inner_pt_movement_over_stop_Continue_mvmt_thres(int pi,double pt_movement,double stop_Continue_mvmt_thres_ij);
        
        /**
         * Determines whether to add new points and retriangulate based on triangle quality.
         */
        void determine_addPt_and_reTria_tria_quality();

        /**
         * Determines whether to add new points and retriangulate based on point movement.
         */
        void determine_addPt_and_reTria_movement();

        /**
         * Determines whether to continue meshing based on point movement.
         */
        void determine_Continue_pt_movement();

        /**
         * Determines whether to continue meshing based on triangle qualities.
         */
        void determine_Continue_tria_quality();

        /**
         * Adds new points to the mesh.
         */
        void add_new_pts();

        /**
         * Checks the Voronoi neighbors' particle requirements and recomputes if necessary.
         *
         * The requirement is for the particle to have two consecutive neighbors that form a triangle. 
         * 
         * @param recompute Flag indicating whether recomputation is required.
         * @param neigh The vector of neighbor indices.
         * @param x The x-coordinate of the point.
         * @param y The y-coordinate of the point.
         */
        void check_voro_nei_particle_requirement(bool &recompute,std::vector<int> neigh,double x,double y);
        
        /**
         * Computes and stores Voronoi information, such as triangle bars, triangle quality, triangle vertices, centroid, etc.
         *
         * @param store_tria_bar Flag indicating whether to store triangle bar information.
         * @param store_triangle_quality Flag indicating whether to store triangle quality information.
         * @param store_tria_vertex Flag indicating whether to store triangle vertex information.
         * @param store_tria_centroid Flag indicating whether to store triangle centroid information.
         * @param store_ccrd_h_ratio Flag indicating whether to store the ratio of circumcircle radius to charasteristic edge length information.
         * @param store_max_tria_ar Flag indicating whether to store the maximum triangle aspect ratio information.
         * @param cvd_pt_update Flag indicating whether to update the points during CVD.
         * @param store_xy_id Flag indicating whether to store xy_id information.
         * @param xy_id_store_array The array to store xy_id information.
         */
        void voro_compute_and_store_info(bool store_tria_bar,
          bool store_triangle_quality,
          bool store_tria_vertex,
          bool store_tria_centroid,
          bool store_ccrd_h_ratio,bool store_max_tria_ar,
          bool cvd_pt_update,
          bool store_xy_id, double *&xy_id_store_array);


        /**
         * Checks if a triangle defined by three vertex IDs is valid based on the geometry shape.
         *
         * @param x0 The x-coordinate of the first vertex.
         * @param y0 The y-coordinate of the first vertex.
         * @param x1 The x-coordinate of the second vertex.
         * @param y1 The y-coordinate of the second vertex.
         * @param x2 The x-coordinate of the third vertex.
         * @param y2 The y-coordinate of the third vertex.
         * @return True if the triangle is valid, false otherwise.
         */
        bool valid_tria(double x0,double y0,double x1,double y1,double x2,double y2);

        /**
         * @brief Extracts triangle information for a triangle defined by vertex ids: vid1, vid2, vid3.
         *
         * @param[in] vid1 First vertex id.
         * @param[in] vid2 Second vertex id.
         * @param[in] vid3 Third vertex id.
         * @param[in] store_tria_vertex Flag to indicate if triangle vertex information should be stored.
         * @param[in] do_tria_ct Flag to indicate if triangle count should be incremented.
         * @param[in] store_tria_centroid Flag to indicate if triangle centroid information should be stored.
         * @param[in] store_ccrd_h_ratio Flag to indicate if circumradius-to-characteristic-edge ratio should be stored.
         * @param[in] store_triangle_quality Flag to indicate if triangle quality (aspect ratio and edge ratio) should be stored.
         * @param[in] store_max_tria_ar Flag to indicate if maximum triangle aspect ratio should be stored.
         * @param[out] tria_vertex_p Vector to store triangle vertex ids.
         * @param[out] tria_centroid_p Vector to store triangle centroid coordinates.
         * @param[out] tria_ccrd_h_p Vector to store circumradius-to-characteristic-edge ratios.
         * @param[out] tria_ar_p Vector to store triangle aspect ratios.
         * @param[out] tria_er_p Vector to store triangle edge ratios.
         * @param[in,out] tria_ct_p Triangle count.
         * @param[in,out] alpha_mean_arsum_p Accumulated sum of square roots of aspect ratios.
         * @param[in,out] alpha_mean_ersum_p Accumulated sum of square roots of edge ratios.
         * @param[in] tria_vertex_p_tct Size of tria_vertex_p vector for thread-safe access.
         * @param[in] tria_centroid_p_tct Size of tria_centroid_p vector for thread-safe access.
         * @param[in] tria_ccrd_h_p_tct Size of tria_ccrd_h_p vector for thread-safe access.
         * @param[in] tria_ar_p_tct Size of tria_ar_p vector for thread-safe access.
         * @param[in] tria_er_p_tct Size of tria_er_p vector for thread-safe access.
         * @param[in,out] max_tria_ar_ Maximum triangle aspect ratio.
         * @param[in] cvd_pt_update If using CVD meshing or not
         * @param[in,out] use_random_point If using CVD meshing: Decide whether to move point to the centroid of the needle type triangle it connects to instead, to eliminate the needle triangle.
         * @param[in,out] random_point_x The x-coordinate of the needle type triangle centroid the point should move to.
         * @param[in,out] random_point_y The y-coordinate of the needle type triangle centroid the point should move to.
         *
         * @return True if the triangle is valid, false otherwise.
         */
        bool extract_tria_info(int vid1,int vid2, int vid3,
    
            bool store_tria_vertex, bool do_tria_ct, 
            bool store_tria_centroid,bool store_ccrd_h_ratio,
            bool store_triangle_quality,bool store_max_tria_ar,

            std::vector<int> &tria_vertex_p,
            std::vector<double> &tria_centroid_p,std::vector<double> &tria_ccrd_h_p,
            std::vector<double> &tria_ar_p,std::vector<double> &tria_er_p,

            int &tria_ct_p, 

            double &alpha_mean_arsum_p, double &alpha_mean_ersum_p,


            int tria_vertex_p_tct,
            int tria_centroid_p_tct,
            int tria_ccrd_h_p_tct,
            int tria_ar_p_tct,
            int tria_er_p_tct, 

            double &max_tria_ar_,

            bool cvd_pt_update, bool &use_random_point, double &random_point_x, double &random_point_y, unsigned int &seed_i);

        /**
         * Perform sphere tracing from an inner point to an outer point lying on the boundary grid.
         * This function finds a boundary point, which is the intersection point of the inner point-outer point line segment with the geometry boundary.
         *
         * @param[out] final_x The final x-coordinate of the boundary point.
         * @param[out] final_y The final y-coordinate of the boundary point.
         * @param[in] px_old The x-coordinate of the inner point.
         * @param[in] py_old The y-coordinate of the inner point.
         * @param[in] px_new The x-coordinate of the outer point lying on the boundary grid.
         * @param[in] py_new The y-coordinate of the outer point lying on the boundary grid.
         */
        void do_sphere_tracing(double &final_x,double &final_y,double px_old,double py_old,double px_new,double py_new);
        
        /**
         * Perform point movement treatment and projection.
         * This function handles the movement of a point and updates its category (inner or boundary).
         * It checks the movement distance, the grid location of the new point, and performs necessary adjustments.
         *
         * @param[in] pid The ID of the point.
         * @param[in] px_old The x-coordinate of the old point.
         * @param[in] py_old The y-coordinate of the old point.
         * @param[in] px_new The x-coordinate of the new point.
         * @param[in] py_new The y-coordinate of the new point.
         * @param[out] final_x The final x-coordinate of the point after treatment and projection.
         * @param[out] final_y The final y-coordinate of the point after treatment and projection.
         */
        void pt_mvmt_treatment_projection(int pid,double px_old,double py_old,double px_new,double py_new,double &final_x, double &final_y);
        
        //output related functions
        bool printOutputs; /**< Flag indicating whether to print outputs. */

        /**
         * Determines whether to print outputs based on the value of printOutputs.
         */
        void determine_printOutputs();

        /**
         * Prints the particle coordinates to a file with the given file name prefix.
         *
         * @param file_name_prefix The prefix of the file name to be used for printing.
         */
        void print_particle_coords(const char *file_name_prefix);
        
        /**
         * Prints the xy_id information to a file with the given file name prefix.
         *
         * @param file_name_prefix The prefix of the file name to be used for printing.
         */
        void print_xy_id(const char *file_name_prefix); 

        /**
         * Boolean variable tracking whether triangle id's are ordered in CCW direction.
         */
        bool tria_order_ccw;

        /**
         * @brief Sorts the vertex IDs of triangles to Counter-Clockwise order
         * 
         * The order is sorted by calcuting the corss product of the two vectors.
         * If the order was A,B,C, calcuate the dot product of vector AB=B-A=(a1,a2) and BC=C-B=(b1,b2). 
         * If the cross product (0,0,a1b2-a2b1) has a1b2-a2b1>0, then A,B,C is the correct order. 
         * Otherwise, the CCW order is A,C,B.
         */
        void sort_tria_vertex_ids_ccw();

        /**
         * Prints the triangle bar ids to a file with the given file name prefix. 
         * Triangle vertex are in CCW order.
         *
         * @param file_name_prefix The prefix of the file name to be used for printing.
         */
        void print_tria_bar_ids(const char *file_name_prefix);

        /**
         * Prints the triangle bar coordinates to a file with the given file name prefix.
         * Triangle vertex are in CCW order.
         *
         * @param file_name_prefix The prefix of the file name to be used for printing.
         */
        void print_tria_bar_coords(const char *file_name_prefix);

        /**
         * Prints the triangle vertex ids to a file with the given file name prefix.
         * Triangle vertex are in CCW order.
         *
         * @param file_name_prefix The prefix of the file name to be used for printing.
         */
        void print_tria_vertex_ids(const char *file_name_prefix);

        /**
         * Prints the triangle vertex coordinates to a file with the given file name prefix.
         * Triangle vertex are in CCW order.
         *
         * @param file_name_prefix The prefix of the file name to be used for printing.
         */
        void print_tria_vertex_coords(const char *file_name_prefix);

        /**
         * Prints the triangle quality statistics to a file with the given file name prefix.
         *
         * @param file_name_prefix The prefix of the file name to be used for printing.
         */
        void print_tria_quality_stat(const char *file_name_prefix);

        /**
         * Prints the triangles with worst qualities (aspect ratio > 2) to a file with the given file name prefix.
         * Print the triangle id, vertices id's, and vertices coordinates.
         *
         * @param file_name_prefix The prefix of the file name to be used for printing.
         */
        void print_worst_tria(const char *file_name_prefix);

        /**
         * Prints the overall triangle mesh quality statistics.
         */
        void print_tria_quality_stat_overall();

        /**
         * Prints the Voronoi diagram to a file with the given file name prefix.
         *
         * @param file_name_prefix The prefix of the file name to be used for printing.
         */
        void print_voro_diagram(const char *file_name_prefix);

        /**
         * Prints all of the above outputs.
         */
        void do_print_outputs(); //print all of the above

        /**
         * @brief Prints various outputs related to particle coordinates, triangle bar coordinates, and Voronoi diagrams.
         *      For the final mesh: Implement Half-edge data strucutre, and print the boundary vertices in CCW order. 
         */
        void do_print_final_outputs();
        

        //virtual: depend on the algorithm using (dm/cvd/hybrid)
        /**
         * Updates the position of a point using the CVD (Centroidal Voronoi Diagram) method.
         *
         * @param i The index of the point.
         * @param c The Voronoicell_neighbor_2d object associated with the point.
         * @param x The current x-coordinate of the point.
         * @param y The current y-coordinate of the point.
         * @param inner_pt_over_stop_Continue_mvmt_thres_ct_ A reference to the variable tracking the count of inner points crossing the stop or continue movement threshold.
         * @param[in,out] use_random_point If using CVD meshing: Decide whether to move point to the centroid of the needle type triangle it connects to instead, to eliminate the needle triangle.
         * @param[in,out] random_point_x The x-coordinate of the needle type triangle centroid the point should move to.
         * @param[in,out] random_point_y The y-coordinate of the needle type triangle centroid the point should move to.
         */
        virtual void update_pt_position_cvd(int i, voronoicell_neighbor_2d &c,double x, double y, int &inner_pt_over_stop_Continue_mvmt_thres_ct_,
            bool use_random_point, double random_point_x, double random_point_y){};
        
        /**
         * Initializes the mesh algorithm for the given geometry grid index.
         *
         * @param ij The index of the geometry grid.
         */
        virtual void mesh_alg_init(int ij){};

        /**
         * Performs computations and information extraction for the Voronoi diagram and the corresponding Delaunay triangulation.
         * This function is specific to the algorithm being used.
         */
        virtual void voro_compute_retria_extract_info(){}; 
        
        /**
         * Determines if the meshing algorithm requires re-triangulation.
         */
        virtual void determine_reTria(){reTria=true;};
        
        /**
         * Updates the position of a point and determines if re-triangulation is required.
         * This function is specific to the algorithm being used.
         */
        virtual void update_pt_position_and_determine_reTria(){}; 
        
        /**
         * Updates the movement threshold for re-triangulation of geometry grid ij.
         *
         * @param ij The index of the geometry grid.
         */
        virtual void update_retria_movement_thres(int ij){};
        
        /**
         * Determines the method switch in the hybrid algorithm based on the relative change in aspect ratio.
         *
         * @param ar_median_relative_change The relative change in median aspect ratio.
         */
        virtual void hybrid_determine_method_switch(double ar_median_relative_change){};

        //virtual: defined in basic_calculation, for projection()
        
        /**
         * Calculates the spatial step size (deps) used in finite difference in projection() value for the given coordinates.
         *
         * @param x The x-coordinate of the point.
         * @param y The y-coordinate of the point.
         * @param deps_prime The default value of deps.
         * @return The calculated value of deps for the point.
         *         Returns the previous valid deps (deps_prime) if the point is not in the boundary grid.
         */
        double deps_func(double x, double y, double deps_prime); 
        
        /**
         * Calculates the geometry tolerance (geps) function for the given coordinates.
         *
         * @param x The x-coordinate of the point.
         * @param y The y-coordinate of the point.
         * @param geps_prime The default value of geps.
         * @return The calculated value of geps for the point.
         *         Returns the previous valid geps (geps_prime) if the point is not in the boundary grid.
         */
        double geps_func(double x, double y, double geps_prime);
        
        /**
         * Calculates the signed distance function (SDF) for the given coordinates.
         *
         * @param x The x-coordinate of the point.
         * @param y The y-coordinate of the point.
         * @return The calculated value of the SDF for the point.
         */
        double sdf_func(double x, double y);

        /**
         * Checks if a point is inside the outer grid.
         *
         * @param x The x-coordinate of the point.
         * @param y The y-coordinate of the point.
         * @return True if the point is inside the outer grid, false otherwise.
         */
        bool pt_in_outer_grid(double x,double y);
        
        


        //meshing framework structure steps
        /**
         * Initializes the meshing framework.
         * This function sets up the necessary data structures and parameters for the meshing process.
         * It also generate the initial meshing points.
         */
        void meshing_init();

        /**
         * Performs iterations of the meshing algorithm until stopping criteria is reached. 
         */
        void meshing_iter();

        /**
         * Finalizes the meshing process.
         * This function performs any necessary cleanup and print outputs after the meshing is complete.
         */
        void meshing_final();


        //timing variables
        double t_meshing_init;
        double t_add_pt;
        double t_add_pt_centroid;
        double t_voro_computation;
        double t_update_pt_position;
        double t_algorithm_remain;
        double t_algorithm_all_non_voro;

        double t_dm_getBarinfo;
        double t_dm_applyBarForce;
        double t_dm_updatePtPos;

 
        
    };


    
}


#endif



