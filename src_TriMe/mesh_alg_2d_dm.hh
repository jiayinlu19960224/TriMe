#ifndef MESH_ALG_2D_DM_HH
#define MESH_ALG_2D_DM_HH


#include "mesh_alg_2d.hh"
//#include "parallel_meshing_2d.hh"

namespace voro {

    /**
     * @class mesh_alg_2d_dm
     * @brief A class for performing meshing algorithms in 2D using the DistMesh method.
     * This class extends the base class mesh_alg_2d and implements additional functionality specific to the DistMesh method.
     * It inherits virtually from mesh_alg_2d to enable polymorphic behavior.
     */
    class mesh_alg_2d_dm : virtual public mesh_alg_2d {
      public:
        /**
         * @brief Constructor for mesh_alg_2d_dm.
         * @param pm2d_ A pointer to the parallel_meshing_2d object.
         */
        mesh_alg_2d_dm(parallel_meshing_2d *pm2d_);
        /**
         * @brief Destructor for mesh_alg_2d_dm.
         */
        ~mesh_alg_2d_dm();
        
        //DM related variables
        double *retria_movement_thres;  /**< -1:dummy; if any point has movement larger than the threshold, reTria=true */
        double *xy_id_previous_tria; /**< Point positions at the previous triangulation */
        double *barinfo;  /**< Bar information */
        int barinfo_barct; /**< Number of bars */
        double sumL2; /**< Sum of squared bar lengths */
        double sumhbar2; /**< Sum of desired bar lengths */
        int pt_mvmt_btw_tria_large_ct; /**< Count of points with large movement between triangulations */

        /**
         * @brief Check if there is large point movement since the previous triangulation.
         * @param i Index of the point to check.
         * @param pt_mvmt_btw_tria_large_ct_ Reference to the variable holding the count of points with large movement.
         */
        void cri_reTria_movement_previous_tria(int i, int &pt_mvmt_btw_tria_large_ct_);

        /**
         * @brief Get triangle mesh bars information.
         */
        void getBarinfo();

        /**
         * @brief Apply bar force based on the DM algorithm.
         */
        void applyBarForce();


        //virtual functions in base class
        /**
         * @brief Initialize the mesh algorithm parameters and initial points in 2D.
         * @param ij The geometry grid index.
         */
        void mesh_alg_init(int ij) override;
        /**
         * @brief Determine if re-triangulation is needed based on certain criteria.
         */
        void determine_reTria() override;

        /**
         * @brief Update the point positions and determine if re-triangulation is needed.
         */
        void update_pt_position_and_determine_reTria() override; 

        /**
         * @brief Perform Voronoi computation and extract information related to Voronoi diagram and the corresponding triangulation.
         */
        void voro_compute_retria_extract_info() override; 

        /**
         * @brief Update the threshold for re-triangulation of the geometry grid ij.
         * @param ij The geometry grid index.
         */
        void update_retria_movement_thres(int ij) override;

    };
}


#endif



