#ifndef MESH_ALG_2D_HYBRID_HH
#define MESH_ALG_2D_HYBRID_HH

#include "mesh_alg_2d_dm.hh"
#include "mesh_alg_2d_cvd.hh"


namespace voro {

    /**
     * @brief A class for hybrid 2D meshing algorithm.
     *
     * The mesh_alg_2d_hybrid class combines the functionality of the mesh_alg_2d_dm and mesh_alg_2d_cvd classes
     * to perform hybrid meshing in 2D. It inherits from both classes and introduces additional variables and functions
     * specific to the hybrid meshing approach.
     */
    class mesh_alg_2d_hybrid : public mesh_alg_2d_dm, public mesh_alg_2d_cvd {
      public:
        /**
         * @brief Constructor for mesh_alg_2d_hybrid class.
         * @param pm2d_ Pointer to parallel_meshing_2d object.
         */
        mesh_alg_2d_hybrid(parallel_meshing_2d *pm2d_);

        /**
         * @brief Destructor for mesh_alg_2d_hybrid class.
         */
        ~mesh_alg_2d_hybrid();
        
        //Hybrid related variables and functions
        bool dm_switch_to_cvd; /**< Flag indicating whether to switch from DM to CVD method. */
        int cvd_iter_ct; /**< Number of iterations for CVD method. */
        bool just_use_cvd; /**< Flag indicating to if CVD is used in the previous iteration. */
        bool just_use_dm; /**< Flag indicating to if DM is used in the previous iteration. */
        /**
         * @brief Function to calculate the count of large movements between triangulations for a particle.
         * @param i Particle index.
         * @param pt_mvmt_btw_tria_large_ct_ Reference to store the count of points with large movements.
         */
        void cri_reTria_movement_previous_tria_cvd_mvmt_fac(int i, int &pt_mvmt_btw_tria_large_ct_);
        
        //virtual functions in base class
        /**
         * @brief Function to initialize mesh algorithm for a grid cell.
         * @param ij Grid cell index.
         */
        void mesh_alg_init(int ij) override;
        /**
         * @brief Function to update particle positions and determine reTria flag.
         */
        void update_pt_position_and_determine_reTria() override; //need to update reTria=true;
        /**
         * @brief Function to compute Voronoi cells and extract information.
         */
        void voro_compute_retria_extract_info() override; //Obtain Voro vertices info
        
        // DM specific functions
        /**
         * @brief Function to determine reTria flag for a grid cell using DM method.
         */
        void determine_reTria() override;
        /**
         * @brief Function to update reTria movement thresholds for a grid cell.
         * @param ij Grid cell index.
         */
        void update_retria_movement_thres(int ij) override;
        
        // Hybrid specific functions
        /**
         * @brief Function to determine the method switch between DM and CVD based on the relative change in alpha-mean of aspect ratio.
         * @param ar_alpha_mean_relative_change Relative change in alpha-mean of aspect ratio.
         */
        void hybrid_determine_method_switch(double ar_alpha_mean_relative_change) override; //determine dm_switch_to_cvd

    };
}



#endif



