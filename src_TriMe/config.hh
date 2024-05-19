/** \file config.hh
 * \brief Global configuration file for setting various compile-time options.
 */

#ifndef PARALLEL_MESHING_CONFIG_HH
#define PARALLEL_MESHING_CONFIG_HH

#include <limits>

namespace voro {
    /**
     * @brief The factor used to construct the geometry grid from the container computation grid.
     * gnx=fac*con.nx; gny=fac*con.ny.
     */
    const int geo_grid_nxy_fac_con_grid=5;

    /**
     * @brief Initialize Ncurrent=pt_init_frac*Ntotal points in the container.
     */
    const double pt_init_frac=0.2; 

    /**
     * @brief Add points into the container: Nadd=add_pt_fac*Ncurrent
     */
    const double add_pt_fac=0.6;

    /**
     * @brief Machine epsilon used in calculating deps(i.e.deltaX,deltaY) in projection() step for numerical differentiation.
     */
    const double eps=std::numeric_limits<double>::epsilon();

    /**
     * @brief Geometry boundary tolerance geps=geps_h_frac*h
     */
    const double geps_h_frac=0.01;

    /**
     * @brief Boundary grid adf quad cell error tolerance: max_err=adf_err_geps_frac*geps
     */
    const double adf_err_geps_frac=0.1;

    /**
     * @brief Boundary grid adf quad cell max depth
     */
    const int adf_max_depth=10;

    /**
     * @brief In sizing_2d.hh, when creating bdry points (to be used in finding the medial axis points),
     * generate Npts_grid_init points in each bdry grid (then project these points onto the geometry boundary, to be boundary points)
     */
    const int Npts_grid_init=5;

    /**
     * @brief In customer_shape_2d, normalize the shape so that the model fits into the domain box.
     * The biggest dimension of the model is normalize_model_length_fac of the size of the minimum length side of the box.
     */
    const double normalize_model_length_fac=0.8;
    

    /**
     * @brief In geometry projection(), projection of outside points onto geometry boundary,
     * maximum allowed Newton iteration steps.
     */
    const int max_projection_step=15;
    
    /**
     * @brief Point movement damping parameter used in damped newton projection().
     */
    const double alpha=1.0;


    /**
     * @brief Boundary pts generated in sizing_2d.cc/.hh should be dis>bdry_pt_dis_fac*hbdry away.
     * i.e., not crowding and nicely spaced.
     */
    const double bdry_pt_dis_fac=0.5;

    /**
     * @brief In trimming medial axis pt, the pts surrounding region to check for neighbors is defined by
     *        region: circular region with radius medial_trim_region_check_r=2*medial_cutoff_nei_ct*bdry_tol;
     *        where bdry_tol=bdry_pt_dis_fac*hbdry; which is the characteristic distance between bdry points
     *        In trimming medial axis pts, if the region is sparse,
     *        that is, a point has less than medial_cutoff_nei_ct pts in its surrounding circular region,
     *        delete the pt
     *        i.e. need at least medial_curoff_nei_ct points nearby as neighbors
     */
    const int medial_cutoff_nei_ct=3;

    /**
     * @brief The factor used to define the radius of the circular region for trimming medial axis points.
     *        It is calculated as medial_trim_region_check_r=2*medial_cutoff_nei_ct*bdry_tol.
     *        bdry_tol=bdry_pt_dis_fac*hbdry, which is the characteristic distance between boundary points.
     */
    const double medial_trim_region_fac=2; //2 works well for poker
    
    /**
     * @brief The relationship between density and sizing field: density=1/(sizing^density_sizing_exp)
     */
    const double density_sizing_exp=2.0;

    /**
     * @brief In calculating the density field from the automatic sizing field,
     *        use the average smoothing sizing value over the region surrounding the grid.
     *        The region is defined by i(j)+/-density_sizing_smoothing_grid_ct.
     */
    const int density_sizing_smoothing_grid_ct=1;

    /**
     * @brief The retriangulation criteria: if any point movement in an iteration has dis> retria_fac*h_i, then reTria=true.
     */
    const double retria_fac=0.1;

    /**
     * @brief Stop Continue based on inner point movement: all inner points have movement < stop_Continue_fac*h_i, Continue=false.
     */
    const double stop_Continue_fac=0.001;

    /**
     * @brief Stop Continue based on the relative change in mesh quality measure.
     *        If the relative change in aspect ratio and edge ratio are both smaller than these values, Continue=false.
     */
    const double stop_Continue_ar_mean_rel_change=0.001;
    const double stop_Continue_er_mean_rel_change=0.001;
    const double stop_Continue_max_ar_rel_change=0.005;

    /**
     * @brief In updating point positions, restrict point movement distance to be <= pt_mvmt_dis_thres_fac*h_i.
     */
    const double pt_mvmt_dis_thres_fac=0.4;

    /**
     * @brief Determine addPt_and_reTtria criteria:
     *        If the relative change in the median aspect ratio is smaller than determine_addPt_relative_change_mean_ar, addPt=true.
     *        If the relative change in the median edge ratio is smaller than determine_addPt_relative_change_mean_er, addPt=true.
     *        If inbtw_tria_iter_ct > inbtw_tria_iter_ct_thres, addPt=true, reTria=true.
     *        If the percentage of the number of inner points' movement over stop_Continue criteria is less than
     *        proportion_inner_pt_mvmt_over_thres_thres, addPt=true, reTria=true.
     */
    const double determine_addPt_relative_change_mean_ar = 0.002;
    const double determine_addPt_relative_change_mean_er = 0.002;
    const int inbtw_tria_iter_ct_thres = 10;
    const double proportion_inner_pt_mvmt_over_thres_criteria = 0.4;


    //factor used in valid triangle criteria: 
    //every edge length of the triangle should have length < tria_length_cri_fac*max(h_va,h_vb,h_vc)
    //const double tria_length_cri_fac=10.0;


    /**
     * @brief The factor used in adding walls to Voronoi cells:
     *        Construct a regular octagon with height and width: voro_wall_fac*hi.
     */
    const double voro_wall_fac=3.0;

    /**
     * @brief The factor used in the valid triangle criteria:
     *        If adf_sdf(circumcenter)/circumradius < thres (0.4), return true.
     */
    const double tria_ccc_cri_fac=0.4;

    /**
     * @brief DM related variables
     */
    const double deltat=0.3; //delta_t used in Forward Euler update
    const double Fscale=1.2; //Scale force larger to let points push against geo bdry

    /**
     * @brief CVD related variables
     *        If true, use adaptive subdivision routine for Voronoi cells to compute the centroid.
     *        If false, use only the first level subdivision: triangle ct = number of Voronoi vertices.
     */
    const bool voro_adaptive_subdivide=true; 
    const double voro_subdivide_etol_eps_min=0.001;
    const double fac_thousandth=1.0/1000.0;

    /**
     * @brief Hybrid related variables
     */
    const double hybrid_switch_ar_alpha_mean_relative_change_threshold_1=0.0025;
    const double hybrid_switch_ar_alpha_mean_relative_change_threshold_2=0.0015;
}

#endif
