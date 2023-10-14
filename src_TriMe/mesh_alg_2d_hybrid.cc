#include "mesh_alg_2d_hybrid.hh"
#include "parallel_meshing_2d.hh"

using namespace voro;


/**
 * @brief Constructor for the `mesh_alg_2d_hybrid` class.
 *
 * @param pm2d_ A pointer to the `parallel_meshing_2d` object.
 */
mesh_alg_2d_hybrid::mesh_alg_2d_hybrid(parallel_meshing_2d *pm2d_)
    :mesh_alg_2d_cvd(pm2d_),mesh_alg_2d_dm(pm2d_),mesh_alg_2d(pm2d_)
    {
        dm_switch_to_cvd=false;
        cvd_iter_ct=0;
        max_tria_ar_prev=-1;
        just_use_cvd=false;
        just_use_dm=false;
    }

/**
 * @brief Destructor for the mesh_alg_2d_hybrid class.
 */
mesh_alg_2d_hybrid::~mesh_alg_2d_hybrid(){
}

/**
 * @brief Updates the point movement threshold for retriangulation at the grid index.
 *
 * This function updates the movement threshold for retriangulation at the provided grid index.
 * If the grid index corresponds to an inner or boundary geometric grid, the movement threshold is set
 * as a product of the retriangulation factor and the characteristic length of the grid cell.
 *
 * @param ij The grid index.
 */
void mesh_alg_2d_hybrid::update_retria_movement_thres(int ij){
    if(geo_grid(ij)<=gnxy){ //inner or bdry geo_grid
        retria_movement_thres[ij]=retria_fac*chrtrt_len_h(ij);
    }
}

/**
 * @brief Determines the method switch between DM and CVD in the hybrid meshing algorithm.
 *
 * This function determines whether to switch between the DM method and the CVD method
 * based on the provided alpha mean relative change in aspect ratio. If `dm_switch_to_cvd` is currently false (indicating the DM
 * method is being used), the function checks if the alpha mean relative change is below a threshold. If so, it sets
 * `dm_switch_to_cvd` to true, indicating a switch to the CVD method. It also sets other related flags accordingly.
 * 
 * If `dm_switch_to_cvd` is currently true (indicating the CVD method is being used), the function checks if the
 * alpha mean relative change is above a threshold. If so, it sets `dm_switch_to_cvd` to false, indicating a switch
 * back to the DM method. It also sets other related flags accordingly.
 *
 * @param ar_alpha_mean_relative_change The alpha mean relative change value in aspect ratio.
 */
void mesh_alg_2d_hybrid::hybrid_determine_method_switch(double ar_alpha_mean_relative_change){
    if(dm_switch_to_cvd==false){
        if(Ncurrent==Ntotal){
            if(ar_alpha_mean_relative_change<hybrid_switch_ar_alpha_mean_relative_change_threshold_2){
                dm_switch_to_cvd=true;
                Continue=true;
                just_use_cvd=true;
                
                addPt=false;
                reTria=false;
            }
        }
        /*
        //Continue=true;
        if(Ncurrent<Ntotal){
            if(ar_alpha_mean_relative_change<hybrid_switch_ar_alpha_mean_relative_change_threshold_1){
                dm_switch_to_cvd=true;
                Continue=true;
                just_use_cvd=true;
                
                addPt=false;
                reTria=false;
            }
        }
        else{
            if(ar_alpha_mean_relative_change<hybrid_switch_ar_alpha_mean_relative_change_threshold_2){
                dm_switch_to_cvd=true;
                Continue=true;
                just_use_cvd=true;
                
                addPt=false;
                reTria=false;
            }
        }
        */
/*
        if(ar_alpha_mean_relative_change<0.0025){
            printf("dm switch to cvd \n");
            dm_switch_to_cvd=true;
            Continue=true;
            just_use_cvd=true;
            
            addPt=false;
            reTria=false;
            
            if(max_tria_ar_prev>0){
                if(max_tria_ar>max_tria_ar_prev){
                    double max_tria_ar_rel_change=fabs(max_tria_ar-max_tria_ar_prev)/max_tria_ar_prev;
                    if(max_tria_ar_rel_change>0.05){
                        printf("dm switch to cvd \n");
                        dm_switch_to_cvd=true;
                        Continue=true;
                        just_use_cvd=true;
                        
                        addPt=false;
                        reTria=false;
                    }
                }
            }
            
        }
*/
        /*
        //see if need to switch to CVD, to stablize worst quality triangle
        if(ar_alpha_mean_relative_change<hybrid_switch_ar_alpha_mean_relative_change_threshold){
            if(max_tria_ar_prev>0){
                if(max_tria_ar>max_tria_ar_prev){
                    double max_tria_ar_rel_change=fabs(max_tria_ar-max_tria_ar_prev)/max_tria_ar_prev;
                    if(max_tria_ar_rel_change>hybrid_switch_max_tria_ar_rel_change_threshold){
                        dm_switch_to_cvd=true;
                        Continue=true;
                        just_use_cvd=true;
                        
                        addPt=false;
                        reTria=false;
                    }
                }
            }
        }
        */
    }
    if(dm_switch_to_cvd==true){
        /*
        if(Ncurrent<Ntotal){
            if(ar_alpha_mean_relative_change>hybrid_switch_ar_alpha_mean_relative_change_threshold_1){
                dm_switch_to_cvd=false;
                just_use_dm=true;

                addPt=false;
                reTria=false;
            }
        }
        else{
            if(ar_alpha_mean_relative_change>hybrid_switch_ar_alpha_mean_relative_change_threshold_2){
                dm_switch_to_cvd=false;
                just_use_dm=true;

                addPt=false;
                reTria=false;
            }
        }
        */
        
/*
        if(ar_alpha_mean_relative_change>0.0025){
            dm_switch_to_cvd=false;
            just_use_dm=true;

            addPt=false;
            reTria=false;

            /*
            if(max_tria_ar_prev>0){
                double max_tria_ar_rel_change=fabs(max_tria_ar-max_tria_ar_prev)/max_tria_ar_prev;
                if(max_tria_ar_rel_change<0.01){
                    dm_switch_to_cvd=false;
                    just_use_dm=true;

                    addPt=false;
                    reTria=false;
                }
            }
            
        }
*/
        
        /*
        //termination only when using DM
        Continue=true;
        
        //see if need to switch back to DM: if worst quality tria is stablized
        if(max_tria_ar<max_tria_ar_prev){
            double max_tria_ar_rel_change=fabs(max_tria_ar-max_tria_ar_prev)/max_tria_ar_prev;
            if(max_tria_ar_rel_change<hybrid_switch_max_tria_ar_rel_change_threshold){
                dm_switch_to_cvd=false;
                just_use_dm=true;

                addPt=false;
                reTria=false;
            }
        }
        */
    }
    
}

/**
 * @brief Initializes the meshing algorithm for the hybrid method.
 *
 * This function initializes the meshing algorithm by setting the re-triangulation threshold for the specified geometry
 * grid index `ij`. The re-triangulation threshold is initialized to -1 and then updated using the
 * `update_retria_movement_thres` function.
 *
 * @param ij The index of the geometry grid.
 */
void mesh_alg_2d_hybrid::mesh_alg_init(int ij){
    //re-triangulation threshold init: adaptive to geometry sizing grid
    retria_movement_thres[ij]=-1;
    update_retria_movement_thres(ij);
}


/**
 * @brief Determines whether re-triangulation is required.
 *
 * This function determines whether re-triangulation is necessary based on the current state of the meshing algorithm.
 * If `reTria` is already true, indicating that re-triangulation has been requested, no further action is taken.
 * Otherwise, if it is the first iteration (`all_iter_ct == 1`), re-triangulation is required to obtain an initial triangulation.
 * For later iterations, the function checks if there has been any large movement of points since the previous triangulation.
 * If any point |p_i - p_i_old| exceeds the threshold `retria_movement_thres[i_old]`, re-triangulation is required.
 *
 * If re-triangulation is triggered, the `just_use_cvd` flag is reset if it was previously set to true.
 */
void mesh_alg_2d_hybrid::determine_reTria(){
    if(reTria==false){
        //first iteration, need an initial triangulation
        if(all_iter_ct==1){ 
            reTria=true;
        }else{ //later iterations, need to test point movement from previous triangulation
            //if there is any large movement of pt since previous triangulation, 
            //i.e. any point |p_i-p_i_old| > retria_movement_thres[i_old], 
            //reTria=true
            if(pt_mvmt_btw_tria_large_ct>0){
                reTria=true;
                if(just_use_cvd==true){
                    just_use_cvd=false;
                }
            }
        }
    }
}



/**
 * @brief Checks if there is any large movement of a point since the previous triangulation.
 *
 * This function compares the movement of a point with its position in the previous triangulation.
 * If the movement exceeds the local retria_movement_thres[i], it increments the pt_mvmt_btw_tria_large_ct_ variable.
 * It also calculates the prev_pt_mvmt_chrtrt_fac for the given point, used in Voronoi centroid approximation tolerance calculation.
 *
 * @param i The index of the point.
 * @param pt_mvmt_btw_tria_large_ct_ Reference to the variable holding the count of points with large movements.
 */
void mesh_alg_2d_hybrid::cri_reTria_movement_previous_tria_cvd_mvmt_fac(int i, int &pt_mvmt_btw_tria_large_ct_){
    int i2=2*i;
    //count number of pts whose movement since previous triangulation
    //exceeds the local retria_move_thres[i]
    double x_old=xy_id_previous_tria[i2];
    double y_old=xy_id_previous_tria[i2+1];
    double movement=d_points(x_old,y_old,pm2d->xy_id[i2],pm2d->xy_id[i2+1]);
    int i_old=(x_old-ax)*inv_gdx;
    int j_old=(y_old-ay)*inv_gdy;
    int ij_old=j_old*gnx+i_old;
    if(movement>retria_movement_thres[ij_old]){
        pt_mvmt_btw_tria_large_ct_+=1;
    }
    prev_pt_mvmt_chrtrt_fac[i]=std::round(movement/chrtrt_len_h(ij_old) * 1000.0) * fac_thousandth;
    
}



/**
 * @brief Updates the point positions and determines whether re-triangulation is needed.
 *
 * This function updates the positions of the points based on the selected algorithm (DM or CVD) and determines whether
 * re-triangulation is required. It uses a hybrid approach where DM is used initially, and then switches to CVD when 
 * Ncurrent==Ntotal and the relative change in median aspect ratio is small.
 */
void mesh_alg_2d_hybrid::update_pt_position_and_determine_reTria(){

    if(just_use_dm==true){
        just_use_dm=false;
        //pt position was already updated in Voro compute part. 
        //update them in xy_id array
        #pragma omp parallel for num_threads(num_t)
        for(int i=Nfixed;i<Ncurrent;i++){
            int i2=2*i;
            pm2d->xy_id[i2]=xy_id_new[i2];
            pm2d->xy_id[i2+1]=xy_id_new[i2+1];
        }

        //determine reTria: CVD always true
        reTria=true;
    }
    else if(just_use_cvd){ //transition iteration DM: next iter will be CVD
        sumL2=0.0; sumhbar2=0.0;

        //if there is any large movement of pt since previous triangulation, 
        //i.e. any point |p_i-p_i_old| > retria_movement_thres[i_old], 
        //reTria=true
        pt_mvmt_btw_tria_large_ct=0;

    double t0=omp_get_wtime();
        #pragma omp parallel num_threads(num_t)
        {
            //compute info related to the triangle bars
            getBarinfo();
        }
    double t1=omp_get_wtime();

        #pragma omp parallel num_threads(num_t)
        {
            //calculate and apply bar force to find new point positions: store in xy_id_new
            applyBarForce();
        }
    double t2=omp_get_wtime();

        #pragma omp parallel num_threads(num_t)
        {
            //Loop through pts and apply pt movement treatment and projections, 
            //and get final pt positions, store in xy_id
            //also, if Ncurrent==Ntotal, count number of inner pts over stop_Continue movement criteria
            #pragma omp for reduction(+:inner_pt_over_stop_Continue_mvmt_thres_ct, pt_mvmt_btw_tria_large_ct)
            for(int i=Nfixed;i<Ncurrent;i++){
                int i2=2*i;
                double px_old=pm2d->xy_id[i2];
                double py_old=pm2d->xy_id[i2+1];
                double px_new=xy_id_new[i2];
                double py_new=xy_id_new[i2+1];

                double final_x,final_y;
                pt_mvmt_treatment_projection(i,px_old,py_old,px_new,py_new,final_x,final_y);
                pm2d->xy_id[i2]=final_x;
                pm2d->xy_id[i2+1]=final_y;


                //test if inner point pi movement larger than stop_Continue_mvmt_thres
                if(inner_pt_movement_over_stop_Continue_mvmt_thres(i,px_old,py_old,final_x,final_y)==true){
                    inner_pt_over_stop_Continue_mvmt_thres_ct+=1;
                }

                cri_reTria_movement_previous_tria_cvd_mvmt_fac(i,pt_mvmt_btw_tria_large_ct);
            }
        }
    double t3=omp_get_wtime();

    t_dm_getBarinfo+=t1-t0;
    t_dm_applyBarForce+=t2-t1;
    t_dm_updatePtPos+=t3-t2;

        reTria=false;
    }
    else if(dm_switch_to_cvd==false){
        sumL2=0.0; sumhbar2=0.0;

        //if there is any large movement of pt since previous triangulation, 
        //i.e. any point |p_i-p_i_old| > retria_movement_thres[i_old], 
        //reTria=true
        pt_mvmt_btw_tria_large_ct=0;

    double t0=omp_get_wtime();
        #pragma omp parallel num_threads(num_t)
        {
            //compute info related to the triangle bars
            getBarinfo();
        }
    double t1=omp_get_wtime();

        #pragma omp parallel num_threads(num_t)
        {
            //calculate and apply bar force to find new point positions: store in xy_id_new
            applyBarForce();
        }
    double t2=omp_get_wtime();

        #pragma omp parallel num_threads(num_t)
        {
            //Loop through pts and apply pt movement treatment and projections, 
            //and get final pt positions, store in xy_id
            //also, if Ncurrent==Ntotal, count number of inner pts over stop_Continue movement criteria
            #pragma omp for reduction(+:inner_pt_over_stop_Continue_mvmt_thres_ct, pt_mvmt_btw_tria_large_ct)
            for(int i=Nfixed;i<Ncurrent;i++){
                int i2=2*i;
                double px_old=pm2d->xy_id[i2];
                double py_old=pm2d->xy_id[i2+1];
                double px_new=xy_id_new[i2];
                double py_new=xy_id_new[i2+1];

                double final_x,final_y;
                pt_mvmt_treatment_projection(i,px_old,py_old,px_new,py_new,final_x,final_y);
                pm2d->xy_id[i2]=final_x;
                pm2d->xy_id[i2+1]=final_y;

                //test if inner point pi movement larger than stop_Continue_mvmt_thres
                if(inner_pt_movement_over_stop_Continue_mvmt_thres(i,px_old,py_old,final_x,final_y)==true){
                    inner_pt_over_stop_Continue_mvmt_thres_ct+=1;
                }

                cri_reTria_movement_previous_tria(i,pt_mvmt_btw_tria_large_ct);
            }
        }
    double t3=omp_get_wtime();

    t_dm_getBarinfo+=t1-t0;
    t_dm_applyBarForce+=t2-t1;
    t_dm_updatePtPos+=t3-t2;

        //determine reTria: DM false
        reTria=false;
    }
    else{
        //pt position was already updated in Voro compute part. 
        //update them in xy_id array
        #pragma omp parallel for num_threads(num_t)
        for(int i=Nfixed;i<Ncurrent;i++){
            int i2=2*i;
            pm2d->xy_id[i2]=xy_id_new[i2];
            pm2d->xy_id[i2+1]=xy_id_new[i2+1];
        }

        //determine reTria: CVD always true
        reTria=true;
        
    }
} 



/**
 * @brief Computes Voronoi information and extracts necessary information for the meshing method.
 *
 * This function computes the Voronoi information and extracts the required information. 
 * Depending on the selected algorithm (DM or CVD), different information is obtained.
 * For DM, bar information is obtained, while for CVD, Voronoi cell vertices information (xy & ct) is obtained.
 */
void mesh_alg_2d_hybrid::voro_compute_retria_extract_info(){

    if(dm_switch_to_cvd==true){
//printf("cvd\n");
        determine_printOutputs();

        bool store_barid=false;
        bool store_triangle_quality=false;
        bool store_tria_vertex=false;
        bool store_tria_centroid=false;
        bool store_ccrd_h_ratio=false;
        bool store_max_tria_ar=true;
        bool cvd_pt_update=true; //indicating to use update_pt_position_cvd(i,...) in each particle i's voronoi computation
        bool store_xy_id=false;
        double *dummy_array=new double[1];
        if(printOutputs==true){
            store_triangle_quality=true;
            store_tria_vertex=true;
            store_barid=true;
        }

        voro_compute_and_store_info(store_barid,store_triangle_quality,
            store_tria_vertex,
            store_tria_centroid,
            store_ccrd_h_ratio, store_max_tria_ar,
            cvd_pt_update,
            store_xy_id, dummy_array);

        delete [] dummy_array;

        //output
        if(printOutputs==true){
            do_print_outputs();
        }
        cvd_iter_ct++;
        
    }
    else{
//printf("dm\n");

        determine_printOutputs();

        bool store_barid=true;
        bool store_triangle_quality=false;
        bool store_tria_vertex=false;
        bool store_tria_centroid=false;
        bool store_ccrd_h_ratio=false;
        bool store_max_tria_ar=true; 
        bool cvd_pt_update=false;
        bool store_xy_id=true;  //store xy_id in xy_id_previous_tria array.
        if(printOutputs==true){
            store_tria_vertex=true;
            store_triangle_quality=true;
        }

        voro_compute_and_store_info(store_barid,store_triangle_quality,
            store_tria_vertex,
            store_tria_centroid,
            store_ccrd_h_ratio, store_max_tria_ar,
            cvd_pt_update,
            store_xy_id,
            xy_id_previous_tria
            );

        if(bar_ct>barinfo_barct){
            delete [] barinfo;
            barinfo=new double[4*bar_ct];
            barinfo_barct=bar_ct;
        }

        //output
        if(printOutputs==true){
            do_print_outputs();
        }

    }
}















