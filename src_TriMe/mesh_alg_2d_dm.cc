#include "mesh_alg_2d_dm.hh"
#include "parallel_meshing_2d.hh"

using namespace voro;

/**
 * @brief Constructor for the `mesh_alg_2d_dm` class.
 * @param pm2d_ A pointer to the `parallel_meshing_2d` object.
 */
mesh_alg_2d_dm::mesh_alg_2d_dm(parallel_meshing_2d *pm2d_)
    :mesh_alg_2d(pm2d_), 
    retria_movement_thres(new double[gnxy]),
    xy_id_previous_tria(new double[2*(pm2d_->Ntotal)]),
    barinfo(new double[1]),pt_mvmt_btw_tria_large_ct(0),
    barinfo_barct(0)
    {}

/**
 * @brief Destructor for the `mesh_alg_2d_dm` class.
 * Cleans up dynamically allocated memory.
 */
mesh_alg_2d_dm::~mesh_alg_2d_dm(){
    delete [] retria_movement_thres;
    delete [] xy_id_previous_tria;
    delete [] barinfo;
}

/**
 * @brief Updates the retriangulation movement threshold for a given geometry grid index.
 * If the geometry grid index is an inner or boundary grid, the threshold is updated based on the characteristic length.
 * 
 * @param ij The geometry grid index.
 */
void mesh_alg_2d_dm::update_retria_movement_thres(int ij){
    if(geo_grid(ij)<=gnxy){ //inner or bdry geo_grid
        retria_movement_thres[ij]=retria_fac*chrtrt_len_h(ij);
    }
}

/**
 * @brief Initializes the mesh algorithm for a given geometry grid index.
 * Initializes the retriangulation threshold and updates it based on the geometry sizing grid.
 * 
 * @param ij The geometry grid index.
 */
void mesh_alg_2d_dm::mesh_alg_init(int ij){
    //re-triangulation threshold init: adaptive to geometry sizing grid
    retria_movement_thres[ij]=-1;
    update_retria_movement_thres(ij);
}


/**
 * @brief Checks if there is any large movement of a point since the previous triangulation.
 * If any point has a movement larger than the threshold, |p_i-p_i_old| > retria_movement_thres[i_old], reTria is set to true.
 * 
 * @param i The index of the point.
 * @param pt_mvmt_btw_tria_large_ct_ The count of points with large movement between triangulations.
 */
void mesh_alg_2d_dm::cri_reTria_movement_previous_tria(int i, int &pt_mvmt_btw_tria_large_ct_){
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
}

/**
 * @brief Determines whether to re-triangulate or not based on various conditions.
 * If reTria is false, it checks if an initial triangulation is needed in the first iteration.
 * In later iterations, it checks for large point movements from the previous triangulation.
 */
void mesh_alg_2d_dm::determine_reTria(){
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
            }
        }
    }
    
}

/**
 * @brief Retrieves information about the bars.
 * This function calculates the length and other properties of the bars and stores them in the barinfo array.
 * It also calculates the sum of squared bar lengths (sumL2) and the sum of squared desired bar lengths (sumhbar2).
 */
void mesh_alg_2d_dm::getBarinfo(){
    #pragma omp for reduction(+ : sumL2,sumhbar2) 
    for(int i=0;i<bar_ct;i++){
        int i2=2*i;
        int i4=4*i;
        int barid_i2_1_2=2*barid[i2+1];
        int barid_i2_2=2*barid[i2];
        
        barinfo[i4]=pm2d->xy_id[barid_i2_1_2]-pm2d->xy_id[barid_i2_2]; //x component
        barinfo[i4+1]=pm2d->xy_id[barid_i2_1_2+1]-pm2d->xy_id[barid_i2_2+1]; //y component
        double barL2=sqr(barinfo[i4])+sqr(barinfo[i4+1]);
        barinfo[i4+2]=sqrt(barL2); //bar length  
        //barinfo[i4+3]=pm2d->sizing((pm2d->xy_id[2*barid[i2+1]]+pm2d->xy_id[2*barid[i2]])*0.5,(pm2d->xy_id[2*barid[i2+1]+1]+pm2d->xy_id[2*barid[i2]+1])*0.5);
        barinfo[i4+3]=get_chrtrt_len_h((pm2d->xy_id[barid_i2_1_2]+pm2d->xy_id[barid_i2_2])*0.5,(pm2d->xy_id[barid_i2_1_2+1]+pm2d->xy_id[barid_i2_2+1])*0.5);  //hbars
        sumL2 += barL2; sumhbar2 += sqr(barinfo[i4+3]);
    }
}

/**
 * @brief Applies bar forces to update the node positions.
 * This function calculates the force vectors based on the bars' properties and applies them to update the node positions.
 */
void mesh_alg_2d_dm::applyBarForce(){ 

    double force_fac=Fscale*sqrt(sumL2/sumhbar2);

    //update xy_id_new to be the same as current pt positions
    #pragma omp for
    for(int i=0;i<Ncurrent;i++){
        int i2=2*i;
        xy_id_new[i2]=pm2d->xy_id[i2];
        xy_id_new[i2+1]=pm2d->xy_id[i2+1];
    }

    //calculate force vectors (x, y components)
    //update xy_id_new
    #pragma omp for
    for(int i=0;i<bar_ct;i++){
        int i4=4*i;
        int i2=2*i;
        double LL = barinfo[i4+3]*force_fac;
        if(LL-barinfo[i4+2]>0.0){
            int id2 = barid[i2+1];
            int id1 = barid[i2];
            int id22=2*id2;
            int id12=2*id1;

            double F= deltat*(LL/barinfo[i4+2]-1.0);
            double ffx = F*barinfo[i4];   //Ftot_x
            double ffy = F*barinfo[i4+1]; //Ftot_y

            //update node positions: by sum of all force vectors from all bars meeting at a node
            #pragma omp atomic
            xy_id_new[id22]+=ffx;   
            #pragma omp atomic
            xy_id_new[id22+1]+=ffy;
            
            #pragma omp atomic
            xy_id_new[id12]+=(-ffx);
            #pragma omp atomic
            xy_id_new[id12+1]+=(-ffy);
            
        }
    }


}

/**
 * @brief Computes Voronoi and the corresponding triangulation, extracts relevant information, and performs output operations.
 * This function computes Voronoi information, extracts the required information based on the specified flags, and performs output operations.
 */
void mesh_alg_2d_dm::voro_compute_retria_extract_info(){

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

/**
 * @brief Updates the point positions and determines whether retriangulation is required.
 * This function updates the point positions based on the calculated bar forces and applies various treatments and projections to get the new point positions.
 * It also determines whether retriangulation is necessary based on the movement of points since the previous triangulation.
 */
void mesh_alg_2d_dm::update_pt_position_and_determine_reTria(){
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


