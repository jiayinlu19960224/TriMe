#include "mesh_alg_2d_cvd.hh"
#include "parallel_meshing_2d.hh"

using namespace voro;

/**
 * @brief Constructs a new 2D CVD meshing algorithm object.
 *
 * This constructor initializes a meshing algorithm object specific to CVD (Centroidal Voronoi Diagram) meshing
 * for 2D meshes. It takes a pointer to a parallel_meshing_2d object as a parameter.
 *
 * For CVD meshing, density information is needed for outer grids as well, since Voronoi vertices can lie on outer grids.
 * This constructor calls the `get_ogrid_sizing_density_grid()` method of the `size_field` object in the `pm2d` parameter
 * to obtain the density information for outer grids. It also dynamically allocates memory for the `prev_pt_mvmt_chrtrt_fac`
 * array based on the total number of meshing point in `pm2d_`.
 *
 * @param pm2d_ A pointer to the parallel_meshing_2d object.
 */
mesh_alg_2d_cvd::mesh_alg_2d_cvd(parallel_meshing_2d *pm2d_)
    :mesh_alg_2d(pm2d_)
    {
        //For CVD meshing, need density info for outer grids too. 
        //Since Voronoi vertices can lie on outer grids.
        pm2d->size_field->get_ogrid_sizing_density_grid();
        prev_pt_mvmt_chrtrt_fac=new double[pm2d_->Ntotal];
    }

/**
 * @brief Destructor for the 2D CVD meshing algorithm object.
 *
 * This destructor is responsible for releasing the memory allocated for the `prev_pt_mvmt_chrtrt_fac` array.
 */
mesh_alg_2d_cvd::~mesh_alg_2d_cvd(){
    delete [] prev_pt_mvmt_chrtrt_fac;
}

/**
 * @brief Calculate the Voronoi cell centroid from triangles.
 *
 * This function calculates the Voronoi cell centroid by iterating over the triangles that make up the cell.
 * It takes the vertices of the triangles, along with the variables `centroid_x_` and `centroid_y_` as references,
 * which will store the computed centroid coordinates.
 *
 * @param tria_make_ct The number of triangles that make up the Voronoi cell.
 * @param tria_vertices_x The x-coordinates of the vertices of the triangles.
 * @param tria_vertices_y The y-coordinates of the vertices of the triangles.
 * @param centroid_x_ Reference to a variable to store the x-coordinate of the centroid.
 * @param centroid_y_ Reference to a variable to store the y-coordinate of the centroid.
 */
void mesh_alg_2d_cvd::calculate_voro_centroid_from_triangles(int tria_make_ct,
    std::vector<double> &tria_vertices_x, std::vector<double> &tria_vertices_y,
    double &centroid_x_, double &centroid_y_){
    double voro_integrand_rho=0.0;
    double voro_integrand_xrho=0.0;
    double voro_integrand_yrho=0.0;
    for(int triai=0;triai<tria_make_ct;triai++){
        int triai_3=3*triai;
        int triai_3_1=triai_3+1;
        int triai_3_2=triai_3+2;
        double x1_=tria_vertices_x[triai_3];
        double y1_=tria_vertices_y[triai_3];
        double x2_=tria_vertices_x[triai_3_1];
        double y2_=tria_vertices_y[triai_3_1];
        double x3_=tria_vertices_x[triai_3_2];
        double y3_=tria_vertices_y[triai_3_2];
        cvd_tria_info* trias=new cvd_tria_info(x1_,y1_,x2_,y2_,x3_,y3_,pm2d);
        voro_integrand_rho+=trias->integrand_rho();
        voro_integrand_xrho+=trias->integrand_xrho();
        voro_integrand_yrho+=trias->integrand_yrho();
        delete trias;
    }

    //centroid calculation
    double inv_voro_integrand_rho=1.0/voro_integrand_rho;
    centroid_x_=voro_integrand_xrho*inv_voro_integrand_rho;
    centroid_y_=voro_integrand_yrho*inv_voro_integrand_rho;
}

/**
 * @brief Calculate the Voronoi cell centroid from triangles.
 *
 * This function calculates the Voronoi cell centroid by iterating over the triangles that make up the cell.
 * It takes the vertices of the triangles, along with the variables `centroid_x_` and `centroid_y_` as references,
 * which will store the computed centroid coordinates. The function also computes the Voronoi cell area and stores
 * it in the `voro_area_` variable.
 *
 * @param tria_make_ct The number of triangles that make up the Voronoi cell.
 * @param tria_vertices_x The x-coordinates of the vertices of the triangles.
 * @param tria_vertices_y The y-coordinates of the vertices of the triangles.
 * @param centroid_x_ Reference to a variable to store the x-coordinate of the centroid.
 * @param centroid_y_ Reference to a variable to store the y-coordinate of the centroid.
 * @param voro_area_ Reference to a variable to store the area of the Voronoi cell.
 */
void mesh_alg_2d_cvd::calculate_voro_centroid_from_triangles(int tria_make_ct,
    std::vector<double> &tria_vertices_x, std::vector<double> &tria_vertices_y,
    double &centroid_x_, double &centroid_y_,double &voro_area_){
    double voro_integrand_rho=0.0;
    double voro_integrand_xrho=0.0;
    double voro_integrand_yrho=0.0;
    voro_area_=0.0;
    for(int triai=0;triai<tria_make_ct;triai++){
        int triai_3=3*triai;
        int triai_3_1=triai_3+1;
        int triai_3_2=triai_3+2;
        double x1_=tria_vertices_x[triai_3];
        double y1_=tria_vertices_y[triai_3];
        double x2_=tria_vertices_x[triai_3_1];
        double y2_=tria_vertices_y[triai_3_1];
        double x3_=tria_vertices_x[triai_3_2];
        double y3_=tria_vertices_y[triai_3_2];
        cvd_tria_info* trias=new cvd_tria_info(x1_,y1_,x2_,y2_,x3_,y3_,pm2d);
        voro_integrand_rho+=trias->integrand_rho();
        voro_integrand_xrho+=trias->integrand_xrho();
        voro_integrand_yrho+=trias->integrand_yrho();
        voro_area_+=trias->_area;
        delete trias;
    }

    //centroid calculation
    double inv_voro_integrand_rho=1.0/voro_integrand_rho;
    centroid_x_=voro_integrand_xrho*inv_voro_integrand_rho;
    centroid_y_=voro_integrand_yrho*inv_voro_integrand_rho;
}


/**
 * @brief Computes the adaptive Voronoi centroid for a particle.
 * 
 * This function calculates the adaptive Voronoi centroid for a particle given its ID (`pid`), and stores the computed
 * centroid coordinates in `cx` and `cy`. The function takes the number of Voronoi vertices (`voro_vertex_ct`), and the
 * X and Y coordinates of the Voronoi vertices (`voro_vertex_xy`). It performs adaptive triangle subdivision of the Voronoi cell
 * based on an error tolerance on the accuracy of the Voronoi centroid approximation. The Voronoi cell is subdivided
 * into triangles, and a quadrature rule is used on each triangle to compute the centroid. The final centroid is updated
 * in `cx` and `cy`.
 * 
 * @param pid Particle ID.
 * @param cx Computed X-coordinate of the adaptive Voronoi centroid.
 * @param cy Computed Y-coordinate of the adaptive Voronoi centroid.
 * @param voro_vertex_ct Number of Voronoi vertices.
 * @param voro_vertex_xy X and Y coordinates of the Voronoi vertices.
 */
void mesh_alg_2d_cvd::adaptive_compute_voro_centroid(int pid,double &cx,double &cy,
    int voro_vertex_ct, std::vector<double> &voro_vertex_xy){

    int pid2=2*pid;
    double x=pm2d->xy_id[pid2];
    double y=pm2d->xy_id[pid2+1];
    cx=x;cy=y;

    if(voro_vertex_ct>0){ //the particle has valid Voronoi cell

        //control variable
        bool continue_subdivide=false; //continue to subdivide or not

        //compute c0 directly from the voronoi cell
        //and compute c1 from the particle-vertices-triagle-subdivision of the voronoi cell
        //and store the first level subdivision triangles in tria_vertices_x(y) arrays
        double c0x,c0y,c1x,c1y;

        std::vector<double> tria_vertices_x;
        std::vector<double> tria_vertices_y;
        int tria_make_ct=0;
        for(int triai=0;triai<voro_vertex_ct;triai++){
            double x1_=voro_vertex_xy[2*triai];
            double y1_=voro_vertex_xy[2*triai+1];
            double x2_=voro_vertex_xy[2*(triai+1)];
            double y2_=voro_vertex_xy[2*(triai+1)+1];
            double x3_=x;
            double y3_=y;
            tria_vertices_x.push_back(x1_);
            tria_vertices_x.push_back(x2_);
            tria_vertices_x.push_back(x3_);
            tria_vertices_y.push_back(y1_);
            tria_vertices_y.push_back(y2_);
            tria_vertices_y.push_back(y3_);
            tria_make_ct++;
        }
        double voro_area;
        calculate_voro_centroid_from_triangles(tria_make_ct,tria_vertices_x,tria_vertices_y,c1x,c1y,voro_area);
        
        if(voro_adaptive_subdivide){

            double voro_subdivide_etol_eps=1.0;
            if(tria_iter_ct>1){
                voro_subdivide_etol_eps=prev_pt_mvmt_chrtrt_fac[pid];
                if(voro_subdivide_etol_eps>1.0){voro_subdivide_etol_eps=1.0;}
                if(voro_subdivide_etol_eps<voro_subdivide_etol_eps_min){voro_subdivide_etol_eps=voro_subdivide_etol_eps_min;}
            }
            double subdivide_etol=sqrt(voro_area)*voro_subdivide_etol_eps; //maybe eps ~ prev_pt_mvmt/chrtrt_i, can bounded by eps_min
            
            double rho_voro=pm2d->size_field->getDensityVal(x,y);
            double fac=voro_area*rho_voro;
            double integrand_rho=fac;
            double integrand_xrho=fac*x;
            double integrand_yrho=fac*y;
            double inv_integrand_rho=1.0/integrand_rho;
            c0x=integrand_xrho*inv_integrand_rho;
            c0y=integrand_yrho*inv_integrand_rho;

            //Compute the error ||c0-c1|| against the error tolerance subdivide_etol
            double centroid_error=d_points(c0x,c0y,c1x,c1y);
            if(centroid_error>subdivide_etol){continue_subdivide=true;}
            c0x=c1x;
            c0y=c1y;

            //if continue_subdivide==true, subdivide and compute the new level's voro centroid
            int subdivide_level=1;
            while(continue_subdivide){
                subdivide_level++;

                //loop through the list of triangles, and subdivide each;
                //and obtain the new list of triangles.
                int tria_make_ct_prev=tria_make_ct;
                for(int ind_current=0;ind_current<tria_make_ct_prev;ind_current++){
                    int ind_current_3=3*ind_current;
                    int ind_current_3_1=ind_current_3+1;
                    int ind_current_3_2=ind_current_3+2;
                    double x1_=tria_vertices_x[ind_current_3];
                    double x2_=tria_vertices_x[ind_current_3_1];
                    double x3_=tria_vertices_x[ind_current_3_2]; //generator x
                    double y1_=tria_vertices_y[ind_current_3];
                    double y2_=tria_vertices_y[ind_current_3_1];
                    double y3_=tria_vertices_y[ind_current_3_2]; //generator y

                    //find the longest edge and partisian into two triangles at its midpoint
                    double len12=d_points(x1_,y1_,x2_,y2_);
                    double len13=d_points(x1_,y1_,x3_,y3_); //generator-vertice-len
                    double len23=d_points(x2_,y2_,x3_,y3_); //generator-vertice-len

                    double Ax1, Ax2, Ax3, Ay1, Ay2, Ay3;
                    double Bx1, Bx2, Bx3, By1, By2, By3;

                    if((len12>=len13) && (len12>=len23)){
                        double mid12x=0.5*(x1_+x2_);
                        double mid12y=0.5*(y1_+y2_);
                        Ax1=x1_; Ax2=x3_; Ax3=mid12x;
                        Ay1=y1_; Ay2=y3_; Ay3=mid12y;
                        Bx1=x2_; Bx2=x3_; Bx3=mid12x;
                        By1=y2_; By2=y3_; By3=mid12y;
                    }
                    else if((len13>=len12) && (len13>=len23)){
                        double mid13x=0.5*(x1_+x3_);
                        double mid13y=0.5*(y1_+y3_);
                        Ax1=x1_; Ax2=x2_; Ax3=mid13x;
                        Ay1=y1_; Ay2=y2_; Ay3=mid13y;
                        Bx1=x2_; Bx2=x3_; Bx3=mid13x;
                        By1=y2_; By2=y3_; By3=mid13y;
                    }else{ //len23 is max
                        double mid23x=0.5*(x2_+x3_);
                        double mid23y=0.5*(y2_+y3_);
                        Ax1=x1_; Ax2=x2_; Ax3=mid23x;
                        Ay1=y1_; Ay2=y2_; Ay3=mid23y;
                        Bx1=x1_; Bx2=x3_; Bx3=mid23x;
                        By1=y1_; By2=y3_; By3=mid23y;
                    }

                    tria_vertices_x[ind_current_3]=Ax1;
                    tria_vertices_x[ind_current_3_1]=Ax2;
                    tria_vertices_x[ind_current_3_2]=Ax3;
                    tria_vertices_y[ind_current_3]=Ay1;
                    tria_vertices_y[ind_current_3_1]=Ay2;
                    tria_vertices_y[ind_current_3_2]=Ay3;

                    tria_vertices_x.push_back(Bx1);
                    tria_vertices_x.push_back(Bx2);
                    tria_vertices_x.push_back(Bx3);
                    tria_vertices_y.push_back(By1);
                    tria_vertices_y.push_back(By2);
                    tria_vertices_y.push_back(By3);

                    tria_make_ct++;
                }
                calculate_voro_centroid_from_triangles(tria_make_ct,tria_vertices_x,tria_vertices_y,c1x,c1y);

                //Compute the error ||c0-c1|| against the error tolerance subdivide_etol
                centroid_error=d_points(c0x,c0y,c1x,c1y);
                if(centroid_error<=subdivide_etol){continue_subdivide=false;}
                c0x=c1x;
                c0y=c1y;
            }
        }

        //update centroid of current cell
        cx=c1x;
        cy=c1y;

//code testing: 
//print out: subdivide_level, tria_make_ct, voro area, centroid_error,  subdivide_etol
//char bug00[256];
//sprintf(bug00,"adaptive_voro_tria_subdivide_test_%d.txt",tria_iter_ct);
//FILE *outFile00 = fopen(bug00, "a");
//fprintf(outFile00,"%d %d %g %g %g \n",
//    subdivide_level, tria_make_ct, voro_area, centroid_error, subdivide_etol);
//fclose(outFile00);

    }
}


/**
 * @brief CVD meshing: Updates the particle positions using the Voronoi centroids and determines reTria.
 * 
 * This function updates the positions of the particles based on the Voronoi centroids, adjusted with point movement rules
 * and projection. The updated positions are stored in the `xy_id` array of the `pm2d` object. The function also sets the
 * `reTria` flag to `true` since the positions have been updated, indicating that retriangulation is required.
 */
void mesh_alg_2d_cvd::update_pt_position_and_determine_reTria(){
    //pt position was already updated in Voro compute part. 
    //update them in xy_id array
    #pragma omp parallel for num_threads(num_t)
    for(int i=0;i<Ncurrent;i++){
        int i2=2*i;
        pm2d->xy_id[i2]=xy_id_new[i2];
        pm2d->xy_id[i2+1]=xy_id_new[i2+1];
    }

    //determine reTria: CVD always true
    reTria=true;
} 


/**
 * @brief Update the position of a point using CVD (Centroidal Voronoi Diagram) calculations.
 *
 * This function updates the position of a point based on its Voronoi cell information. It calculates the centroids
 * of the Voronoi cell using the provided point's coordinates, applies point treatments and projections, and obtains
 * the final point position. Additionally, it counts the number of inner points that exceed the stop_Continue movement
 * criteria. It also updates reTria=true;
 *
 * @param i The index of the point to update.
 * @param c The Voronoi cell of the point.
 * @param x The current x-coordinate of the point.
 * @param y The current y-coordinate of the point.
 * @param inner_pt_over_stop_Continue_mvmt_thres_ct_ Reference to a variable to store the count of inner points that exceed the movement criteria.
 * @param use_random_point If using CVD meshing: Decide whether to move point to the centroid of the needle type triangle it connects to instead, to eliminate the needle triangle.
 * @param random_point_x The x-coordinate of the needle type triangle centroid the point should move to.
 * @param random_point_y The y-coordinate of the needle type triangle centroid the point should move to.
 */
void mesh_alg_2d_cvd::update_pt_position_cvd(int i, voronoicell_neighbor_2d &c,
    double x, double y, int &inner_pt_over_stop_Continue_mvmt_thres_ct_,
    bool use_random_point, double random_point_x, double random_point_y){

    double px_new,py_new;
    double final_x,final_y;
    if(use_random_point){
        final_x=random_point_x;
        final_y=random_point_y;
        if(pm2d->pt_ctgr[i]!=-1){
            pm2d->pt_ctgr[i]=-1;
            #pragma omp atomic
            pm2d->inner_pt_ct++;
        }
    }
    else{
        //For each point and its Voronoi cell,
        //calculate the centroids of the Voronoi cell
        //then apply pt treatment and projections,
        //and get final pt positions, store in xy_id
        //also, count number of inner pts over stop_Continue movement criteria
        int voro_vertex_ct=0;
        std::vector<double> voro_vertex_xy;
        //collect voronoi cell vertex
        if(c.p!=0){
            int k=0;
            do {
                voro_vertex_xy.push_back(x+0.5*c.pts[2*k]);   //i<<2 is the same as i*4
                voro_vertex_xy.push_back(y+0.5*c.pts[2*k+1]);
                voro_vertex_ct++;
                k=c.ed[2*k];
            } while (k!=0);
            voro_vertex_xy.push_back(x+0.5*c.pts[0]); 
            voro_vertex_xy.push_back(y+0.5*c.pts[1]); //closed polygon, last=first vertex
        }

        //calculate centroid of Voronoi cell of the current point
        adaptive_compute_voro_centroid(i,px_new,py_new,voro_vertex_ct,voro_vertex_xy);
        pt_mvmt_treatment_projection(i,x,y,px_new,py_new,final_x,final_y);
    }

    double pt_movement=d_points(x,y,final_x,final_y);
    int i_old=(x-ax)*inv_gdx;
    int j_old=(y-ay)*inv_gdy;
    int ij_old=j_old*gnx+i_old;
    prev_pt_mvmt_chrtrt_fac[i]=std::round(pt_movement/chrtrt_len_h(ij_old) * 1000.0) * fac_thousandth;
    //test if inner point pi movement larger than stop_Continue_mvmt_thres
    if(inner_pt_movement_over_stop_Continue_mvmt_thres(i,pt_movement,stop_Continue_mvmt_thres[ij_old])==true){
        inner_pt_over_stop_Continue_mvmt_thres_ct_+=1;
    }

    int i2=2*i;
    xy_id_new[i2]=final_x;
    xy_id_new[i2+1]=final_y;
}



/**
 * @brief Obtain information about Voronoi cells' vertices for use in CVD (Centroidal Voronoi Diagram) meshing calculations.
 *
 * This function computes and extracts information about the vertices of Voronoi cells to be used in triangle subdivision for Voronoi centroid calculation. 
 * It determines whether to store various types of information, such as bar IDs, triangle qualities, triangle vertices, triangle centroids,
 * cell coordinate to cell size ratios, and maximum triangle areas. It also enables the update of point positions using
 * the `update_pt_position_cvd` function during Voronoi computation for each particle. Additionally, it provides an option
 * to store the coordinates and IDs of the cells.
 *
 * The function determines the information to store based on the `printOutputs` flag and performs the necessary computations.
 * Afterward, it outputs the computed information if `printOutputs` is set to true.
 */
void mesh_alg_2d_cvd::voro_compute_retria_extract_info(){

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
}


/**
 * @brief Constructs a `cvd_tria_info` object representing a subdivision triangle of a Voronoi cell and computes the density-triangle-area-weighted centroid.
 *
 * This constructor initializes a `cvd_tria_info` object with the given vertex coordinates and a pointer to the `parallel_meshing_2d` object.
 * It represents a subdivision triangle of a Voronoi cell. The constructor also computes the density-triangle-area-weighted centroid of the triangle using the provided density information.
 *
 * @param x1_ The x-coordinate of the first vertex.
 * @param y1_ The y-coordinate of the first vertex.
 * @param x2_ The x-coordinate of the second vertex.
 * @param y2_ The y-coordinate of the second vertex.
 * @param x3_ The x-coordinate of the third vertex.
 * @param y3_ The y-coordinate of the third vertex.
 * @param pm2d_ A pointer to the `parallel_meshing_2d` object.
 */
cvd_tria_info::cvd_tria_info(double x1_,double y1_,double x2_,double y2_,double x3_,double y3_,parallel_meshing_2d *pm2d_):
pm2d(pm2d_),_x1(x1_),_y1(y1_),_x2(x2_),_y2(y2_),_x3(x3_),_y3(y3_)
{
    compute_info();
}

/**
 * @brief Calculates the density at a given point (x, y) using the size field of the associated `parallel_meshing_2d` object.
 *
 * This function computes the density at the specified point (x, y) by calling the `getDensityVal` function of the size field in the associated `parallel_meshing_2d` object.
 *
 * @param x The x-coordinate of the point.
 * @param y The y-coordinate of the point.
 * @return The density at the given point.
 */
double cvd_tria_info::density(double x, double y){
    return pm2d->size_field->getDensityVal(x,y);
}


/**
 * @brief Computes the information of the triangle.
 *
 * This function calculates various properties of the triangle, including its area, centroid, and integrands for density calculation.
 * The area is computed using the formula for the area of a triangle.
 * The centroid is determined as the average of the coordinates of the triangle's vertices.
 * The density at the centroid is obtained using the `density` function of the associated `parallel_meshing_2d` object.
 * The integrands for density calculation are then computed by multiplying the area with the density at the centroid and its coordinates.
 */
void cvd_tria_info::compute_info(){
    _area=fabs(_x1*(_y2-_y3)+_x2*(_y3-_y1)+_x3*(_y1-_y2))*0.5;

    double one_third=1.0/3.0;
    double centroid_x=(_x1+_x2+_x3)*one_third;
    double centroid_y=(_y1+_y2+_y3)*one_third;
    double rho_centroid=density(centroid_x, centroid_y);
    
    double fac=_area*rho_centroid;
    
    _integrand_rho=fac;
    _integrand_xrho=fac*centroid_x;
    _integrand_yrho=fac*centroid_y;
}














