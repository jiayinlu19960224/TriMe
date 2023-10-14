#include "parallel_meshing_2d.hh"

using namespace voro;

/**
 * @brief Constructor for the parallel_meshing_2d class.
 * Initializes the parallel_meshing_2d object with the provided parameters.
 *
 * @param con_ Pointer to the container_2d object.
 * @param shp_ Pointer to the shape_2d object.
 * @param size_field_ Pointer to the sizing_2d object.
 * @param num_t_ Number of threads.
 * @param output_interval_ Interval for output.
 * @param file_name_prefix_ Prefix for file names.
 */
parallel_meshing_2d::parallel_meshing_2d(container_2d* con_, shape_2d* shp_, sizing_2d* size_field_, int num_t_, int output_interval_, const char *file_name_prefix_)
	: con(con_), shp(shp_), size_field(size_field_), 

	file_name_prefix(file_name_prefix_),
	num_t(num_t_), 

	ax(shp_->ax), bx(shp_->bx), ay(shp_->ay), by(shp_->by),
	lx(bx-ax), ly(by-ay),

	gnx(shp_->gnx),gny(shp_->gny),gnxy(shp_->gnxy),
	gdx(shp_->gdx),gdy(shp_->gdy),diag_gdxy(shp_->diag_gdxy),
	inv_gdx(shp_->inv_gdx),inv_gdy(shp_->inv_gdy), gdxy(shp_->gdx * shp_->gdy),
	geo_bgrid_ct(shp_->geo_bgrid_ct),geo_igrid_ct(shp_->geo_igrid_ct),

	xy_id(new double[1]), pt_ctgr(new int[1]),chrtrt_len_h(new double[1]),
	bgrid_geps(new double[1]), bgrid_deps(new double[1]),bdry_adf_construction(false),
	Ntotal(0),Ncurrent(0),Nremain(0),Nfixed(0),inner_pt_ct(0),

	output_interval(output_interval_)
	{
	}


/**
 * @brief Destructor for the parallel_meshing_2d class.
 * Cleans up dynamically allocated memory.
 */
parallel_meshing_2d::~parallel_meshing_2d(){

	delete [] xy_id;
	delete [] pt_ctgr;
	delete [] chrtrt_len_h;
	delete [] bgrid_geps;
	delete [] bgrid_deps;

	if(bdry_adf_construction){
		#pragma omp parallel for num_threads(num_t)
		for(int i=0;i<geo_bgrid_ct;i++){
			delete bgrid_adf[i];
			delete bgrid_adf_stat[i];
		}
		delete [] bgrid_adf;
		delete [] bgrid_adf_stat;
	}
	
}


/**
 * Returns the geometry tolerance (geps) associated with the boundary grid that the given coordinates lie in.
 *
 * @param x The x-coordinate of the point.
 * @param y The y-coordinate of the point.
 * @return The geometry tolerance (geps) associated with the boundary grid the point lies in.
 * @throws std::exception if the point is not within the geometry boundary grid.
 */
double parallel_meshing_2d::get_bgrid_geps(double x, double y){
	int ind=geo_grid(x,y)-1;
	if(ind<0 || ind>=gnxy){
		printf("ERROR 1: pt not in geometry bdry grid.\n");
    	throw std::exception();
    }
	return bgrid_geps[ind];
}

/**
 * Calculates the geometry tolerance (geps) function for the given coordinates.
 *
 * @param x The x-coordinate of the point.
 * @param y The y-coordinate of the point.
 * @param geps_prime The default value of geps.
 * @return The calculated value of geps for the point.
 */
double parallel_meshing_2d::geps_func(double x, double y, double geps_prime){
	int ind=geo_grid(x,y)-1;
	if(ind<0 || ind>=gnxy){ //pt not in geometry bdry grid.
		return geps_prime; //dummy
    }
	return bgrid_geps[ind];
}

/**
 * Calculates the spatial step size (deps) value for the given coordinates.
 *
 * @param x The x-coordinate of the point.
 * @param y The y-coordinate of the point.
 * @param deps_prime The default value of deps.
 * @return The calculated value of deps for the point.
 */
double parallel_meshing_2d::deps_func(double x, double y, double deps_prime){
	int ind=geo_grid(x,y)-1;
	if(ind<0 || ind>=gnxy){ //pt not in geometry bdry grid.
		return deps_prime; //dummy
    }
	return bgrid_deps[ind];
}



/**
 * Calculates the normalized density arrays and the characteristic edge length (chrtrt_len_h) for all inner and boundary grids.
 *
 * @param bgrid_ingeo_ct Number of boundary grids in/on the geometry.
 * @param bgrid_outgeo_ct Number of boundary grids outside the geometry.
 * @param bgrid_ingeo_ind Vector containing the indices of the boundary grids in/on the geometry.
 * @param bgrid_outgeo_ind Vector containing the indices of the boundary grids outside the geometry.
 * @param density_igrid Array to store the normalized densities of the inner grids.
 * @param density_bgrid_ingeo Array to store the normalized densities of the boundary grids in/on the geometry.
 * @param Ncurrent_rho_igrid Array to store the expected number of points in the inner grids.
 * @param Ncurrent_rho_bgrid_ingeo Array to store the expected number of points in the boundary grids in/on the geometry.
 */
void parallel_meshing_2d::get_normalized_density_chrtrt_h(
	int bgrid_ingeo_ct,int bgrid_outgeo_ct,
	std::vector<int> &bgrid_ingeo_ind,std::vector<int> &bgrid_outgeo_ind, 
    double *&density_igrid, double *&density_bgrid_ingeo,
	double *&Ncurrent_rho_igrid,double *&Ncurrent_rho_bgrid_ingeo)
{
	//sum of density from S=S1+S2
	double density_sum=0.0;
	//obtain density_sum for density normalization later
	#pragma omp parallel for reduction(+:density_sum) num_threads(num_t) 
	for(int gi=0;gi<geo_igrid_ct;gi++){
		//get inner grid ij, and its midpoint, and calcuate its density
		int ij=geo_igrid_ij(gi);
		int i=ij%gnx;
        int j=ij/gnx;
        double x=ax+(0.5+i)*gdx; 
        double y=ay+(0.5+j)*gdy;
        double rho_ij=density(x,y); 
        density_igrid[gi]=rho_ij;
        density_sum+=rho_ij;
	}

	#pragma omp parallel for reduction(+:density_sum) num_threads(num_t) 
	for(int gii=0;gii<bgrid_ingeo_ct;gii++){
		//get inner grid ij, and its midpoint, and calcuate its density
		int gi=bgrid_ingeo_ind[gii];
		int ij=geo_bgrid_ij(gi);
		int i=ij%gnx;
        int j=ij/gnx;
        double x=ax+(0.5+i)*gdx; 
        double y=ay+(0.5+j)*gdy;
    	double rho_ij=density(x,y); 
        density_bgrid_ingeo[gii]=rho_ij;
        density_sum+=rho_ij;	
	}

	chrtrt_len_h_avg=0;
	//Calculate the normalized density arrays, and the chrtrt_len_h for all inner and bdry grids
	#pragma omp parallel num_threads(num_t)
	{
		//Firstly, inner grids
		#pragma omp for reduction(+:chrtrt_len_h_avg) 
		for(int gi=0;gi<geo_igrid_ct;gi++){
	        density_igrid[gi]=density_igrid[gi]/density_sum;
	        //calculate the expected number of points in the inner grid: Ncurrent*rho_normalized
	        Ncurrent_rho_igrid[gi]=density_igrid[gi]*Ncurrent;
	        //calcualte characteristic element size lengthscale h
	        int ij=geo_igrid_ij(gi);
	        chrtrt_len_h[ij]=sqrt((gdxy)/Ncurrent_rho_igrid[gi]);
	        chrtrt_len_h_avg+=chrtrt_len_h[ij];
		}

		//Next, ingeo bdry grids
		#pragma omp for reduction(+:chrtrt_len_h_avg)
		for(int gii=0;gii<bgrid_ingeo_ct;gii++){
			int gi=bgrid_ingeo_ind[gii];
			int ij=geo_bgrid_ij(gi);

	        density_bgrid_ingeo[gii]=density_bgrid_ingeo[gii]/density_sum;
	        //calculate the expected number of points in the bdry grid: Ncurrent*rho_normalized
	        Ncurrent_rho_bgrid_ingeo[gii]=density_bgrid_ingeo[gii]*Ncurrent;
	        //calcualte characteristic element size lengthscale h
	        chrtrt_len_h[ij]=sqrt((gdxy)/Ncurrent_rho_bgrid_ingeo[gii]);
	        chrtrt_len_h_avg+=chrtrt_len_h[ij];
		}

		//Lastly, outgeo bdry grids
		#pragma omp for
		for(int gii=0;gii<bgrid_outgeo_ct;gii++){
			int gi=bgrid_outgeo_ind[gii];
			int ij=geo_bgrid_ij(gi);
			int i=ij%gnx;
	    	int j=ij/gnx;

	    	//Use the same procedure to calculate chrtrt_len_h for bdry grids (using the same normalizing constancts), 
	    	//so that bdry grid's hi and inner grid's hi satisfy their relative density relationship:
	    	//hi/hj=sqrt(rhoj/rhoi)
	    	//this approach is good for capturing cases when bdry grid hi is much smaller than nearby inner grid's hi (i.e. narrow elongated band/corners)
	        double x=ax+(0.5+i)*gdx; 
	        double y=ay+(0.5+j)*gdy;
    		double density_bgrid_nor=density(x,y)/density_sum;
		    //calcualte characteristic element size lengthscale h
		    chrtrt_len_h[ij]=sqrt((gdxy)/(density_bgrid_nor*Ncurrent));
		}
	}
	chrtrt_len_h_avg=chrtrt_len_h_avg/(geo_igrid_ct+bgrid_ingeo_ct);
}


/**
 * Get adaptive arrays for boundary grids.
 *
 * @param scale_adf_error_tol Scaling factor for the Adaptive Distance Field (ADF) quadtree error tolerance.
 */
void parallel_meshing_2d::get_bdry_grid_adaptive_arrays(double scale_adf_error_tol)
{
	#pragma omp parallel for num_threads(num_t) schedule(dynamic)
	for(int gi=0;gi<geo_bgrid_ct;gi++){
		int ij=geo_bgrid_ij(gi);
		int i=ij%gnx;
    	int j=ij/gnx;
    	double grid_ax=ax+(i)*gdx; 
    	double grid_bx=ax+(1+i)*gdx; 
    	double grid_ay=ay+(j)*gdy;
    	double grid_by=ay+(1+j)*gdy;

        //calculate geometric adaptive arrays based on bdry grid hi
		//conservatively restricted to smaller value of hi or hi_avg
		bgrid_geps[gi]=geps_h_frac*min(chrtrt_len_h[ij],chrtrt_len_h_avg); //geometry boundry tolerance
		bgrid_deps[gi]=sqrt(eps)*min(chrtrt_len_h[ij],chrtrt_len_h_avg); //deltaX,deltaY used in projection() 

		//construct ADF quad tree cells for the bdry grids, adaptive to local geometry bdry tolerance
		int current_level=0;
		double adf_err_tol=adf_err_geps_frac*bgrid_geps[gi]*scale_adf_error_tol;
		bgrid_adf_stat[gi]=new adf_stat_2d();
		bgrid_adf[gi]=new adf_2d(bgrid_adf_stat[gi],shp,
			grid_ax,grid_bx,grid_ay,grid_by,
			adf_err_tol,adf_max_depth,current_level);
	}
}


/**
 * Collects valid neighbors for error diffusion using Floyd-Steinberg dithering.
 *
 * @param nei_ij The grid index of the neighbor.
 * @param nei_ct Reference to the count of valid neighbors.
 * @param inner_nei_gi Vector to store the indices of inner grid neighbors.
 * @param ingeo_bdry_nei_gii Vector to store the indices of boundary grid neighbors inside the geometry.
 * @param bgrid_ingeo_outgeo_ctgr Array indicating the category of each boundary grid (ingeo_bdry_grid or outgeo_bdry_grid).
 */
void parallel_meshing_2d::error_diffusion_collect_nei(
	int nei_ij, int &nei_ct,
	std::vector<int> &inner_nei_gi,
	std::vector<int> &ingeo_bdry_nei_gii,
	int *&bgrid_ingeo_outgeo_ctgr)
{

	if(nei_ij>=0 && nei_ij<gnxy){
		int grid_ctgr=geo_grid(nei_ij);	
		int gi=-1;

		//bdry grid
		if(grid_ctgr>0 && grid_ctgr<=gnxy){ 
			gi=grid_ctgr-1;
    		if(bgrid_ingeo_outgeo_ctgr[gi]>0){ //a ingeo_bdry_grid
        		int gii=bgrid_ingeo_outgeo_ctgr[gi]-1;
        		ingeo_bdry_nei_gii.push_back(gii);
        		nei_ct++;
    		}
		}

		//inner grid
		if(grid_ctgr<0){ 
			gi=abs(grid_ctgr)-1;
			inner_nei_gi.push_back(gi);
			nei_ct++;
		}
	}
}


/**
 * Checks the category of each point and prints an error message if the category is incorrect.
 * This function is used for debugging purposes.
 */
void parallel_meshing_2d::check_pt_ctgr(){
	int inner_pt_ct_temp=0;
	for(int i=0;i<Ncurrent;i++){
		int i2=2*i;
		double x=xy_id[i2];
		double y=xy_id[i2+1];
		int g_ind0=geo_grid(x,y);
		if(g_ind0<0){
			inner_pt_ct_temp++;
			if(pt_ctgr[i]!=-1){
				printf("pt %d (%g %g) in I,I pt but ctgr wrong \n",i,x,y);
			}
		}
		else if(g_ind0>0 && g_ind0<=gnxy){ 
			double new_xy_adf=geo_bgrid_adf(x,y);
			double new_xy_bgrid_geps=get_bgrid_geps(x,y);
    		if(abs(new_xy_adf)<=new_xy_bgrid_geps){ //bdry pt
    			if(pt_ctgr[i]!=1){
					printf("pt %d (%g %g) in B,B pt but ctgr wrong \n",i,x,y);
				}
    		}
    		else if(new_xy_adf<-new_xy_bgrid_geps){
    			inner_pt_ct_temp++;
    			if(pt_ctgr[i]!=-1){
					printf("pt %d (%g %g) in B,I pt but ctgr wrong \n",i,x,y);
				}
    		}
    		else{
    			printf("ERROR: check pt i (%g %g) in geo bdry grid but not in/on geometry.\n",i,x,y);
    		}
				    		
		}
		else{
			printf("ERROR: check pt i (%g %g) in outer grid.\n",i,x,y);
		}
	}
	if(inner_pt_ct_temp!=inner_pt_ct){
		printf("check: inner_pt_ct_temp: %d inner_pt_ct: %d wrong\n",inner_pt_ct_temp,inner_pt_ct);
	}
}


/**
 * Generates Npt points in the specified grid (i, j).
 *
 * @param generate_pt Flag indicating whether to generate points.
 * @param Npt Number of points to generate.
 * @param i The x-coordinate of the grid.
 * @param j The y-coordinate of the grid.
 * @param ij The index of the grid.
 * @param gi The index of the boundary/inner grid in their associated arrays.
 * @param grid_ctgr The category of the grid.
 * @param pt_counter Reference to the point counter.
 * @param xy_id_temp Pointer to the array to store the generated points.
 * @param pt_ctgr_temp Pointer to the array to store the categories of the generated points.
 * @param seed Reference to the random seed.
 * @param inner_pt_ct_temp Reference to the counter for inner points.
 */
void parallel_meshing_2d::generate_pt_in_grid(bool generate_pt,int Npt,
    int i, int j, int ij, int gi, int grid_ctgr,int &pt_counter,
    double *&xy_id_temp,int *&pt_ctgr_temp,
    unsigned int &seed, int &inner_pt_ct_temp)
{
    if(generate_pt && Npt>0){

        double grid_ax=ax+(i)*gdx; 
        double grid_ay=ay+(j)*gdy;
        while(Npt>0){
        	Npt--;
            double x=grid_ax+rnd_r(&seed)*gdx;
            double y=grid_ay+rnd_r(&seed)*gdy;
            int m;
            #pragma omp atomic capture
            m=pt_counter++;

            pt_ctgr_temp[m]=0;

            if(grid_ctgr<0){ //inner grid, store directly as inner pt
                
                xy_id_temp[2*m]=x;
                xy_id_temp[2*m+1]=y;
                pt_ctgr_temp[m]=-1; //-1: inner pt
                inner_pt_ct_temp+=1;

            }
            else{ //ingeo_bdry grid, need to project pt if pt outside; inner/bdry pt
                double finalx=x; double finaly=y;
                double pt_sdf_temp=sdf_adf(finalx,finaly);
                double final_xy_bgrid_geps=bgrid_geps[gi];
                bool find_new_pt=false;

                if(pt_sdf_temp>final_xy_bgrid_geps){
                	find_new_pt=true;
                }

                while(find_new_pt){

                    bool project_pt_outside_geo=true;
                    if(pt_projection(finalx, finaly,x,y,bgrid_deps[gi],bgrid_geps[gi],project_pt_outside_geo)==false){

                        x=grid_ax+rnd_r(&seed)*gdx;
                        y=grid_ay+rnd_r(&seed)*gdy;
                        finalx=x; finaly=y;
                        pt_sdf_temp=sdf_adf(finalx,finaly);
                        if(pt_sdf_temp<=final_xy_bgrid_geps){
		                	find_new_pt=false;
		                }
                    }
                    else{
                    	find_new_pt=false;
                    	pt_sdf_temp=sdf_adf(finalx,finaly);
                    }
                }     
  
                xy_id_temp[2*m]=finalx;
                xy_id_temp[2*m+1]=finaly;

                int final_xy_g_ind=geo_grid(finalx,finaly);
                if(final_xy_g_ind>0 && final_xy_g_ind<=gnxy){ //bdry grid
					final_xy_bgrid_geps=get_bgrid_geps(finalx,finaly);
					if(fabs(pt_sdf_temp)<=final_xy_bgrid_geps){ //bdry pt
						pt_ctgr_temp[m]=1; //1: bdry pt
					}
					else if(pt_sdf_temp<-final_xy_bgrid_geps){ //inner pt
						pt_ctgr_temp[m]=-1; //-1: inner pt
                    	inner_pt_ct_temp+=1;
					}
					else{ //outer pt in bdry grid
						printf("ERROR: pt in pt_init generated in ingeo_bdry_grid pt outside geometry. Should not happen.\n");
						throw std::exception();
					}
				}
				else if(final_xy_g_ind<0){ //inner grid
					pt_ctgr_temp[m]=-1; //-1: inner pt
                    inner_pt_ct_temp+=1;

				}
				else{ //outer grid: should not happen, since projection() return false if pt in outer grid
					printf("ERROR: pt in pt_init generated in ingeo_bdry_grid pt outside geometry. Should not happen.\n");
					throw std::exception();
				}
            }
        }
    }
}


/**
 * Performs error diffusion by diffusing the error to the surrounding neighbors.
 *
 * @param generate_pt Flag indicating whether to generate points.
 * @param diff The amount of error to be diffused.
 * @param ij The grid index.
 * @param bgrid_ingeo_outgeo_ctgr Array indicating the category of each boundary grid (ingeo_bdry_grid or outgeo_bdry_grid).
 * @param Ncurrent_rho_igrid Array storing the expected number of points in the inner grids.
 * @param Ncurrent_rho_bgrid_ingeo Array storing the expected number of points in the boundary grids in/on the geometry.
 */
void parallel_meshing_2d::error_diffusion(bool generate_pt, double diff,
	int ij, int *&bgrid_ingeo_outgeo_ctgr,
	double *&Ncurrent_rho_igrid,double *&Ncurrent_rho_bgrid_ingeo)
{

	if(generate_pt && diff!=0){
		int nei_ct=0;
		std::vector<int> inner_nei_gi;
		std::vector<int> ingeo_bdry_nei_gii;

		int nei_r_ij=ij+1; error_diffusion_collect_nei(nei_r_ij,nei_ct,inner_nei_gi,ingeo_bdry_nei_gii,bgrid_ingeo_outgeo_ctgr);
		int nei_u_ij=ij+gnx; error_diffusion_collect_nei(nei_u_ij,nei_ct,inner_nei_gi,ingeo_bdry_nei_gii,bgrid_ingeo_outgeo_ctgr);
		int nei_ur_ij=nei_u_ij+1; error_diffusion_collect_nei(nei_ur_ij,nei_ct,inner_nei_gi,ingeo_bdry_nei_gii,bgrid_ingeo_outgeo_ctgr);
		int nei_ul_ij=nei_u_ij-1; error_diffusion_collect_nei(nei_ul_ij,nei_ct,inner_nei_gi,ingeo_bdry_nei_gii,bgrid_ingeo_outgeo_ctgr);

		//diffuse error to neighbors
		if(nei_ct>0){
			//equal weight in error diffusion
			double nei_weight=1.0/(1.0*nei_ct); 
			double error_diffuse_per_nei=diff*nei_weight;
			//diffuse to inner grids
			for(int nei_i=0;nei_i<(signed int)inner_nei_gi.size();nei_i++){
				int nei_gi=inner_nei_gi[nei_i];
				Ncurrent_rho_igrid[nei_gi]+=error_diffuse_per_nei;
			}
			//diffuse to ingeo_bdry grids
			for(int nei_i=0;nei_i<(signed int)ingeo_bdry_nei_gii.size();nei_i++){
				int nei_gii=ingeo_bdry_nei_gii[nei_i];
				Ncurrent_rho_bgrid_ingeo[nei_gii]+=error_diffuse_per_nei;
			}
		}
	}
}



/** 
 * Add fixed points into the mesh. Rescale the points using the same scaling as the scaling of shp.
 * 
 * @param Nfixed The number of fixed points
 * @param shp A custom_shape_2d model, whose scaling parameters we will use for normalizing the fixed points
 * @param fixed_pt_list the x and y coordinations of the fixed points
 */
void parallel_meshing_2d::add_fixed_points_normailze(int Nfixed_, shape_2d* shp_, double* fixed_pt_list){
	Nfixed=Nfixed_;
	xy_id_fixed_pt_init=new double[2*Nfixed];
	#pragma omp parallel for num_threads(num_t)
    for(int i=0;i<Nfixed;i++){
    	int i2=2*i;
    	double ptx=fixed_pt_list[i2];
    	double pty=fixed_pt_list[i2+1];

    	if(shp_->shape_scaling==true){
			//scale and center x component of point
            ptx=shp_->scale_min_domain_range*normalize_model_length_fac*(ptx-shp_->scale_xmid)/shp_->scale_max_range+0.5*lx;
            //scale and center y component of point
            pty=shp_->scale_min_domain_range*normalize_model_length_fac*(pty-shp_->scale_ymid)/shp_->scale_max_range+0.5*ly;
		}

    	xy_id_fixed_pt_init[i2]=ptx;
    	xy_id_fixed_pt_init[i2+1]=pty;
    }
}
/** 
 * Add fixed points into the mesh.
 * 
 * @param Nfixed The number of fixed points
 * @param fixed_pt_list the x and y coordinations of the fixed points
 */
void parallel_meshing_2d::add_fixed_points(int Nfixed_, double* fixed_pt_list){
	Nfixed=Nfixed_;
	xy_id_fixed_pt_init=new double[2*Nfixed];
	#pragma omp parallel for num_threads(num_t)
    for(int i=0;i<Nfixed;i++){
    	int i2=2*i;
    	xy_id_fixed_pt_init[i2]=fixed_pt_list[i2];
    	xy_id_fixed_pt_init[i2+1]=fixed_pt_list[i2+1];
    }
}



/**
 * Adding fixed points to the mesh. Used in pt_init() function.
 */
void parallel_meshing_2d::add_fixed_points_procedure(){
	if(Nfixed>0){
		//print fixed points
		char bug0[256];
	     sprintf(bug0,"%s/fixed_pt_coords.par",file_name_prefix);
	     FILE *outFile0 = fopen(bug0, "a");
	     for(int i=0;i<Nfixed;i++){
	        fprintf(outFile0,"%g %g \n",xy_id_fixed_pt_init[2*i], xy_id_fixed_pt_init[2*i+1]);
	     }
	     fclose(outFile0);

		printf("start adding %d fixed pt \n",Nfixed);
	//---------------------------------------------------------------------
		//the first Nfixed points are fixed points
	    #pragma omp parallel for num_threads(num_t) schedule(guided) 
	    for(int i=0;i<Nfixed;i++){
	    	int i2=2*i;
	    	xy_id[i2]=xy_id_fixed_pt_init[i2];
	    	xy_id[i2+1]=xy_id_fixed_pt_init[i2+1];
	    }

		Nremain-=Nfixed;
		int Ncurrent_old=Ncurrent;
		Ncurrent+=Nfixed;

	//---------------------------------------------------------------------
	    //update adaptive geometric arrays
	    double adaptive_scale_fac=1.0/(sqrt(1.0*(Ncurrent/Ncurrent_old)));
	    chrtrt_len_h_avg*=adaptive_scale_fac;
		#pragma omp parallel num_threads(num_t)
		{
			#pragma omp for
			for(int gi=0;gi<geo_bgrid_ct;gi++){
				int ij=geo_bgrid_ij(gi);
				chrtrt_len_h[ij]*=adaptive_scale_fac;
				bgrid_geps[gi]*=adaptive_scale_fac;
				bgrid_deps[gi]*=adaptive_scale_fac;
			}
			#pragma omp for
			for(int gi=0;gi<geo_igrid_ct;gi++){
				int ij=geo_igrid_ij(gi);
				chrtrt_len_h[ij]*=adaptive_scale_fac;
			}
		}

		printf("Ncurrrent is now %d \n", Ncurrent);
	//---------------------------------------------------------------------
	    //loop through points against new geps, 
	    //re-ctgr them to pt_ctgr
	    //and project any outside points back to the new geps
	    int inner_pt_ct_temp=0;   
	    #pragma omp parallel for num_threads(num_t) schedule(dynamic) reduction(+: inner_pt_ct_temp)
		for(int i=0;i<Ncurrent;i++){
			int i2=2*i;
			double x=xy_id[i2];
			double y=xy_id[i2+1];
			int g_ind0=geo_grid(x,y);
	    	if(g_ind0<0){ //inner grid
	    		pt_ctgr[i]=-1; //-1: inner pt
	    		inner_pt_ct_temp+=1;
	    	}
	    	else if(g_ind0>0 && g_ind0<=gnxy){  //bdry grid
	    		double pt_adf=geo_bgrid_adf(x,y);
	    		double pt_bgrid_geps=get_bgrid_geps(x,y);
	    		if(abs(pt_adf)<=pt_bgrid_geps){ //bdry pt
	    			pt_ctgr[i]=1;
	    		}
	    		else if(pt_adf<-pt_bgrid_geps){
					pt_ctgr[i]=-1; //-1: inner pt
	    			inner_pt_ct_temp+=1;
				}
				else{ //pt_adf>pt_bgrid_geps: outside of geometry with new geps
					if(i<Nfixed){
						printf("ERROR: Pt init procedure: fixed pt (%g %g) in geo bdry grid but not in/on geometry.\n",x,y);
				    	throw std::exception();
					}

					double new_x, new_y;
					if(pt_projection(new_x,new_y,x,y,-1,-1)==true){ //(new_x,new_y) on geo bdry/inner grid now
						xy_id[i2]=new_x;
						xy_id[i2+1]=new_y;

						//decide if new pt is bdry pt or inner pt
				    	int g_ind=geo_grid(new_x,new_y);
				    	if(g_ind>0 && g_ind<=gnxy){ //bdry grid
				    		double new_xy_adf=geo_bgrid_adf(new_x,new_y);
				    		double new_xy_bgrid_geps=get_bgrid_geps(new_x,new_y);
				    		if(abs(new_xy_adf)<=new_xy_bgrid_geps){ //bdry pt
				    			pt_ctgr[i]=1;
				    		}
				    		else if(new_xy_adf<-new_xy_bgrid_geps){
				    			pt_ctgr[i]=-1; //-1: inner pt
				    			inner_pt_ct_temp+=1;
				    		}
				    		else{
				    			printf("ERROR: add_fixed_pt: projected bdry pt now (%g %g) in geo bdry grid but not in/on geometry.\n",new_x,new_y);
				    			throw std::exception();
				    		}
				    	}
				    	else if(g_ind<0){ //inner grid
				    		pt_ctgr[i]=-1; //-1: inner pt
				    		inner_pt_ct_temp+=1;
				    	}
				    	else{
				    		printf("ERROR: add_fixed_pt: projected bdry pt now (%g %g) not in geometry inner/bdry grid.\n",new_x,new_y);
				    		throw std::exception();
				    	}

					}
					else{ //fail to project pt onto bdry: this should happen very scarsely!
						//loop outwards of grids, and if find a inner grid, put pt on its midpt instead
						//and tag pt inner pt
						printf("bdry pt with new geps projection fail: this should be very rare! \n");
						int incre=0;
						int pgridi=(x-ax)*inv_gdx;
	                	int pgridj=(y-ay)*inv_gdy;
	                	int pgridij=gnx*pgridj+pgridi;
	                	bool continue_search=true;

	                	std::vector<int> inner_grid_found;
            	        int inner_grid_found_ct=0;
	                	
	                	while(continue_search){
	                		int il=pgridi-incre;int ih=pgridi+incre;
	                		int jl=pgridj-incre;int jh=pgridj+incre;
	                		if(il<0&&ih>=gnx&&jl<0&&jh>=gny){
	                			printf("ERROR: searching inner grid but out of bound, shouldn't happen.\n");
	    						throw std::exception();
	                		}
	                		for(int gridi=il;gridi<=ih;gridi++){
	                			for(int gridj=jl;gridj<=jh;gridj++){
	                				//only search the new layer
	                				if(gridi==il||gridi==ih||gridj==jl||gridj==jh){
	                					//only search if grid is valid
	                					if(gridi>=0&&gridi<gnx&&gridj>=0&&gridj<gny){
	                						int gridij=gnx*gridj+gridi;
	                						if(geo_grid(gridij)<0){ //inner grid
	                							continue_search=false;
	                							inner_grid_found.push_back(gridi);
		            							inner_grid_found.push_back(gridj);
		            							inner_grid_found_ct++;
	                						}
	                					}
	                				}
	                			}
	                		}
	                		incre++;
	                	}

	                	unsigned int seed_i=omp_get_thread_num()+100;
	                	//Select a random inner grid in the outer layer found
		            	int random_igrid_found_i2=2* (rand_r(&seed_i) % inner_grid_found_ct); //random int from 0 to inner_grid_found_ct-1
		            	int gridi=inner_grid_found[random_igrid_found_i2];
		            	int gridj=inner_grid_found[random_igrid_found_i2+1];
		            	double random_point_x=ax+(rnd_r(&seed_i)+gridi)*gdx; 
		            	double random_point_y=ay+(rnd_r(&seed_i)+gridj)*gdy;
						xy_id[i2]=random_point_x;
						xy_id[i2+1]=random_point_y;
						pt_ctgr[i]=-1; //tag inner pt
	                	inner_pt_ct_temp+=1;
					}
				}
			}
			else{
				if(i<Nfixed){
					printf("ERROR: Pt init procedure: Fixed pt (%g %g) not in inner nor geo bdry grid.\n",x,y);
	    			throw std::exception();
				}
				else{
					printf("ERROR: Initialized pt %d (%g %g) not in inner nor geo bdry grid.\n",i,x,y);
	    			throw std::exception();
				}
	    	}
			
		}
		inner_pt_ct=inner_pt_ct_temp;

		delete [] xy_id_fixed_pt_init;

		printf("finished adding fixed points \n");
		//---------------------------------------------------------------------
	}

}



/**
 * @brief Initialize the points in the container. It first initialize int(pt_init_frac*Ntotal) points. It then adds the inputted Nfixed fixed points into the mesh.
 * 
 * It initializes various variables and arrays used in the point initialization process.
 * It calculates the normalized density arrays and characteristic length arrays for inner and boundary grids.
 * It constructs adaptive arrays for boundary grids based on the error tolerance.
 * It performs error diffusion and calculates the number of points to generate for each grid.
 * It generates random points within each inner and ingeo_bdry grid according to the density field and performs error diffusion.
 * It stores the initially generated points in temporary arrays and handles cases where too many or too few points were generated.
 * It generates the remaining points needed to reach the desired total number of points.
 * It assigns the generated points to the xy_id and pt_ctgr arrays and adds them to the container.
 * It cleans up temporary arrays and releases memory.
 *
 * @param Ntotal_ Total number of points to generate. The function initializes Ncurrent = int(pt_init_frac * Ntotal) points to start with in the container.
 *  			  It then adds the inputted Nfixed fixed points into the mesh.
 */
void parallel_meshing_2d::pt_init(int Ntotal_){

	//-------------A. Init general set up--------------//

	Ntotal=Ntotal_;
	Ncurrent=min(int(pt_init_frac*Ntotal),Ntotal-Nfixed);
	inner_pt_ct=0;
	Nremain=Ntotal-Ncurrent;
	int pt_counter=0; //after the end of the loop, this is the number of points generated for inner grids.
	//bdry ADF cells should have error tolerance based on geps of Ntotal, 
	//rather than geps of Ncurrent, 
	//since they should only be constructed once and can be used throughout.
	double scale_adf_error_tol=1.0/sqrt(1.0*Ntotal/Ncurrent);

	printf("start pt init routine, geneating %d pts\n",Ncurrent);

	delete [] xy_id;
	xy_id=new double[Ntotal*2];
	delete [] pt_ctgr;
	pt_ctgr=new int[Ntotal];
	delete [] chrtrt_len_h;
	chrtrt_len_h=new double[gnxy];

	//initialize array values to undefined dummy values
	#pragma omp parallel for num_threads(num_t)
	for(int i=0;i<Ntotal;i++){
		xy_id[2*i]=ax+bx; //x=ax+bx: undefined
		xy_id[2*i+1]=ay+by; //y=ay+by: undefined
		pt_ctgr[i]=0; //0: undefined
	}
	#pragma omp parallel for num_threads(num_t)
	for(int i=0;i<gnxy;i++){
		chrtrt_len_h[i]=-1; //-1: undefined
	}

	//For bdry grids: 
	//Compute adaptive geometry related arrays: geps, deps, adf cells
	delete [] bgrid_geps;
	bgrid_geps=new double[geo_bgrid_ct];
	delete [] bgrid_deps;
	bgrid_deps=new double[geo_bgrid_ct];
	bdry_adf_construction=true;
	bgrid_adf=new adf_2d*[geo_bgrid_ct];
	bgrid_adf_stat=new adf_stat_2d*[geo_bgrid_ct];

	std::vector<int> bgrid_ingeo_ind; //index of bdry grids whose center sdf(cx,cy)<=0.5*min(gdx,gdy)
	std::vector<int> bgrid_outgeo_ind; //index of bdry grids whose center sdf(cx,cy)>0.5*min(gdx,gdy)
	int *bgrid_ingeo_outgeo_ctgr=new int[geo_bgrid_ct];
	int bgrid_ingeo_ct=0;
	int bgrid_outgeo_ct=0;
	for(int gi=0;gi<geo_bgrid_ct;gi++){
		int ij=geo_bgrid_ij(gi);
		int i=ij%gnx;
        int j=ij/gnx;
        double x=ax+(0.5+i)*gdx; 
        double y=ay+(0.5+j)*gdy;
        if(sdf(x,y)<=0.5*min(gdx,gdy)){
        	bgrid_ingeo_ind.push_back(gi);
        	bgrid_ingeo_ct++;
        	bgrid_ingeo_outgeo_ctgr[gi]=bgrid_ingeo_ct; //ingeo: >0; 
        }
        else{
        	bgrid_outgeo_ind.push_back(gi);
        	bgrid_outgeo_ct++;
        	bgrid_ingeo_outgeo_ctgr[gi]=-bgrid_outgeo_ct; //outgeo: <0; 
        }
	}

	//-------------B. Calculate the normalized density arrays, and the chrtrt_len_h for all inner and bdry grids--------------//
	//set S1:density array for every inner grid and in_geo bdry grid
	double *density_igrid=new double[geo_igrid_ct];
	//set S2:density array for bdry grids whose sdf(cx,cy)<0
	double *density_bgrid_ingeo=new double[bgrid_ingeo_ct];
	//expected number of points in a inner grid
	double *Ncurrent_rho_igrid=new double[geo_igrid_ct];
	//expected number of points in a bdry grid in set S2
	double *Ncurrent_rho_bgrid_ingeo=new double[bgrid_ingeo_ct];

	get_normalized_density_chrtrt_h(bgrid_ingeo_ct,bgrid_outgeo_ct,
	bgrid_ingeo_ind,bgrid_outgeo_ind,density_igrid,density_bgrid_ingeo,
	Ncurrent_rho_igrid,Ncurrent_rho_bgrid_ingeo);

	//-------------C. get adaptive arrays for bdry grids--------------//
	get_bdry_grid_adaptive_arrays(scale_adf_error_tol);

	//-------------D. generate random points in inner and ingeo_bdry grids, according to density field, and do error diffusion--------------//
	
	//Serial process: Error diffusion and calculate Npt to generate for each grid
    int Npt_total_temp=0;
    int *grid_Npt=new int[gnxy];
    for(int j=0;j<gny;j++){
        for(int i=0;i<gnx;i++){
            int ij=gnx*j+i; 
            int grid_ctgr=geo_grid(ij);
            bool generate_pt=false;
            int Npt=0; //number of points to generate
            double diff=0; //error residual
            int gi=-1;

            if(grid_ctgr>0 && grid_ctgr<=gnxy){ //bdry grid
                gi=grid_ctgr-1;
                if(bgrid_ingeo_outgeo_ctgr[gi]>0){ //a ingeo_bdry_grid
                    generate_pt=true;
                    int gii=bgrid_ingeo_outgeo_ctgr[gi]-1;
                    //generate Npt random points in the bdry grid
                    Npt=round(Ncurrent_rho_bgrid_ingeo[gii]); if(Npt<0){Npt=0;} //because of error diffusion, it could be negative number. We restrict Npt to have min 0.
                    //error diffusion to surrounding grids with equal weights
                    diff=Ncurrent_rho_bgrid_ingeo[gii]-1.0*Npt;
                }
            }
            if(grid_ctgr<0){ //inner grid
            	generate_pt=true;
                gi=abs(grid_ctgr)-1;
                Npt=round(Ncurrent_rho_igrid[gi]); if(Npt<0){Npt=0;}  //because of error diffusion, it could be negative number. We restrict Npt to have min 0.
                diff=Ncurrent_rho_igrid[gi]-1.0*Npt;
            }
            Npt_total_temp+=Npt;
            grid_Npt[ij]=Npt;

            //error diffusion
            error_diffusion(generate_pt,diff,ij,bgrid_ingeo_outgeo_ctgr,Ncurrent_rho_igrid,Ncurrent_rho_bgrid_ingeo);
        }
    }

    //Store initial points generated via error diffusion in temp arrays. 
    //If pt_counter!=Ncurrent, delete/add points of temp arrays. 
    //Then put the final set of initial points into the xy_id and pt_ctgr arrays.
    int array_size_temp=max(Npt_total_temp,Ncurrent);
    double *xy_id_temp=new double[array_size_temp*2];
    int *pt_ctgr_temp=new int[array_size_temp];

	int inner_pt_ct_temp=inner_pt_ct;

double t0=omp_get_wtime();
    //parallel process: For each grid, generate Npt initial pts
    #pragma omp parallel num_threads(num_t) 
    {
        unsigned int seed=omp_get_thread_num()+50;
        #pragma omp for schedule(guided) reduction(+:inner_pt_ct_temp)
        for(int ij=0;ij<gnxy;ij++){
            int i=ij%gnx;
            int j=ij/gnx;
            int grid_ctgr=geo_grid(ij);
            bool generate_pt=false;
            int gi=-1;

            if(grid_ctgr>0 && grid_ctgr<=gnxy){ //bdry grid
                gi=grid_ctgr-1;
                if(bgrid_ingeo_outgeo_ctgr[gi]>0){ //a ingeo_bdry_grid
                    generate_pt=true;
                }
            }
            if(grid_ctgr<0){ //inner grid
                gi=abs(grid_ctgr)-1;
                generate_pt=true;
            }

            generate_pt_in_grid(generate_pt,grid_Npt[ij],i,j,ij,gi,grid_ctgr,pt_counter,
                xy_id_temp,pt_ctgr_temp,seed, inner_pt_ct_temp);

        }
    }
t_generate_pt_in_grid=omp_get_wtime()-t0;

    delete [] grid_Npt;
    delete [] bgrid_ingeo_outgeo_ctgr;

    inner_pt_ct=inner_pt_ct_temp;

	//-------------E. Generate the rest Ncurrent-pt_counter points-------------
	//Generate the rest Ncurrent-pt_counter points in the domain's inner/ingeo_bdry grids
	//Then, we will end up with exactly Ncurrent inner/bdry points
	int num_more_pt_needed=Ncurrent-pt_counter;
	printf("pt init generated %d pt \n",pt_counter);
	
	if(num_more_pt_needed<0){ //too many points were generated
		printf("Ncurrent: %d pt_counter: %d ERROR: pt init too many! Now delete excess points\n",Ncurrent,pt_counter);
		while(num_more_pt_needed<0){
			int pi=rand() % pt_counter; //select a random point in the range of [0,pt_counter) to delete
			if(pt_ctgr_temp[pi]==-1){ //if delete pt is inner point, decrease inner_pt_ct
				inner_pt_ct--;
			}
			int replace_pi=pt_counter-1; //replace the pt slot with the last point
			xy_id_temp[2*pi]=xy_id_temp[2*replace_pi];
			xy_id_temp[2*pi+1]=xy_id_temp[2*replace_pi+1];
			pt_ctgr_temp[pi]=pt_ctgr_temp[replace_pi];
			pt_counter--;
			num_more_pt_needed=Ncurrent-pt_counter;
		}
	}

	if(num_more_pt_needed>0){ //not enough points were generated
		printf("Ncurrent: %d pt_counter: %d ERROR: pt init too few! Now add more points\n",Ncurrent,pt_counter);
		
		while(num_more_pt_needed>0){

			//pick a random inner/ingeo_bdry grid
 			int ij=-1;
 			int gii=rand() % (geo_igrid_ct+bgrid_ingeo_ct);
 			int gi=gii;
 			if(gii<geo_igrid_ct){
 				ij=geo_igrid_ij(gi);
 			}
 			else{
 				gi=bgrid_ingeo_ind[gii-geo_igrid_ct];
 				ij=geo_bgrid_ij(gi);
 			}

 			int m;
 			int i=ij%gnx;
 	    	int j=ij/gnx;
 	    	double grid_ax=ax+(i)*gdx; 
 	        double grid_ay=ay+(j)*gdy;
 	    	double x=grid_ax+rnd()*gdx;
 		    double y=grid_ay+rnd()*gdy;

 		    if(gii<geo_igrid_ct){ //inner grid 
 		    	m=pt_counter++;
 		    	num_more_pt_needed--;

 		    	xy_id_temp[2*m]=x;
 	    		xy_id_temp[2*m+1]=y;

 	    		inner_pt_ct++;
 	    		pt_ctgr_temp[m]=-1; //-1: inner pt
 		    	
 		    }
 		    else{ //bdry grid
 		    	double finalx=x;double finaly=y;
 		    	bool store=false;
 		    	double pt_sdf_temp=sdf_adf(x,y);
 		    	if(pt_sdf_temp>bgrid_geps[gi]){ //pt outside
 		    		bool project_pt_outside_geo=true;
 		    		if(pt_projection(finalx, finaly,x,y,bgrid_deps[gi],bgrid_geps[gi],project_pt_outside_geo)){
 		    			int final_i=(finalx-ax)*inv_gdx;
						int final_j=(finaly-ay)*inv_gdy;
						int final_ij=gnx*final_j+final_i;
						if(final_ij==ij)
						{
							store=true;
							pt_sdf_temp=sdf_adf(finalx,finaly);
						}
 		    		}
 		    	}
 		    	else{ //pt in/on geometry
 		    		store=true;
 		    	}
 		    	if(store==true){
 		    		m=pt_counter++;
 		    		num_more_pt_needed--;

 		    		xy_id_temp[2*m]=finalx;
 	    			xy_id_temp[2*m+1]=finaly;

 	    			if(fabs(pt_sdf_temp)<=bgrid_geps[gi]){//bdry pt
 	    				pt_ctgr_temp[m]=1; //1: bdry pt
 	    			}
 	    			else if(pt_sdf_temp<-bgrid_geps[gi]){
 	    				inner_pt_ct++;
 	    				pt_ctgr_temp[m]=-1; //-1: inner pt
 	    			}
 	    			else{//should not happen
                		printf("ERROR: pt generated not inner nor bdry pt! \n");
                		throw std::exception();
                	}

 		    	}
 		    }
		}
	 }

	//exactly Ncurrent points were generated and is now being put into container; shifted rightwards of Nfixed space
	#pragma omp parallel for num_threads(num_t)
	for(int pi=0;pi<Ncurrent;pi++){
		int pii=Nfixed+pi;
		xy_id[2*pii]=xy_id_temp[2*pi];
		xy_id[2*pii+1]=xy_id_temp[2*pi+1];
		pt_ctgr[pii]=pt_ctgr_temp[pi];
	}

	//Add Nfixed fixed points into the intial points, updating xy_id and pt_ctgr. 
	//The first Nfixed spaces in xy_id and pt_ctgr are the fixed points.
	//This also updates Ncurrent and Nremain,
	//And updates chrtrt_len_h[ij], bgrid_geps[gi], bgrid_deps[gi].
	add_fixed_points_procedure();

	con->add_parallel(xy_id, Ncurrent, num_t);
    con->put_reconcile_overflow();

    delete [] xy_id_temp;
    delete [] pt_ctgr_temp;
	delete [] density_igrid;
	delete [] density_bgrid_ingeo;
	delete [] Ncurrent_rho_igrid;
	delete [] Ncurrent_rho_bgrid_ingeo;

/*
//check and print out only adf_bdry with a finer grid
    char fadf[256];
    sprintf(fadf,"%s/adf_field.txt",file_name_prefix);
    FILE *fadfout=fopen(fadf,"a");

    for(int i=0;i<10*gnx;i++){
        for(int j=0;j<10*gny;j++){
        	double grid_mx=ax+(i+0.5)*0.1*gdx; 
        	double grid_my=ay+(j+0.5)*0.1*gdy; 
        	int ind=shp->pt_geo_grid_val(grid_mx,grid_my)-1;
        	if(ind>=0 && ind<gnxy){ //bdry grid
        		fprintf(fadfout,"%g ",geo_bgrid_adf(grid_mx,grid_my));
        	}
        	else{
        		fprintf(fadfout,"%g ",0.0);
        	}

        }
        fprintf(fadfout,"\n");
        
    }
    fclose(fadfout);
*/

}

/**
 * @brief Calculate the value from the adaptive distance field (ADF) grid at a given point within the geometry boundary grid.
 * 
 * It determines the grid index ind based on the point's location within the geometry boundary grid using shp->pt_geo_grid_val(x, y).
 * If the grid index ind is less than zero or greater than or equal to gnxy (number of grid cells), it throws an exception.
 * It retrieves and returns the value from the ADF grid at the calculated index using bgrid_adf[ind]->getVal(x, y).
 * 
 * @param x The x-coordinate of the point.
 * @param y The y-coordinate of the point.
 * @return The value from the adaptive distance field (ADF) grid at the given point within the geometry boundary grid.
 * @throws std::exception if the point is not within the geometry boundary grid.
 */
double parallel_meshing_2d::geo_bgrid_adf(double x, double y){
	int ind=shp->pt_geo_grid_val(x,y)-1;
	if(ind<0 || ind>=gnxy){
		printf("ERROR 3: pt (%g, %g) not in geometry bdry grid, grid ind %d.\n",x,y,ind);
    	throw std::exception();
	}
	return bgrid_adf[ind]->getVal(x,y);
}


/**
 * @brief Calculate the signed distance field (SDF) value at a given point.
 * 
 * It determines the grid index ind based on the point's location within the geometry boundary grid.
 * If the point is within the boundary grid, it returns the value obtained from the adaptive distance field (ADF) grid at the calculated index.
 * If the point is outside the boundary grid, it invokes the sdf function to compute the SDF value using the point's coordinates (x, y).
 * 
 * @param x The x-coordinate of the point.
 * @param y The y-coordinate of the point.
 * @return The signed distance field (SDF) value at the given point.
 */
double parallel_meshing_2d::sdf_adf(double x, double y){
	int ind=shp->pt_geo_grid_val(x,y)-1;
	if(ind>=0 && ind<gnxy){ //bdry grid
		return bgrid_adf[ind]->getVal(x,y);
	}
	else{
		return sdf(x,y);
	}
}


/**
 * @brief Get the characteristic edge length value at a given point.
 * 
 * It calculates the indices i and j based on the x and y coordinates.
 * It computes the grid index ij based on i and j.
 * It returns the characteristic edge length value at the computed index ij.
 * 
 * @param x The x-coordinate of the point.
 * @param y The y-coordinate of the point.
 * @return The characteristic edge length value at the given point.
 */
double parallel_meshing_2d::get_chrtrt_len_h(double x, double y){
	int i=(x-ax)*inv_gdx;
	int j=(y-ay)*inv_gdy;
	int ij=gnx*j+i;
	return chrtrt_len_h[ij];
}


/**
 * @brief Prints the initial points to a file.
 *
 * It appends the coordinates of the initial points to a file named "initial_pts.txt" in the specified directory.
 *
 * @param case_name The name of the directory where the file will be created.
 */
void parallel_meshing_2d::print_init_pts_to_file(const char *case_name){
	char fp[256];
    sprintf(fp,"%s/initial_pts.txt",case_name);
    FILE *fpout=fopen(fp,"a");
    for(int i=0;i<Ncurrent;i++){
        double x=xy_id[2*i];
        double y=xy_id[2*i+1];
        fprintf(fpout,"%g %g \n",x,y);
    }
    fclose(fpout);
}




