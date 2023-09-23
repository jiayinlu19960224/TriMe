#include "sizing_2d.hh"

using namespace voro;

/**
 * @brief Constructor for the `sizing_2d` class.
 *
 * @param shp_ A pointer to the `shape_2d` object providing the underlying geometry grid information.
 */
sizing_2d::sizing_2d(shape_2d *shp_)
	: num_t(shp_->num_t), shp(shp_), 

	ax(shp->ax),bx(shp->bx),ay(shp->ay),by(shp->by),
	gnx(shp->gnx),gny(shp->gny),gnxy(shp->gnxy),
	gdx(shp->gdx),gdy(shp->gdy),diag_gdxy(shp->diag_gdxy),
    inv_gdx(shp->inv_gdx),inv_gdy(shp->inv_gdy),
	geo_bgrid_ct(shp->geo_bgrid_ct),geo_igrid_ct(shp->geo_igrid_ct),
    geo_ogrid_ct(shp->geo_ogrid_ct)

	{}

/**
 * @brief Destructor for the `sizing_2d_automatic` class.
 *
 * Cleans up dynamically allocated memory used by the class.
 */
sizing_2d_automatic::~sizing_2d_automatic(){
	delete [] igrid_sizing;
	delete [] igrid_density;
	delete [] bgrid_sizing;
	delete [] bgrid_density;
    delete [] ogrid_sizing;
    delete [] ogrid_density;
	delete [] bdry_pt_xy;
	delete [] bdry_pt_lfs;
	delete [] medial_pt_xy;
}

/**
 * @brief Checks if a point is inside the outer grid.
 *
 * Determines whether a given point (x, y) is located inside the outer grid.
 *
 * @param x The x-coordinate of the point.
 * @param y The y-coordinate of the point.
 * @return True if the point is inside the outer grid, false otherwise.
 */
bool sizing_2d::pt_in_outer_grid(double x,double y){
    int geo_grid_val=shp->pt_geo_grid_val(x,y);
    if(geo_grid_val>gnxy){ //outer grid
        return true;
    }
    else{
        return false;
    }
}


/**
 * @brief Retrieves the boundary points for sizing in 2D automatic mode.
 *
 * This function collects the boundary points required for automatic sizing field calculation for a 2D shape.
 * It calculates the characteristic lengthscale for the boundary edge segments and performs
 * necessary calculations to determine the boundary points that are evenly spaced based on the lengthscales along the boundary.
 *
 */
void sizing_2d_automatic::get_bdry_pts(){

    //characteristic lengthscale for the bdry edge segments: 
    double hbdry=min(gdx,gdy); 
    //imagine bdry grids on the inner and outer sides and the one through the bdry itself, the points all collapes 
    //on the same line segment approximated by diag_gdxy 
    deps_temp=sqrt(eps)*hbdry; //deltaX,deltaY: used in projection step
    geps_temp=geps_h_frac*hbdry; //geometry boundary tolerance

    //collect bdry pt cts and xy's from parallel threads
    int *bdry_ct_private=new int[num_t]; for(int i=0;i<num_t;i++){bdry_ct_private[i]=0;}
    std::vector<double> *bdry_xy_private=new vector<double>[num_t];

    //if shape is custom shape from contour lines, 
    //use contour lines bdry pts directly as bdry pts.
    if(shp->is_custom_shape_contour==true){

        //loop through all boundaries
        for(int i=0;i<shp->b_ct();i++){

            #pragma omp parallel num_threads(num_t)
            {
                #ifdef _OPENMP
                    const int ithread = omp_get_thread_num();
                #else
                    const int ithread=0;
                #endif

                //loop through all line segments of the current boundary
                #pragma omp for schedule(dynamic)
                for(int j=0;j<shp->seg_ct(i);j++){
                    //get pt at begining of line segment
                    double x0=shp->b_pts_x(i,j);
                    double y0=shp->b_pts_y(i,j);
                    bdry_xy_private[ithread].push_back(x0);
                    bdry_xy_private[ithread].push_back(y0);
                    bdry_ct_private[ithread]++;
                    //get pt at end of line segment
                    double x1=shp->b_pts_x(i,j+1);
                    double y1=shp->b_pts_y(i,j+1);
                    bdry_xy_private[ithread].push_back(x1);
                    bdry_xy_private[ithread].push_back(y1);
                    bdry_ct_private[ithread]++;
                    //calculate length of current line segment
                    double l_len=d_points(x0,y0,x1,y1);
                    //if the two points are too far away
                    if(l_len>hbdry){
                        //generate N_more points evenly spaced inside line segment
                        //collect them as bdry pts
                        //so each sub-line-segment is less than hbdry away
                        int N_more=int(l_len/hbdry);
                        double abdx=(x1-x0)/(N_more+1);
                        double abdy=(y1-y0)/(N_more+1);
                        for (int N_more_i=1;N_more_i<=N_more;N_more_i++){
                            double xn=x0+abdx*N_more_i;
                            double yn=y0+abdy*N_more_i;
                            bdry_xy_private[ithread].push_back(xn);
                            bdry_xy_private[ithread].push_back(yn);
                            bdry_ct_private[ithread]++;
                        }
                    }
                }
            }

        }
    }
    else{  //SDF is function defined, no bdry points readily available.
        bool project_pt_outside_geo=false; //since bdry grid may generate points inside geo/on bdry as well
    	#pragma omp parallel num_threads(num_t)
    	{
    		#ifdef _OPENMP
                const int ithread = omp_get_thread_num();
            #else
                const int ithread=0;
            #endif
    		unsigned int seed=ithread+5;

    		//loop through bdry grids
    		#pragma omp for schedule(dynamic)
    		for(int gridi=0;gridi<geo_bgrid_ct;gridi++){
    			int ij=geo_bgrid_ij(gridi);
    			int i=ij%gnx;
    	    	int j=ij/gnx;

                double grid_ax=ax+(i)*gdx; 
                double grid_bx=ax+(1+i)*gdx; 
                double grid_ay=ay+(j)*gdy;
                double grid_by=ay+(1+j)*gdy;

                //generate in the bdry grid, Npts_grid_init random points
                for(int pi=0;pi<Npts_grid_init;){

                    //random point in the bdry grid
                    double x=grid_ax+rnd_r(&seed)*gdx;
                    double y=grid_ay+rnd_r(&seed)*gdy;
                    double new_x,new_y;

                    //project the point onto the geometry bdry
                    //if projection successful, collect the point as a bdry pt
                    if(pt_projection(new_x,new_y,x,y,-1,-1,project_pt_outside_geo)){
                        bdry_ct_private[ithread]++;
                        bdry_xy_private[ithread].push_back(new_x);
                        bdry_xy_private[ithread].push_back(new_y);
                        pi++;
                    }
                }
    	    }
    	}
    }

        
    //to store the initial set of bdry points before trimming
    int bdry_pt_ct_temp=0;
    double *bdry_pt_xy_temp;
    //store the pt ind in a temporary geometry grid, 
    //which will help with finding closest pts efficiently later
    std::vector<int> *grid_pt_ind=new vector<int>[gnxy];

	//gather all bdry points results from all threads, 
    //store into variables bdry_pt_xy_temp and bdry_pt_ct_temp:
    for(int i=0;i<num_t;i++){bdry_pt_ct_temp+=bdry_ct_private[i];}
    bdry_pt_xy_temp=new double[bdry_pt_ct_temp*2];
    double *p=bdry_pt_xy_temp;
    int pt_ind=0;
    for(int i=0;i<num_t;i++){
        for(int j=0;j<bdry_ct_private[i];j++){
            double px=bdry_xy_private[i][2*j];
            double py=bdry_xy_private[i][2*j+1];
            int gridi=(px-ax)*inv_gdx;
            int gridj=(py-ay)*inv_gdy;
            grid_pt_ind[gridj*gnx+gridi].push_back(pt_ind);
            *p=px;
            p++;
            *p=py;
            p++;
            pt_ind++;
        }
    }
    delete [] bdry_ct_private;
    delete [] bdry_xy_private;


    char fb[256];
    sprintf(fb,"bdry_pts_old.txt");
    FILE *fbout=fopen(fb,"a");
    
    for(int j=0;j<bdry_pt_ct_temp;j++){
        double bx=bdry_pt_xy_temp[2*j];
        double by=bdry_pt_xy_temp[2*j+1];
        fprintf(fbout,"%g %g \n",bx,by);
    }
    fclose(fbout);

    //===========================Process points: So they are not too close together =============
    //Now, we have initial set of bdry points, stored in:
        //bdry_pt_xy_temp
        //bdry_pt_ct_temp
    //They should approximate bdry fine features well and are not too far apart, since for 
        //shp contour input: It uses contour line pts, and add in pts in between if too far away
        //sdf function input: the initial pt generation procedure would generate bdry pts "crowding" along the bdry
    //Now, we will loop through and trim pts, so pts on the same bdry are at least bdry_pt_dis_fac*hbdry apart.
    //After this, the bdry pts will be evenly/nicely spaced. 
    
    double max_len_dummy=(bx-ax)+(by-ay);
    int *pt_tag=new int[bdry_pt_ct_temp];
    #pragma omp parallel for num_threads(num_t)
    for(int pi=0;pi<bdry_pt_ct_temp;pi++){pt_tag[pi]=0;}
    bool Continue=true;
    int iter_ct=0;
    while(Continue){
        Continue=false;
        iter_ct++;

        
        for(int pi=0;pi<bdry_pt_ct_temp;pi++){
            bool check=false;
            if(iter_ct==1){
                if(pt_tag[pi]!=-1){check=true;}
            }
            else{
                if(pt_tag[pi]==1){check=true;}
            }
            if(check){
                bool continue_find=true;
                //current pt info
                double px=bdry_pt_xy_temp[2*pi];
                double py=bdry_pt_xy_temp[2*pi+1];
                int pgridi=(px-ax)*inv_gdx;
                int pgridj=(py-ay)*inv_gdy;
                //closest nei pt and dis
                double closest_dis=max_len_dummy;
                int closest_pt_ind=-1;
                //outward grid search procedure
                int incre=0;
                bool continue_search_local=true;
                bool found_nei=false;
                //current pt id: grid_pt_ind[current_pt_grid_ind][current_pt_grid_ind_location]
                int current_pt_grid_ind=-1; 
                int current_pt_grid_ind_location=-1;
                //closest nei id: grid_pt_ind[closest_pt_grid_ind][closest_pt_grid_ind_location]
                int closest_pt_grid_ind=-1;
                int closest_pt_grid_ind_location=-1;
                while(continue_search_local){
                    //if have found a neighboring pt,
                    //just need to check one more layer grids outward, 
                    //then can stop
                    if(found_nei==true){continue_search_local=false;}
                    int searchli=pgridi-incre;
                    int searchhi=pgridi+incre;
                    int searchlj=pgridj-incre;
                    int searchhj=pgridj+incre;

                    for(int searchi=searchli;searchi<=searchhi;searchi++){
                        for(int searchj=searchlj;searchj<=searchhj;searchj++){
                            //just search the new layer of grids
                            if(searchi==searchli || searchi==searchhi || searchj==searchlj || searchj==searchhj){

                                //grid not out of bound
                                if(searchi>=0 && searchi <gnx && searchj>=0 && searchj<gny){
                                    //if there are pts in the grid
                                    int searchij=searchj*gnx+searchi;

                                    int grid_num_pt=(signed int)grid_pt_ind[searchij].size();
                                    for(int neii=0;neii<grid_num_pt;neii++){
                                        int nei_id=grid_pt_ind[searchij][neii];
                                        if(nei_id!=-1){ //invalid/deleted point
                                            if(nei_id!=pi){ //exclude pt itself
                                                double neix=bdry_pt_xy_temp[2*nei_id];
                                                double neiy=bdry_pt_xy_temp[2*nei_id+1];
                                                double nei_dis=d_points(px,py,neix,neiy);
                                                if(nei_dis<closest_dis){
                                                    closest_dis=nei_dis;
                                                    closest_pt_ind=nei_id;
                                                    found_nei=true;
                                                    closest_pt_grid_ind=searchij;
                                                    closest_pt_grid_ind_location=neii;
                                                }
                                            }
                                            else{
                                                current_pt_grid_ind=searchij;
                                                current_pt_grid_ind_location=neii;

                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    //have exhausted grid, stop search
                    if(searchli<0 && searchhi>=gnx && searchlj<0 && searchhj>=gny){
                        continue_search_local=false;
                    }
                    incre++;
                }

                while(continue_find){

                    if(closest_dis<bdry_pt_dis_fac*hbdry){

                        //find the midpoint between the bdry pt and its current closest pt
                        double mid_x=0.5*(px+bdry_pt_xy_temp[2*closest_pt_ind]);
                        double mid_y=0.5*(py+bdry_pt_xy_temp[2*closest_pt_ind+1]);
                        if(abs(sdf(mid_x,mid_y))<geps_temp){

                            //move pt to midpoint
                            bdry_pt_xy_temp[2*pi]=mid_x;
                            bdry_pt_xy_temp[2*pi+1]=mid_y;

                            //tag pt [1]
                            pt_tag[pi]=1;

                            //update pt_grid_ind array: 
                            //delete the deleted pt in the grid: indicate with "-1"
                            grid_pt_ind[closest_pt_grid_ind][closest_pt_grid_ind_location]=-1;


                            //and delete the current pt old position in the grid: indicate with "-1"
                            //update current pt position in grid
                            int new_pgridi=(mid_x-ax)*inv_gdx;
                            int new_pgridj=(mid_y-ay)*inv_gdy;

                            if(new_pgridi!=pgridi || new_pgridj!=pgridj){

                                grid_pt_ind[current_pt_grid_ind][current_pt_grid_ind_location]=-1;
                                grid_pt_ind[new_pgridj*gnx+new_pgridi].push_back(pi);
                            }

                            //delete nei pt
                            //tag pt [-1]
                            pt_tag[closest_pt_ind]=-1;

                            continue_find=false;
                            Continue=true;
                        }

                    }
                    else{
                        //done with the point: No need to further worry about being too close to other bdry pts
                        continue_find=false;
                        pt_tag[pi]=0;
                    }
                    if(continue_find==true){
                        //the next closest pt should be farther away than old_closest_dis
                        double old_closest_dis=closest_dis;
                        int old_closest_ind=closest_pt_ind;
                        px=bdry_pt_xy_temp[2*pi];
                        py=bdry_pt_xy_temp[2*pi+1];
                        pgridi=(px-ax)*inv_gdx;
                        pgridj=(py-ay)*inv_gdy;
                        //reset closest pt variables
                        closest_dis=max_len_dummy;
                        closest_pt_ind=-1;
                        //find next closest pt & dis
                        //outward grid search procedure
                        incre=0;
                        continue_search_local=true;
                        found_nei=false;
                        //closest nei id: grid_pt_ind[closest_pt_grid_ind][closest_pt_grid_ind_location]
                        closest_pt_grid_ind=-1;
                        closest_pt_grid_ind_location=-1;
                        while(continue_search_local){
                            //if have found a new close neighboring pt,
                            //just need to check one more layer grids outward, 
                            //then can stop
                            if(found_nei==true){continue_search_local=false;}
                            int searchli=pgridi-incre;
                            int searchhi=pgridi+incre;
                            int searchlj=pgridj-incre;
                            int searchhj=pgridj+incre;
                            for(int searchi=searchli;searchi<=searchhi;searchi++){
                                for(int searchj=searchlj;searchj<=searchhi;searchj++){
                                    //just search the new layer of grids
                                    if(searchi==searchli || searchi==searchhi || searchj==searchlj || searchj==searchhj){
                                        //grid not out of bound
                                        if(searchi>=0 && searchi <gnx && searchj>=0 && searchj<gny){
                                            //if there are pts in the grid
                                            int searchij=searchj*gnx+searchi;
                                            int grid_num_pt=(signed int)grid_pt_ind[searchij].size();
                                            for(int neii=0;neii<grid_num_pt;neii++){
                                                int nei_id=grid_pt_ind[searchij][neii];
                                                if(nei_id!=-1){
                                                    if(nei_id!=pi){ //exclude pt itself
                                                        double neix=bdry_pt_xy_temp[2*nei_id];
                                                        double neiy=bdry_pt_xy_temp[2*nei_id+1];
                                                        double nei_dis=d_points(px,py,neix,neiy);
                                                        if(nei_dis<closest_dis && nei_dis>old_closest_dis+(1e-12) && nei_id!=old_closest_ind){
                                                            closest_dis=nei_dis;
                                                            closest_pt_ind=nei_id;
                                                            found_nei=true;
                                                            closest_pt_grid_ind=searchij;
                                                            closest_pt_grid_ind_location=neii;
                                                        }
                                                    }
                                                    else{
                                                        current_pt_grid_ind=searchij;
                                                        current_pt_grid_ind_location=neii;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            //have exhausted grid, stop search
                            if(searchli<0 && searchhi>=gnx && searchlj<0 && searchhj>=gny){
                                continue_search_local=false;
                            }
                            incre++;
                        }

                    }
                    
                }
            }
        }
    }


    //now, get the final set of bdry pts, by keeping pts with tag [0]
    bdry_pt_ct=0;
    for(int j=0;j<bdry_pt_ct_temp;j++){
        if(pt_tag[j]==0){
            bdry_pt_ct++;
        }
    }
    bdry_pt_xy=new double[2*bdry_pt_ct];
    int ct_temp=0;
    #pragma omp parallel for num_threads(num_t)
    for(int j=0;j<bdry_pt_ct_temp;j++){
        if(pt_tag[j]==0){
            double x=bdry_pt_xy_temp[2*j];
            double y=bdry_pt_xy_temp[2*j+1];
            int m;
            #pragma omp atomic capture
            m=ct_temp++;
            bdry_pt_xy[2*m]=x;
            bdry_pt_xy[2*m+1]=y;
        }
    }

    delete [] pt_tag;
    delete [] bdry_pt_xy_temp;
    delete [] grid_pt_ind;

}



/**
 * @brief Computes the medial axis points of a shape using a Voronoi diagram of evenly spaced shape boundary points.
 */
void sizing_2d_automatic::get_medial_pts(){

    //characteristic lengthscale for the bdry edge segments: 
    double hbdry=min(gdx,gdy); 
    geps_temp=geps_h_frac*hbdry; //geometry boundary tolerance

	//A. Create Voronoi diagram using the boundary edge points of the shape
    //create container, using Nopt for localized particle case
    int cnx=sqrt(bdry_pt_ct/30.0); int cny=cnx;
    container_2d con(ax,bx,ay,by,cnx,cny,false,false,16,num_t);

    con.add_parallel(bdry_pt_xy, bdry_pt_ct, num_t);
    con.put_reconcile_overflow();
    
    container_2d::iterator cli;
    
    int *medial_ct_private=new int[num_t]; 
    std::vector<double> *medial_xy_private=new vector<double>[num_t];
    double min_len_dummy=-((bx-ax)+(by-ay));
    //bdry_tol: medial pts should be at least bdry_tol away from geo bdry
    double bdry_tol=bdry_pt_dis_fac*hbdry;   
    double medial_trim_region_check_r=medial_trim_region_fac*medial_cutoff_nei_ct*bdry_tol;
    int medial_region_check_grid_ct=ceil(medial_trim_region_check_r/hbdry);
printf("medial_trim_region_check_r: %g, medial_region_check_grid_ct: %d \n", medial_trim_region_check_r,medial_region_check_grid_ct);
    
    #pragma omp parallel num_threads(num_t)
    {
        #ifdef _OPENMP
            const int ithread = omp_get_thread_num();
        #else
            const int ithread=0;
        #endif
        //create thread-private variables here
        voronoicell_neighbor_2d c(con);
        medial_ct_private[ithread]=0;

        //Loop through every particle in the container, by
        //looping through iterators of particle in container
        #pragma omp for schedule(dynamic)
        for(cli=con.begin(); cli<con.end(); cli++){
            //if voronoi cell true:
            if(con.compute_cell(c, cli)) {
                //getting (x, y) position and id of the particle (ij, q)
                int ijk=(*cli).ijk;
                int q=(*cli).q;    
                int id;double x,y;
                double *pp=con.p[ijk]+con.ps*q;
                x=*(pp++);y=*pp;
                id=con.id[ijk][q];

                //get all the vertices of the Voronoi cell
                std::vector<double> voro_vertices_x;
                std::vector<double> voro_vertices_y;
                if(c.p!=0){
                    int k=0;
                    do {
                        voro_vertices_x.push_back(x+0.5*c.pts[2*k]);   //i<<2 is the same as i*4
                        voro_vertices_y.push_back(y+0.5*c.pts[2*k+1]);
                        k=c.ed[2*k];
                    } while (k!=0);
                    voro_vertices_x.push_back(x+0.5*c.pts[0]);  //back to the first vertice, making a closed cell
                    voro_vertices_y.push_back(y+0.5*c.pts[1]);
                }

                int voro_vertice_ct=(signed int)voro_vertices_x.size()-1;
                
                double dis_max=min_len_dummy;
                double vx_max, vy_max;
                double *vdis=new double[voro_vertice_ct];
                //find the vertex with greatest distance to particle
                for(int vi=0;vi<voro_vertice_ct;vi++){
                    double vx=voro_vertices_x[vi];
                    double vy=voro_vertices_y[vi];
                    vdis[vi]=d_points(x,y,vx,vy);
                    if(vdis[vi]>dis_max){
                        dis_max=vdis[vi];
                        vx_max=vx; 
                        vy_max=vy;
                    }
                }

                //if the vertex found is not on domain boundary
                //nor too close to the shape bdry, 
                //store it as a medial axis point
                if(vx_max>ax+geps_temp && vy_max>ay+geps_temp 
                    && vx_max<bx-geps_temp && vy_max<by-geps_temp){
                    if(fabs(sdf(vx_max,vy_max))>bdry_tol){
                    	//thread private copy of counting and storing space
                        medial_ct_private[ithread]++;
                        medial_xy_private[ithread].push_back(vx_max);
                        medial_xy_private[ithread].push_back(vy_max);
                    }
                }

                //find the vertex with greatest distance to particle on the other side of
                //the line through the particle orthogonal to line defined by (x,y) and (vx_max, vy_max) 
                //We can use obtuse triangle: the other vertex must form 
                //an obtuse triangle with the current point (x,y) and the previously found (vx_max, vy_max)
                //so should satisfy c^2>a^2+b^2
                dis_max=min_len_dummy;
                double dis_a=d_points(x,y,vx_max,vy_max);
                double vx_max_1=vx_max; double vy_max_1=vy_max;
                for(int vi=0;vi<voro_vertice_ct;vi++){
                    double vx=voro_vertices_x[vi];
                    double vy=voro_vertices_y[vi];
                    double dis_b=d_points(x,y,vx,vy);
                    double dis_c=d_points(vx_max_1,vy_max_1,vx,vy);
                    if(sqr(dis_c)>sqr(dis_a)+sqr(dis_b)){ //obtuse triangle
                        if(vdis[vi]>dis_max){
                            dis_max=vdis[vi];
                            vx_max=vx; 
                            vy_max=vy;
                        }
                    }
                }

                //if the vertex found is not on domain boundary
                //nor too close to the shape bdry, 
                //store it as a medial axis point
                if(vx_max>ax+geps_temp && vy_max>ay+geps_temp 
                    && vx_max<bx-geps_temp && vy_max<by-geps_temp){
                    if(fabs(sdf(vx_max,vy_max))>bdry_tol){
                    	//thread private copy of counting and storing space
                        medial_ct_private[ithread]++;
                        medial_xy_private[ithread].push_back(vx_max);
                        medial_xy_private[ithread].push_back(vy_max);
                    }
                }

                delete [] vdis;
            }
        }
    }
    
    //gather all medial axis points results from all threads, into variables medial_pt_xy and medial_pt_ct:
    //and sort medial points into a grid structure, grid_medial_xy_temp
    int medial_pt_ct_temp=0;
    for(int i=0;i<num_t;i++){medial_pt_ct_temp+=medial_ct_private[i];}
    double *medial_pt_xy_temp=new double[medial_pt_ct_temp*2];
    
    //store medial pt xy's in the medial_pt_xy_temp array
    int ct_temp=0;
    for(int i=0;i<num_t;i++){
        #pragma omp parallel for num_threads(num_t)
        for(int j=0;j<medial_ct_private[i];j++){
            double px=medial_xy_private[i][2*j];
            double py=medial_xy_private[i][2*j+1];
            int m;
            #pragma omp atomic capture
            m=ct_temp++;
            medial_pt_xy_temp[2*m]=px;
            medial_pt_xy_temp[2*m+1]=py;
        }
    }
    delete [] medial_ct_private;
    delete [] medial_xy_private;

    //Used for trimming medial pts later: 
    //grid_medial_xy_temp[0]: stores x0,y0,x1,y1, etc, medial pts xy's in this geo grid.
    //to test if there are other neighbor medial pts of the current medial pt.
    //criteria: within medial_trim_region_check_r there should be at least medial_cutoff_nei_ct neighbors
    std::vector<double> *grid_medial_xy_temp=new std::vector<double>[gnxy];
    std::vector<int> *grid_medial_id_temp=new std::vector<int>[gnxy];

    for(int j=0;j<medial_pt_ct_temp;j++){
        double px=medial_pt_xy_temp[2*j];
        double py=medial_pt_xy_temp[2*j+1];

        //find out grid the pt is in; increase the pt ct in the grid.
        int gridi=(px-ax)*inv_gdx;
        int gridj=(py-ay)*inv_gdy;
        int gridij=gnx*gridj+gridi;
        grid_medial_id_temp[gridij].push_back(j);
        grid_medial_xy_temp[gridij].push_back(px);
        grid_medial_xy_temp[gridij].push_back(py);
    }


char bug0[256];
sprintf(bug0,"bdry_pts_voro_diagram.gnu");
con.draw_cells_gnuplot(bug0);


    con.clear();


    
    char fm[256];
    sprintf(fm,"medial_pts_old.txt");
    FILE *fmout=fopen(fm,"a");
    
    for(int j=0;j<medial_pt_ct_temp;j++){
        double mx=medial_pt_xy_temp[2*j];
        double my=medial_pt_xy_temp[2*j+1];
        fprintf(fmout,"%g %g \n",mx,my);
    }
    fclose(fmout);
    


    //tag: 0 keep; -1 discard
    int *pt_tag=new int[medial_pt_ct_temp];
    #pragma omp parallel for num_threads(num_t)
    for(int i=0;i<medial_pt_ct_temp;i++){pt_tag[i]=0;}

    //trim medial axis pts:
    //if pt's surrounding region has less than medial_cutoff_nei_ct medial axis pts, delete the pt
    //because we assume medial axis pts should be close to some other medial axis points, to form the skeleton
    bool continue_medial_trim=true;
    while(continue_medial_trim){
        continue_medial_trim=false;
        int discard_pt_ct=0;

        #pragma omp parallel for num_threads(num_t) schedule(dynamic)
        for(int j=0;j<medial_pt_ct_temp;j++){
            if(pt_tag[j]!=-1){
                double mx=medial_pt_xy_temp[2*j];
                double my=medial_pt_xy_temp[2*j+1];
                int gridi=(mx-ax)*inv_gdx;
                int gridj=(my-ay)*inv_gdy;
                int gridij=gnx*gridj+gridi;

                //find the range of geo grids to check
                //the circular region to check for neighbors has medial_trim_region_check_r=medial_cutoff_nei_ct*bdry_tol
                int gridli=gridi-medial_region_check_grid_ct; int gridhi=gridi+medial_region_check_grid_ct;
                int gridlj=gridj-medial_region_check_grid_ct; int gridhj=gridj+medial_region_check_grid_ct;
                
                int nei_ct=0;
                //loop through the local region including the surrounding grid
                //count number of neiboring pts 
                for(int gridii=gridli;gridii<=gridhi;gridii++){
                    for(int gridjj=gridlj;gridjj<=gridhj;gridjj++){
                        if(gridii>=0 && gridii<gnx && gridjj>=0 && gridjj<gny){
                            int gridiijj=gnx*gridjj+gridii;

                            for(int grid_pt_i=0;grid_pt_i<(signed int)grid_medial_id_temp[gridiijj].size();grid_pt_i++){
                                if(grid_medial_id_temp[gridiijj][grid_pt_i]!=-1){
                                    //medial point trmming criteria
                                    //criteria: within medial_trim_region_check_r there should be at least medial_cutoff_nei_ct neighbors
                                    double grid_pt_x=grid_medial_xy_temp[gridiijj][2*grid_pt_i];
                                    double grid_pt_y=grid_medial_xy_temp[gridiijj][2*grid_pt_i+1];
                                    double grid_pt_dis=d_points(mx,my,grid_pt_x,grid_pt_y);
                                    //check in circular region and not overlapping
                                    if(grid_pt_dis<medial_trim_region_check_r
                                        && grid_pt_dis>1e-14){
                                        nei_ct++;
                                    }
                                }
                                
                            }
                        }
                    }
                }

                //if number of neighbors less than cutoff value, discard the pt
                if(nei_ct<medial_cutoff_nei_ct){
                    pt_tag[j]=-1;
                    #pragma omp atomic
                    discard_pt_ct++;
                    //update temp medial pts arrays, and remove the current pt
                    bool erase_continue=true;
                    int grid_pt_i=0;

                    while(erase_continue==true && grid_pt_i<(signed int)grid_medial_id_temp[gridij].size()){
                        if(grid_medial_id_temp[gridij][grid_pt_i]==j){
                            grid_medial_id_temp[gridij][grid_pt_i]=-1;
                            erase_continue=false;
                        }
                        grid_pt_i++;
                    }
                }
            }
        }
        if(discard_pt_ct>0){
            continue_medial_trim=true;
        }
    }
    

    //gather the final set of medial axis pts, after trimming
    medial_pt_ct=0;
    for(int j=0;j<medial_pt_ct_temp;j++){
        if(pt_tag[j]!=-1){
            medial_pt_ct++;
        }
    }
    medial_pt_xy=new double[2*medial_pt_ct];
    ct_temp=0;
    #pragma omp parallel for num_threads(num_t)
    for(int j=0;j<medial_pt_ct_temp;j++){
        if(pt_tag[j]!=-1){
            double mx=medial_pt_xy_temp[2*j];
            double my=medial_pt_xy_temp[2*j+1];
            int m;
            #pragma omp atomic capture
            m=ct_temp++;
            medial_pt_xy[2*m]=mx;
            medial_pt_xy[2*m+1]=my;
        }
    }

    delete [] pt_tag;
    delete [] medial_pt_xy_temp;
    delete [] grid_medial_xy_temp;
    delete [] grid_medial_id_temp;
}

/**
 * @brief Calculates the lengths from boundary points to the closest medial axis points: local feature size values.
 * 
 * This function calculates the lengths from each boundary point to the closest
 * medial axis point.
 */
void sizing_2d_automatic::get_bdry_pt_lfs(){
	bdry_pt_lfs=new double[bdry_pt_ct];
	double max_len_dummy=(bx-ax)+(by-ay);
    //loop through bdry pts
    #pragma omp parallel for num_threads(num_t)
    for(int i=0;i<bdry_pt_ct;i++){
        double x=bdry_pt_xy[2*i];
        double y=bdry_pt_xy[2*i+1];
        double min_dis=max_len_dummy;
        for(int j=0;j<medial_pt_ct;j++){
            double mx=medial_pt_xy[2*j];
            double my=medial_pt_xy[2*j+1];
            double dis=d_points(x,y,mx,my);
            if(dis<min_dis){min_dis=dis;}
        }
        bdry_pt_lfs[i]=min_dis;
    }
}

/**
 * @brief Calculates the sizing grid for both boundary and inner grids.
 */
void sizing_2d_automatic::get_sizing_grid(){
	double max_len_dummy=(bx-ax)+(by-ay);

	//sizing grid for bdry grids
	#pragma omp parallel for num_threads(num_t)
	for(int gridi=0;gridi<geo_bgrid_ct;gridi++){
		int ij=geo_bgrid_ij(gridi);
		int i=ij%gnx;
    	int j=ij/gnx;
        //get the central point of the grid cell
        double x=ax+(0.5+i)*gdx; 
        double y=ay+(0.5+j)*gdy;

        double min_size=max_len_dummy;
        for(int k=0;k<bdry_pt_ct;k++){
            double ex=bdry_pt_xy[2*k];
            double ey=bdry_pt_xy[2*k+1];
            double size=K*d_points(x,y,ex,ey)+bdry_pt_lfs[k];
            if(size<min_size){min_size=size;}
        }
        bgrid_sizing[gridi]=min_size;
	}

	//sizing grid for inner grids
	#pragma omp parallel for num_threads(num_t)
	for(int gridi=0;gridi<geo_igrid_ct;gridi++){
		int ij=geo_igrid_ij(gridi);
		int i=ij%gnx;
    	int j=ij/gnx;
        //get the central point of the grid cell
        double x=ax+(0.5+i)*gdx; 
        double y=ay+(0.5+j)*gdy;

        double min_size=max_len_dummy;
        for(int k=0;k<bdry_pt_ct;k++){
            double ex=bdry_pt_xy[2*k];
            double ey=bdry_pt_xy[2*k+1];
            double size=K*d_points(x,y,ex,ey)+bdry_pt_lfs[k];
            if(size<min_size){min_size=size;}
        }
        igrid_sizing[gridi]=min_size;
	}

}


/**
* @brief Calculates the density grid based on the sizing grid.
*/
void sizing_2d_automatic::get_density_grid(){

    if(K!=0){
    double t0=omp_get_wtime();
    	get_bdry_pts();
    double t1=omp_get_wtime(); printf("get bdry pt takes %g time\n",t1-t0);
        get_medial_pts();
    double t2=omp_get_wtime(); printf("get medial pt takes %g time\n",t2-t1);
        get_bdry_pt_lfs();
    double t3=omp_get_wtime(); printf("get bdry pt lfs takes %g time\n",t3-t2);
        get_sizing_grid();
    double t4=omp_get_wtime(); printf("get sizing grid takes %g time\n",t4-t3);

        //calculate density field from sizing field
        //using average smoothing over the region surrounding the grid. 
        //density grid for bdry grids
    	#pragma omp parallel for num_threads(num_t)
    	for(int gridi=0;gridi<geo_bgrid_ct;gridi++){
            double sizing_mean=0.0;
            int sizing_ct=0;
            int ij=geo_bgrid_ij(gridi);
            int i=ij%gnx;
            int j=ij/gnx;
            int il=i-density_sizing_smoothing_grid_ct; int ih=i+density_sizing_smoothing_grid_ct; 
            int jl=j-density_sizing_smoothing_grid_ct; int jh=j+density_sizing_smoothing_grid_ct;
            for(int ii=il;ii<=ih;ii++){
                for(int jj=jl;jj<=jh;jj++){
                    if(ii>=0 && ii<gnx && jj>=0 && jj<gny){
                        int iijj=gnx*jj+ii;
                        int grid_ctgr=geo_grid(iijj);
                        if(grid_ctgr<0){ //inner grid
                            sizing_mean+=igrid_sizing[abs(grid_ctgr)-1];
                            sizing_ct++;
                        }
                        if(grid_ctgr>0 && grid_ctgr<=gnxy){ //bdry grid
                            sizing_mean+=bgrid_sizing[grid_ctgr-1];
                            sizing_ct++;
                        }
                    }
                }
            }
            sizing_mean=sizing_mean/sizing_ct;
            bgrid_density[gridi]=1.0/pow(sizing_mean,density_sizing_exp);
    	}
    	//density grid for inner grids
    	#pragma omp parallel for num_threads(num_t)
    	for(int gridi=0;gridi<geo_igrid_ct;gridi++){
            double sizing_mean=0.0;
            int sizing_ct=0;
            int ij=geo_igrid_ij(gridi);
            int i=ij%gnx;
            int j=ij/gnx;
            int il=i-density_sizing_smoothing_grid_ct; int ih=i+density_sizing_smoothing_grid_ct; 
            int jl=j-density_sizing_smoothing_grid_ct; int jh=j+density_sizing_smoothing_grid_ct;
            for(int ii=il;ii<=ih;ii++){
                for(int jj=jl;jj<=jh;jj++){
                    if(ii>=0 && ii<gnx && jj>=0 && jj<gny){
                        int iijj=gnx*jj+ii;
                        int grid_ctgr=geo_grid(iijj);
                        if(grid_ctgr<0){ //inner grid
                            sizing_mean+=igrid_sizing[abs(grid_ctgr)-1];
                            sizing_ct++;
                        }
                        if(grid_ctgr>0 && grid_ctgr<=gnxy){ //bdry grid
                            sizing_mean+=bgrid_sizing[grid_ctgr-1];
                            sizing_ct++;
                        }
                    }
                }
            }
            sizing_mean=sizing_mean/sizing_ct;
            igrid_density[gridi]=1.0/pow(sizing_mean,density_sizing_exp);
    	}
    }
}


/**
 * @brief Calculates the sizing and density grids for the outer grids.
 * 
 * This function calculates the sizing and density grids for the outer grids.
 * Only called in CVD meshing.
 */
void sizing_2d_automatic::get_ogrid_sizing_density_grid(){
    if(K!=0){
        if(ogrid_sizing_density==false){

            ogrid_sizing_density=true;
            delete [] ogrid_sizing; delete [] ogrid_density;
            ogrid_sizing=new double[geo_ogrid_ct];
            ogrid_density=new double[geo_ogrid_ct];

            //calculate ogrid_sizing
            double max_len_dummy=(bx-ax)+(by-ay);

            //sizing grid for outer grids
            #pragma omp parallel for num_threads(num_t)
            for(int gridi=0;gridi<geo_ogrid_ct;gridi++){
                int ij=geo_ogrid_ij(gridi);
                int i=ij%gnx;
                int j=ij/gnx;
                //get the central point of the grid cell
                double x=ax+(0.5+i)*gdx; 
                double y=ay+(0.5+j)*gdy;

                double min_size=max_len_dummy;
                for(int k=0;k<bdry_pt_ct;k++){
                    double ex=bdry_pt_xy[2*k];
                    double ey=bdry_pt_xy[2*k+1];
                    double size=K*d_points(x,y,ex,ey)+bdry_pt_lfs[k];
                    if(size<min_size){min_size=size;}
                }
                ogrid_sizing[gridi]=min_size;
            }

            //density grid for outer grids
            #pragma omp parallel for num_threads(num_t)
            for(int gridi=0;gridi<geo_ogrid_ct;gridi++){
                double sizing_mean=0.0;
                int sizing_ct=0;
                int ij=geo_ogrid_ij(gridi);
                int i=ij%gnx;
                int j=ij/gnx;
                int il=i-density_sizing_smoothing_grid_ct; int ih=i+density_sizing_smoothing_grid_ct; 
                int jl=j-density_sizing_smoothing_grid_ct; int jh=j+density_sizing_smoothing_grid_ct;
                for(int ii=il;ii<=ih;ii++){
                    for(int jj=jl;jj<=jh;jj++){
                        if(ii>=0 && ii<gnx && jj>=0 && jj<gny){
                            int iijj=gnx*jj+ii;
                            int grid_ctgr=geo_grid(iijj);
                            if(grid_ctgr<0){ //inner grid
                                sizing_mean+=igrid_sizing[abs(grid_ctgr)-1];
                                sizing_ct++;
                            }
                            if(grid_ctgr>0 && grid_ctgr<=gnxy){ //bdry grid
                                sizing_mean+=bgrid_sizing[grid_ctgr-1];
                                sizing_ct++;
                            }
                            if(grid_ctgr>gnxy){ //outer grid
                                sizing_mean+=ogrid_sizing[grid_ctgr-gnxy-1];
                                sizing_ct++;
                            }
                        }
                    }
                }
                sizing_mean=sizing_mean/sizing_ct;
                ogrid_density[gridi]=1.0/pow(sizing_mean,density_sizing_exp);
            }

        }

    }

}

/**
 * @brief Calculates the sizing value for a given (x, y) coordinate in 2D.
 *
 * This function determines the sizing value for a specified (x, y) coordinate in 2D.
 * It checks the category of the grid where the point lies and returns the corresponding sizing value.
 * If the point lies in the inner grid, the sizing value from the `igrid_sizing` array is returned.
 * If the point lies in the boundary grid, the sizing value from the `bgrid_sizing` array is returned.
 * If the point lies in the outer grid, the sizing value from the `ogrid_sizing` array is returned.
 *
 * @param x The x-coordinate of the point.
 * @param y The y-coordinate of the point.
 * @return The sizing value for the specified (x, y) coordinate.
 */
double sizing_2d_automatic::getSizingVal(double x, double y){
    if(K==0){
        return 1.0;
    }
    else{
        //locate (x,y) in which cell of the density grid that it locates at. 
    	int grid_ctgr=shp->pt_geo_grid_val(x,y);
    	//if lie in inner grid
    	if(grid_ctgr<0){
    		int gridi=abs(grid_ctgr)-1;
    		return igrid_sizing[gridi];
    	}
    	//if lie in bdry grid
    	if(grid_ctgr>0 && grid_ctgr<=gnxy){
    		int gridi=grid_ctgr-1;
    		return bgrid_sizing[gridi];
    	}
        if(grid_ctgr>gnxy){ //in outer grid
            if(ogrid_sizing_density==false){
                //printf("Point in outer grid. ogrid sizing field unknown.\n");
                return 0;
            }
            else{
                int gridi=grid_ctgr-gnxy-1;
                return ogrid_sizing[gridi];
            }
        }
    	//printf("Point not in geometry grid! Sizing field unknown.\n");
    	return 0;
    }
}

/**
 * @brief Calculates the density value for a given (x, y) coordinate in 2D.
 *
 * This function determines the density value for a specified (x, y) coordinate in 2D.
 * It checks the category of the grid where the point lies and returns the corresponding density value.
 * If the point lies in the inner grid, the density value from the `igrid_density` array is returned.
 * If the point lies in the boundary grid, the density value from the `bgrid_density` array is returned.
 * If the point lies in the outer grid, the density value from the `ogrid_density` array is returned.
 *
 * @param x The x-coordinate of the point.
 * @param y The y-coordinate of the point.
 * @return The density value for the specified (x, y) coordinate.
 */
double sizing_2d_automatic::getDensityVal(double x, double y){
    if(K==0){
        return 1.0;
    }
    else{
        //locate (x,y) in which cell of the density grid that it locates at. 
    	int grid_ctgr=shp->pt_geo_grid_val(x,y);
    	//if lie in inner grid
    	if(grid_ctgr<0){
    		int gridi=abs(grid_ctgr)-1;
    		return igrid_density[gridi];
    	}
    	//if lie in bdry grid
    	if(grid_ctgr>0 && grid_ctgr<=gnxy){
    		int gridi=grid_ctgr-1;
    		return bgrid_density[gridi];
    	}
        if(grid_ctgr>gnxy){ //in outer grid
            if(ogrid_sizing_density==false){
                //printf("Point in outer grid. ogrid density field unknown.\n");
                return 0;
            }
            else{
                int gridi=grid_ctgr-gnxy-1;
                return ogrid_density[gridi];
            }
        }
        //printf("Point not in geometry grid! Density field unknown.\n");
    	return 0;
    }
}



/**
 * @brief Prints the sizing and density fields to files.
 *
 * This function writes the sizing and density fields to separate files.
 *
 * @param case_name The name of the case or directory where the files will be stored.
 */
void sizing_2d::print_fields_to_file(const char *case_name){
    char fs[256];
    sprintf(fs,"%s/sizing_field.txt",case_name);
    FILE *fsout=fopen(fs,"a");
    
    char fd[256];
    sprintf(fd,"%s/density_field.txt",case_name);
    FILE *fdout=fopen(fd,"a");
    
    for(int i=0;i<gnx;i++){
        for(int j=0;j<gny;j++){

        	double x=ax+(0.5+i)*gdx; 
        	double y=ay+(0.5+j)*gdy;

            fprintf(fsout,"%g ",getSizingVal(x,y));
            fprintf(fdout,"%g ",getDensityVal(x,y));
        }
        fprintf(fsout,"\n");
        fprintf(fdout,"\n");
    }
    fclose(fsout);
    fclose(fdout);
}


/**
 * @brief Prints the medial axis and boundary points to files.
 *
 * This function writes the medial and boundary points to separate files.
 *
 * @param case_name The name of the case or directory where the files will be stored.
 */
void sizing_2d_automatic::print_pts_to_file(const char *case_name){
    char fm[256];
    sprintf(fm,"%s/medial_pts.txt",case_name);
    FILE *fmout=fopen(fm,"a");
    
    char fe[256];
    sprintf(fe,"%s/bdry_pts.txt",case_name);
    FILE *feout=fopen(fe,"a");
    
    for(int i=0;i<bdry_pt_ct;i++){
        double x=bdry_pt_xy[2*i];
        double y=bdry_pt_xy[2*i+1];
        fprintf(feout,"%g %g \n",x,y);
    }
    for(int j=0;j<medial_pt_ct;j++){
        double mx=medial_pt_xy[2*j];
        double my=medial_pt_xy[2*j+1];
        fprintf(fmout,"%g %g \n",mx,my);
    }
    fclose(fmout);
    fclose(feout);
}




