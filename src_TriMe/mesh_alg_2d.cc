#include "mesh_alg_2d.hh"
#include "parallel_meshing_2d.hh"

using namespace voro;

/**
 * @brief Constructor for the mesh_alg_2d class.
 * @param pm2d_ A pointer to the parallel_meshing_2d object.
 */
mesh_alg_2d::mesh_alg_2d(parallel_meshing_2d *pm2d_)
	: pm2d(pm2d_), wis(wall_is_2d(pm2d_)),

	num_t(pm2d_->num_t), ax(pm2d_->ax), bx(pm2d_->bx), ay(pm2d_->ay), by(pm2d_->by),
	lx(pm2d_->lx), ly(pm2d_->ly), 

	gnx(pm2d_->gnx),gny(pm2d_->gny),gnxy(pm2d_->gnxy),
	gdx(pm2d_->gdx),gdy(pm2d_->gdy),diag_gdxy(pm2d_->diag_gdxy),
	inv_gdx(pm2d_->inv_gdx),inv_gdy(pm2d_->inv_gdy),
	geo_bgrid_ct(pm2d_->geo_bgrid_ct),geo_igrid_ct(pm2d_->geo_igrid_ct),
	Ntotal(pm2d_->Ntotal),Ncurrent(pm2d_->Ncurrent),Nremain(pm2d_->Nremain),Nfixed(pm2d_->Nfixed),

	xy_id_new(new double[2*(pm2d_->Ntotal)]),

	output_interval(pm2d_->output_interval),
	Continue(true),reTria(false),addPt(false),
	all_iter_ct(0),tria_iter_ct(0),inbtw_tria_iter_ct(0),
	stop_Continue_mvmt_thres(new double[1]),
	pt_mvmt_dis_thres(new double[1]),

	tria_ct(0),tria_vertex(new int[1]),tria_ccrd_h(new double [1]),
	bar_ct(0),barid(new int[1]),tria_centroid(new double[1]),
	tria_ar(new double[1]),tria_er(new double[1]),
	previous_alpha_mean_tria_ar(1.0), previous_alpha_mean_tria_er(1.0),

	Faces(new HE_face[1]),Vertices(new HE_vert[1]),HE_exist(false),

	tria_order_ccw(false)
	{
		tria_ct_private=new int[num_t];
		tria_vertex_private=new std::vector<int>[num_t];
		tria_centroid_private=new std::vector<double>[num_t];
		tria_ccrd_h_private=new std::vector<double>[num_t];
		tria_ar_private=new std::vector<double>[num_t];
		tria_er_private=new std::vector<double>[num_t];
		barid_private=new std::vector<int>[num_t];
		bar_ct_private=new int[num_t];
		max_tria_ar_private=new double[num_t];

		alpha_mean_arsum_private=new double[num_t];
		alpha_mean_ersum_private=new double[num_t];

		for(int i=0;i<num_t;i++){
			tria_ct_private[i]=0;
			bar_ct_private[i]=0;
		}
		tria_vertex_tct=0; tria_ccrd_h_tct=0; tria_centroid_tct=0;
		tria_ar_tct=0; tria_er_tct=0; barid_barct=0;
	}

/**
 * @brief Destructor for the mesh_alg_2d class.
 */
mesh_alg_2d::~mesh_alg_2d(){
	delete [] xy_id_new;

	delete [] tria_ct_private;
	delete [] tria_vertex_private;
	delete [] tria_centroid_private;
	delete [] tria_ccrd_h_private;
	delete [] tria_ar_private;
	delete [] tria_er_private;
	delete [] barid_private;
	delete [] bar_ct_private;
	delete [] max_tria_ar_private;

	delete [] alpha_mean_arsum_private;
	delete [] alpha_mean_ersum_private;


	delete [] stop_Continue_mvmt_thres;
	delete [] pt_mvmt_dis_thres;
	delete [] tria_vertex;
	delete [] tria_ccrd_h;
	delete [] tria_centroid;
	delete [] tria_ar;
	delete [] tria_er;
	delete [] barid;

	delete [] Faces;
	delete [] Vertices;
}

/**
 * @brief Changes the number of threads used in the mesh_alg_2d code.
 * @param new_num_t The new number of threads.
 */
void mesh_alg_2d::change_number_thread(int new_num_t){
	num_t=new_num_t;
	pm2d->con->change_number_thread(new_num_t);

	delete [] tria_ct_private;
	delete [] tria_vertex_private;
	delete [] tria_centroid_private;
	delete [] tria_ccrd_h_private;
	delete [] tria_ar_private;
	delete [] tria_er_private;
	delete [] barid_private;
	delete [] bar_ct_private;
	delete [] max_tria_ar_private;

	delete [] alpha_mean_arsum_private;
	delete [] alpha_mean_ersum_private;

	tria_ct_private=new int[new_num_t];
	tria_vertex_private=new std::vector<int>[new_num_t];
	tria_centroid_private=new std::vector<double>[new_num_t];
	tria_ccrd_h_private=new std::vector<double>[new_num_t];
	tria_ar_private=new std::vector<double>[new_num_t];
	tria_er_private=new std::vector<double>[new_num_t];
	barid_private=new std::vector<int>[new_num_t];
	bar_ct_private=new int[new_num_t];
	max_tria_ar_private=new double[new_num_t];

	alpha_mean_arsum_private=new double[new_num_t];
	alpha_mean_ersum_private=new double[new_num_t];

	for(int i=0;i<new_num_t;i++){
		tria_ct_private[i]=0;
		bar_ct_private[i]=0;
	}


}

/**
 * @brief Returns the category of a point identified by its ID.
 * @param id The ID of the point.
 * @return The category of the point.
 */
int mesh_alg_2d::pt_ctgr(int id){
	return pm2d->pt_ctgr[id];
}

/**
 * @brief Returns the characteristic edge length at a given location (x, y).
 * @param x The x-coordinate of the location.
 * @param y The y-coordinate of the location.
 * @return The characteristic edge length at the specified location.
 */
double mesh_alg_2d::get_chrtrt_len_h(double x, double y){
	return pm2d->get_chrtrt_len_h(x,y);
}

/**
 * @brief Returns the characteristic edge length at a given geometry grid with index ij.
 * @param ij The index of the geometry grid.
 * @return The characteristic length at the specified grid.
 */
double mesh_alg_2d::chrtrt_len_h(int ij){
	return pm2d->chrtrt_len_h[ij];
}

/**
 * Get the geometry boundary tolerance (geps) for a given point (x, y) lie in boundary grid.
 * @param x The x-coordinate of the point.
 * @param y The y-coordinate of the point.
 * @return The geometry boundary tolerance (geps) for the point.
 */
double mesh_alg_2d::get_bgrid_geps(double x, double y){
	return pm2d->get_bgrid_geps(x,y);
}

/**
 * Get the spatial step size (deps) for finite difference calculation in the point projection step,
 * for a given point (x, y) lie in boundatry grid.
 * @param x The x-coordinate of the point.
 * @param y The y-coordinate of the point.
 * @return The spatial step size (deps) for the point.
 * @throws std::exception if the point is not in the geometry boundary grid.
 */
double mesh_alg_2d::get_bgrid_deps(double x, double y){
	int ind=geo_grid(x,y)-1;
	if(ind<0 || ind>=gnxy){
		printf("ERROR 2: pt not in geometry bdry grid.\n");
    	throw std::exception();
    }
	return pm2d->bgrid_deps[ind];
}

/**
 * @brief Returns the category of the grid cell at the given index in the underlying geometry grid.
 * @param ij The index of the grid cell.
 * @return The status of the grid cell: -1, -2, -3 for inside; 1, 2, 3 for boundary; gnxy+1, gnxy+2, ... for outside.
 */
int mesh_alg_2d::geo_grid(int ij){
	return pm2d->geo_grid(ij);
}

/**
 * @brief Returns the category of the underlying geometry grid cell that the point (x, y) lies in.
 * @param x The x-coordinate of the point.
 * @param y The y-coordinate of the point.
 * @return The status of the grid cell: -1, -2, -3 for inside; 1, 2, 3 for boundary; gnxy+1, gnxy+2, ... for outside.
 */
int mesh_alg_2d::geo_grid(double x, double y){
	return pm2d->geo_grid(x,y);
}

/**
 * @brief Returns the current x-coordinate of a particle with the given id.
 * @param id The id of the particle.
 * @return The current x-coordinate of the particle.
 */
double mesh_alg_2d::current_x(int id){
	return pm2d->xy_id[2*id];
}

/**
 * @brief Returns the current y-coordinate of a particle with the given id.
 * @param id The id of the particle.
 * @return The current y-coordinate of the particle.
 */
double mesh_alg_2d::current_y(int id){
	return pm2d->xy_id[2*id+1];
}

/**
 * @brief Updates the x-coordinate of a particle with the given id to a new value.
 * @param id The id of the particle.
 * @param new_x The new x-coordinate value.
 */
void mesh_alg_2d::update_x(int id, double new_x){
	pm2d->xy_id[2*id]=new_x;
}

/**
 * @brief Updates the y-coordinate of a particle with the given id to a new value.
 * @param id The id of the particle.
 * @param new_y The new y-coordinate value.
 */
void mesh_alg_2d::update_y(int id, double new_y){
	pm2d->xy_id[2*id+1]=new_y;
}

/**
 * Calculates the geometry tolerance (geps) function for the given coordinates.
 *
 * @param x The x-coordinate of the point.
 * @param y The y-coordinate of the point.
 * @param geps_prime The default value of geps.
 * @return The calculated value of geps for the point.
 */
double mesh_alg_2d::geps_func(double x, double y, double geps_prime){
	return pm2d->geps_func(x,y,geps_prime);
}

/**
 * Calculates the finite differencing spatial step size (deps) function for the given coordinates.
 *
 * @param x The x-coordinate of the point.
 * @param y The y-coordinate of the point.
 * @param deps_prime The default value of deps.
 * @return The calculated value of deps for the point.
 */
double mesh_alg_2d::deps_func(double x, double y, double deps_prime){
	return pm2d->deps_func(x,y,deps_prime);
}

/**
 * @brief Computes the signed distance function (SDF) value at a given point.
 * @param x The x-coordinate of the point.
 * @param y The y-coordinate of the point.
 * @return The computed SDF value at the given point.
 */
double mesh_alg_2d::sdf_func(double x, double y){
	return pm2d->sdf_func(x,y);
}

/**
 * @brief Checks if a point is inside the outer grid.
 * @param x The x-coordinate of the point.
 * @param y The y-coordinate of the point.
 * @return True if the point is inside the outer grid, false otherwise.
 */
bool mesh_alg_2d::pt_in_outer_grid(double x,double y){
	return pm2d->pt_in_outer_grid(x,y);
}


/**
 * @brief Prints the particle coordinates to a file in a specific format.
 * @param file_name_prefix The file name prefix used for the output file.
 */
void mesh_alg_2d::print_particle_coords(const char *file_name_prefix){
	char bug0[256];
    sprintf(bug0,"%s_particle_coords.par",file_name_prefix);
    pm2d->con->draw_particles(bug0);
}


/**
 * @brief Prints the particle IDs and their corresponding coordinates to a text file.
 * @param file_name_prefix The file name prefix used for the output file.
 */
void mesh_alg_2d::print_xy_id(const char *file_name_prefix){
    char bug0[256];
    sprintf(bug0,"%s_xy_id.txt",file_name_prefix);
    FILE *outFile0 = fopen(bug0, "a");
    for(int pi=0;pi<Ncurrent;pi++){
        fprintf(outFile0,"%g %g \n",current_x(pi),current_y(pi));
    }
    fclose(outFile0);
}

/**
 * @brief Sorts the vertex IDs of triangles to Counter-Clockwise order
 * 
 * The order is sorted by calcuting the corss product of the two vectors.
 * If the order was A,B,C, calcuate the dot product of vector AB=B-A=(a1,a2) and BC=C-B=(b1,b2). 
 * If the cross product (0,0,a1b2-a2b1) has a1b2-a2b1>0, then A,B,C is the correct order. 
 * Otherwise, the CCW order is A,C,B.
 */
void mesh_alg_2d::sort_tria_vertex_ids_ccw(){
	#pragma omp parallel for num_threads(num_t)
	for(int triai=0;triai<tria_ct;triai++){
		int Aid=tria_vertex[3*triai];
		int Bid=tria_vertex[3*triai+1];
		int Cid=tria_vertex[3*triai+2];
		double Ax=current_x(Aid);
		double Ay=current_y(Aid);
		double Bx=current_x(Bid);
		double By=current_y(Bid);
		double Cx=current_x(Cid);
		double Cy=current_y(Cid);

		double a1=Bx-Ax; double a2=By-Ay;
		double b1=Cx-Bx; double b2=Cy-By;

		if(a1*b2-a2*b1<0){
			tria_vertex[3*triai+1]=Cid;
			tria_vertex[3*triai+2]=Bid;
		}
    }
    tria_order_ccw=true;
}


/**
 * @brief Prints the IDs of triangle bar vertices to a text file.
 * @param file_name_prefix The file name prefix used for the output file.
 */
void mesh_alg_2d::print_tria_bar_ids(const char *file_name_prefix){
	char bug0[256];
    sprintf(bug0,"%s_tria_bar_ids.txt",file_name_prefix);
    FILE *outFile0 = fopen(bug0, "a");

    if(HE_exist==true){
		//Loop through all edges, and print out the edge coords if edge ID (u,v) is u<v
		for(std::unordered_map<size_t, HE_edge*>::iterator itr_e=Edges.begin();itr_e!=Edges.end();itr_e++){
			int eu=itr_e->second->eu;
			int ev=itr_e->second->ev;
			if(eu<ev){
				fprintf(outFile0, "%d %d\n", eu,ev);
			}
		}
	}
	else{
	    for(int barii=0; barii<bar_ct; barii++){
	        int idd=barid[2*barii]; int nidd=barid[2*barii+1];
	        fprintf(outFile0, "%d %d\n", idd,nidd);
	    }
	}

    fclose(outFile0);
}

/**
 * @brief Prints the coordinates of triangle bar vertices to a text file.
 * @param file_name_prefix The file name prefix used for the output file.
 */
void mesh_alg_2d::print_tria_bar_coords(const char *file_name_prefix){
	char bug0[256];
    sprintf(bug0,"%s_tria_bar_coords.txt",file_name_prefix);
    FILE *outFile0 = fopen(bug0, "a");

	if(HE_exist==true){
		//Loop through all edges, and print out the edge coords if edge ID (u,v) is u<v
		for(std::unordered_map<size_t, HE_edge*>::iterator itr_e=Edges.begin();itr_e!=Edges.end();itr_e++){
			int eu=itr_e->second->eu;
			int ev=itr_e->second->ev;
			if(eu<ev){
				fprintf(outFile0, "%g %g \n%g %g\n\n", current_x(eu),current_y(eu), current_x(ev),current_y(ev));
			}
		}
	}
	else{
	    for(int barii=0; barii<bar_ct; barii++){
	        int idd=barid[2*barii]; int nidd=barid[2*barii+1];
	        fprintf(outFile0, "%g %g \n%g %g\n\n", current_x(idd),current_y(idd), current_x(nidd),current_y(nidd));
	    }
	}
	fclose(outFile0);
}


/**
 * @brief Prints the vertex IDs of triangles to a text file.
 * @param file_name_prefix The file name prefix used for the output file.
 */
void mesh_alg_2d::print_tria_vertex_ids(const char *file_name_prefix){
	if(tria_order_ccw==false){sort_tria_vertex_ids_ccw();}
	char bug0[256];
    sprintf(bug0,"%s_tria_vertex_ids.txt",file_name_prefix);
    FILE *outFile0 = fopen(bug0, "a");
    for(int triai=0;triai<tria_ct;triai++){
        fprintf(outFile0,"%d %d %d %d \n",triai,tria_vertex[3*triai],tria_vertex[3*triai+1],tria_vertex[3*triai+2]);
    }
    fclose(outFile0);
}

/**
 * @brief Prints the coordinates of triangle vertices to a text file.
 * @param file_name_prefix The file name prefix used for the output file.
 */
void mesh_alg_2d::print_tria_vertex_coords(const char *file_name_prefix){
	if(tria_order_ccw==false){sort_tria_vertex_ids_ccw();}
	char bug0[256];
    sprintf(bug0,"%s_tria_vertex_coords.txt",file_name_prefix);
    FILE *outFile0 = fopen(bug0, "a");
    for(int triai=0;triai<tria_ct;triai++){
    	int idd=tria_vertex[3*triai];
    	int nidd=tria_vertex[3*triai+1];
    	int nidd2=tria_vertex[3*triai+2];
        fprintf(outFile0,"%d %g %g %g %g %g %g \n",triai,current_x(idd),current_y(idd), current_x(nidd),current_y(nidd),current_x(nidd2),current_y(nidd2));
    }
    fclose(outFile0);
}

/**
 * @brief Prints the triangle quality statistics to text files.
 * @param file_name_prefix The file name prefix used for the output files.
 */
void mesh_alg_2d::print_tria_quality_stat(const char *file_name_prefix){
	char bug0[256];
    sprintf(bug0,"%s_tria_quality_stat_ar.txt",file_name_prefix);
    FILE *outFile0 = fopen(bug0, "a");

    char bug1[256];
    sprintf(bug1,"%s_tria_quality_stat_er.txt",file_name_prefix);
    FILE *outFile1 = fopen(bug1, "a");

    for(int triai=0;triai<tria_ct;triai++){
        fprintf(outFile0,"%g ",tria_ar[triai]);
        fprintf(outFile1,"%g ",tria_er[triai]);
    }

    fclose(outFile0);
    fclose(outFile1);
}

/**
 * @brief Prints the worst triangles (vertices IDs) to text files.
 * @param file_name_prefix The file name prefix used for the output files.
 */
void mesh_alg_2d::print_worst_tria(const char *file_name_prefix){
	char bug0[256];
    sprintf(bug0,"%s_worst_tria.txt",file_name_prefix);
    FILE *outFile0 = fopen(bug0, "a");

    for(int triai=0;triai<tria_ct;triai++){
		int idd=tria_vertex[3*triai];
    	int nidd=tria_vertex[3*triai+1];
    	int nidd2=tria_vertex[3*triai+2];

    	double tax=current_x(idd);
    	double tay=current_y(idd);
    	double tbx=current_x(nidd);
    	double tby=current_y(nidd);
    	double tcx=current_x(nidd2);
    	double tcy=current_y(nidd2);

    	double aspect_ratio=aspect_ratio_tria(tax,tay,tbx,tby,tcx,tcy);
        double edge_ratio=edge_ratio_tria(tax,tay,tbx,tby,tcx,tcy);
        if(aspect_ratio>2 || edge_ratio > 2){
        	fprintf(outFile0,"%d %d %d %d %g %g %g %g %g %g %g %g \n",triai, idd,nidd,nidd2, current_x(idd),current_y(idd), current_x(nidd),current_y(nidd),current_x(nidd2),current_y(nidd2),aspect_ratio,edge_ratio);

        }
    }

    fclose(outFile0);
}

/**
 * @brief Prints the overall statistics of triangle quality to a text file.
 */
void mesh_alg_2d::print_tria_quality_stat_overall(const char *file_name_prefix){
    char bug2[256];
    sprintf(bug2,"%s_tria_quality_stat_overall.txt",file_name_prefix);
    FILE *outFile2 = fopen(bug2, "a");

    fprintf(outFile2,"%d %d %g %g %g %g %g %g %g %g %g %g \n",
    	tria_iter_ct, tria_ct,
    	max_tria_ar,max_tria_er,
    	median_tria_ar, median_tria_er,
    	mean_tria_ar,mean_tria_er,
    	alpha_mean_tria_ar,alpha_mean_tria_er,
    	std_tria_ar,std_tria_er);

    fclose(outFile2);
}

/**
 * @brief Prints the Voronoi diagram to a GNUplot file.
 *
 * @param file_name_prefix The prefix for the output file name.
 */
void mesh_alg_2d::print_voro_diagram(const char *file_name_prefix){
	char bug0[256];
    sprintf(bug0,"%s_voro_diagram.gnu",file_name_prefix);
    pm2d->con->draw_cells_gnuplot(bug0);
}

/**
 * @brief Determines whether to print outputs based on the specified output interval and current iteration count.
 */
void mesh_alg_2d::determine_printOutputs(){
	printOutputs=false;
	//int pm2d->output_interval: 0, no output; -1, last final output; 10, every 10 triangulations output  
	if(pm2d->output_interval==-1){
		if(Continue==false || tria_iter_ct==1){
			printOutputs=true;
		}
	}
	if(pm2d->output_interval>0){
		if(Continue==false){
			printOutputs=true;
		}
		else{
			if(tria_iter_ct%pm2d->output_interval==0 || tria_iter_ct==1){
				printOutputs=true;
			}
		}
	}
}

/**
 * @brief Print out the boundary vertices in CCW order; v1, v2, ..., v10 : no wrap around at the end.
 * 		  Each line correspond to a boundary. 
 *
 */
void mesh_alg_2d::print_bdry_CCW(const char *file_name_prefix){

	for(int i=0;i<(signed int) Bdry_Edges_Start.size();i++){
		char bug0[256];
		sprintf(bug0,"%s_bdry_vertices_ids_CCW_%d.txt",pm2d->file_name_prefix,i);
    	FILE *outFile0 = fopen(bug0, "a");

    	char bug1[256];
		sprintf(bug1,"%s_bdry_vertices_coords_CCW_%d.txt",pm2d->file_name_prefix,i);
    	FILE *outFile1 = fopen(bug1, "a");

		bool Continue_bdry=true;
		std::pair<unsigned int, unsigned int> uv0=Bdry_Edges_Start[i];
		HE_edge* edge_i=Edges[KKey(uv0)];
		while(Continue_bdry){
			fprintf(outFile0, "%d \n", edge_i->vert->HE_vert_id);
			fprintf(outFile1, "%g %g \n", current_x(edge_i->vert->HE_vert_id),current_y(edge_i->vert->HE_vert_id));

			//CCW order, go in "prev" direction
			edge_i=edge_i->prev;
			if(edge_i==Edges[KKey(uv0)]){
				Continue_bdry=false;
			}
		}
		fclose(outFile0);
		fclose(outFile1);
	}

}


/**
 * @brief Prints various outputs related to particle coordinates, triangle bar coordinates, and Voronoi diagrams.
 */
void mesh_alg_2d::do_print_outputs(){
	char bug0[256];
    sprintf(bug0,"%s/%s_%d",pm2d->file_name_prefix,pm2d->file_name_prefix,tria_iter_ct);

	
	print_xy_id(bug0);   //Print out vertices ID's and coordinates 

	print_tria_bar_coords(bug0);  //Print out the unique triangulation bar coordinates
	print_tria_bar_ids(bug0);   //Print out the unique triangulation bar ID's
	print_tria_vertex_ids(bug0);  //Print out the triangle vertex ID's
	print_tria_vertex_coords(bug0);  //Print out the triangle vertex coordinates
	
	print_tria_quality_stat(bug0);  //Print out the statistics of each triangle's quality measures


	char bug1[256];
    sprintf(bug1,"%s/%s",pm2d->file_name_prefix,pm2d->file_name_prefix);
	print_tria_quality_stat_overall(bug1); //Print out the overall triangulation quality statistics

	//print_worst_tria(bug0);   //Print out the triangles with worst qualities
	//print_voro_diagram(bug0);   //Print out the associated Voronoi diagram
}

/**
 * @brief Prints various outputs related to particle coordinates, triangle bar coordinates, and Voronoi diagrams.
 * 		  For the final mesh: Implement Half-edge data strucutre, and print the boundary vertices in CCW order. 
 */
void mesh_alg_2d::do_print_final_outputs(){
	if(HE_exist==false){
		printf("ERROR: Half-edge data strucutre not constructed to print out data of the final mesh.\n");
    	throw std::exception();
    }

	char bug0[256];
    sprintf(bug0,"%s/%s_final",pm2d->file_name_prefix,pm2d->file_name_prefix);
	
	print_xy_id(bug0);   //Print out vertices ID's and coordinates 
	
	print_tria_bar_ids(bug0);   //Print out the unique triangulation bar ID's
	print_tria_bar_coords(bug0);  //Print out the unique triangulation bar coordinates
	print_tria_vertex_ids(bug0);  //Print out the triangle vertex ID's
	print_tria_vertex_coords(bug0);  //Print out the triangle vertex coordinates
	
	print_tria_quality_stat(bug0);  //Print out the statistics of each triangle's quality measures
	
	char bug1[256];
    sprintf(bug1,"%s/%s",pm2d->file_name_prefix,pm2d->file_name_prefix);
	print_tria_quality_stat_overall(bug1); //Print out the overall triangulation quality statistics

	print_bdry_CCW(bug0); //Print out boundary vertices in CCW order; Each line correspond to a separate bdry

}

/**
 * Checks whether a triangle defined by three vertices is a valid triangle within the shape.
 *
 * @param x0 The x-coordinate of the first vertex.
 * @param y0 The y-coordinate of the first vertex.
 * @param x1 The x-coordinate of the second vertex.
 * @param y1 The y-coordinate of the second vertex.
 * @param x2 The x-coordinate of the third vertex.
 * @param y2 The y-coordinate of the third vertex.
 * @return True if the triangle is valid, false otherwise.
 */
bool mesh_alg_2d::valid_tria(double x0,double y0,double x1,double y1,double x2,double y2){
	//test 1:Centroid
	double cen_x, cen_y;
	centroid_tria(x0,y0,x1,y1,x2,y2,cen_x,cen_y);
	//if centroid in inner grid, return true
	if(geo_grid(cen_x,cen_y)<0){return true;}
	//if centroid lie on outer grid, return false
	if(geo_grid(cen_x,cen_y)>gnxy){return false;}
	//if centroid on bdry grid: 
		//if centroid is not in/on geometry:fd(cx, cy) <=geps, return false;
	if(geo_grid(cen_x,cen_y)>0 && geo_grid(cen_x,cen_y)<=gnxy){
		if(pm2d->sdf_adf(cen_x,cen_y)>get_bgrid_geps(cen_x,cen_y)){return false;}
	}


	//test 2: after test 1, we know the remaining trias to test all have all vertices and centroid in/on geometry.
	//Next, we test if at least two edges midpoints <=geps, return true
	double e01x=0.5*(x0+x1);double e01y=0.5*(y0+y1);
	double e02x=0.5*(x0+x2);double e02y=0.5*(y0+y2);
	double e12x=0.5*(x1+x2);double e12y=0.5*(y1+y2);
	int edge_in_geo_ct=0;
	if(geo_grid(e01x,e01y)<0){edge_in_geo_ct++;}
	else if(geo_grid(e01x,e01y)>0 && geo_grid(e01x,e01y)<=gnxy){
		if(pm2d->sdf_adf(e01x,e01y)<=get_bgrid_geps(e01x,e01y)){
			edge_in_geo_ct++;
		}
	}

	if(geo_grid(e02x,e02y)<0){edge_in_geo_ct++;}
	else if(geo_grid(e02x,e02y)>0 && geo_grid(e02x,e02y)<=gnxy){
		if(pm2d->sdf_adf(e02x,e02y)<=get_bgrid_geps(e02x,e02y)){
			edge_in_geo_ct++;
		}
	}

	if(geo_grid(e12x,e12y)<0){edge_in_geo_ct++;}
	else if(geo_grid(e12x,e12y)>0 && geo_grid(e12x,e12y)<=gnxy){
		if(pm2d->sdf_adf(e12x,e12y)<=get_bgrid_geps(e12x,e12y)){
			edge_in_geo_ct++;
		}
	}
	if(edge_in_geo_ct>=2){
		return true;
	}
	
	//test 4: 
		//find circumcenter of triangle
			//if circumcenter in outer grid, return false
			//else
				//if adf_sdf(circumcenter)<0: return true
				//if adf_sdf(circumcenter)/circumradius<thres(0.4), return true
				//else: return false
	//printf("still need to test triangle: use circumcenter and circumradius testing \n");
	double ccc_x,ccc_y;
	circumcenter_tria(x0,y0,x1,y1,x2,y2,ccc_x,ccc_y);
	double sdf_ccc=pm2d->sdf_adf(ccc_x,ccc_y);
	double ccr=circumradius_tria(x0,y0,x1,y1,x2,y2);
	if(sdf_ccc/ccr<tria_ccc_cri_fac){return true;}
	return false;
	
}

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
//extract tria info for a tria by vertex ids: vid1, vid2, vid3
bool mesh_alg_2d::extract_tria_info(int vid1,int vid2, int vid3,
    
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

	bool cvd_pt_update, bool &use_random_point, double &random_point_x, double &random_point_y,  unsigned int &seed_i){

	int vid1_2=2*vid1; int vid2_2=2*vid2; int vid3_2=2*vid3;
	double tax=pm2d->xy_id[vid1_2];double tay=pm2d->xy_id[vid1_2+1];
    double tbx=pm2d->xy_id[vid2_2];double tby=pm2d->xy_id[vid2_2+1];
    double tcx=pm2d->xy_id[vid3_2];double tcy=pm2d->xy_id[vid3_2+1];
    bool tria_is_valid=false;

    double ar_alpha_mean_rel_change=fabs(alpha_mean_tria_ar-previous_alpha_mean_tria_ar)/previous_alpha_mean_tria_ar;

    if(valid_tria(tax,tay,tbx,tby,tcx,tcy)==true){
        tria_is_valid=true;

        double edge_ratio=0.0;
        double aspect_ratio=0.0;

        
        if(cvd_pt_update==true && Ncurrent==Ntotal && ar_alpha_mean_rel_change<0.002 && tria_iter_ct%5==0){
        	double a=d_points(tax,tay,tbx,tby);
			double b=d_points(tbx,tby,tcx,tcy);
			double c=d_points(tax,tay,tcx,tcy);
        	edge_ratio=edge_ratio_tria(a,b,c);
        	double s=s_tria(a,b,c);
		    aspect_ratio=aspect_ratio_tria(a,b,c,s);
    		if(edge_ratio>4.0){
        		//a (vid1,vid2) is the smallest length edge 
        		if(a<b && a<c){
        			//move the current vertex vid1 to a random point in one of the nearby inner grid
        			if(vid1<vid2){
                		use_random_point=true;
        			}
        			if(vid2<Nfixed){
        				use_random_point=true;
        			}
        		}
        		//c (vid1,vid3) is the smallest length edge
        		if(c<a && c<b){
        			//move the current vertex vid1 to a random point in one of the nearby inner grid
        			if(vid1<vid3){
        				use_random_point=true;
        			}
        			if(vid3<Nfixed){
        				use_random_point=true;
        			}
        		}
        	}
        	
        	else if(aspect_ratio>6.0){
        		double tria_centroid_x, tria_centroid_y;
        		centroid_tria(tax,tay,tbx,tby,tcx,tcy,tria_centroid_x,tria_centroid_y);
        		double ac=d_points(tax,tay,tria_centroid_x,tria_centroid_y);
        		double bc=d_points(tbx,tby,tria_centroid_x,tria_centroid_y);
        		double cc=d_points(tcx,tcy,tria_centroid_x,tria_centroid_y);

        		if(ac<bc && ac<cc){
        			use_random_point=true;
        		}
        	}
			

        	if(use_random_point==true){
        		//set values for random_point_x, random_point_y
        		int incre=2;
				int pgridi=(tax-ax)*inv_gdx;
            	int pgridj=(tay-ay)*inv_gdy;
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
            							inner_grid_found.push_back(gridi);
            							inner_grid_found.push_back(gridj);
            							inner_grid_found_ct++;
            							continue_search=false;
            						}
            					}
            				}
            			}
            		}
            		incre++;
            	}

            	//Select a random inner grid in the outer layer found
            	int random_igrid_found_i2=2* (rand_r(&seed_i) % inner_grid_found_ct); //random int from 0 to inner_grid_found_ct-1
            	int gridi=inner_grid_found[random_igrid_found_i2];
            	int gridj=inner_grid_found[random_igrid_found_i2+1];
            	random_point_x=ax+(rnd_r(&seed_i)+gridi)*gdx; 
            	random_point_y=ay+(rnd_r(&seed_i)+gridj)*gdy;
        	}
        }

        if(vid3>vid2 && vid2>vid1){ //to store triangles related quantities uniquely. Only count triangles in vertex ascending order, vid1<vid2<vid3

            if(store_tria_vertex){
                if(tria_ct_p<tria_vertex_p_tct){
                    int store_ind_3=tria_ct_p*3;
                    tria_vertex_p[store_ind_3]=vid1;
                    tria_vertex_p[store_ind_3+1]=vid2;
                    tria_vertex_p[store_ind_3+2]=vid3;
                }
                else{
                    tria_vertex_p.push_back(vid1);
                    tria_vertex_p.push_back(vid2);
                    tria_vertex_p.push_back(vid3);
                }
                
            }
            
            if(store_tria_centroid || store_ccrd_h_ratio){
                double centroid_x; double centroid_y;
                centroid_tria(tax,tay,tbx,tby,tcx,tcy,centroid_x,centroid_y);

                if(store_tria_centroid){
                    if(tria_ct_p<tria_centroid_p_tct){
                        int store_ind_2=tria_ct_p*2;
                        tria_centroid_p[store_ind_2]=centroid_x;
                        tria_centroid_p[store_ind_2+1]=centroid_y;
                    }
                    else{
                        tria_centroid_p.push_back(centroid_x);
                        tria_centroid_p.push_back(centroid_y);
                    }
                    
                }

                if(store_ccrd_h_ratio){
                    double circumradius=circumradius_tria(tax,tay,tbx,tby,tcx,tcy);
                    double h_i=get_chrtrt_len_h(centroid_x,centroid_y);
                    if(tria_ct_p<tria_ccrd_h_p_tct){
                        tria_ccrd_h_p[tria_ct_p]=circumradius/h_i;
                    }
                    else{
                        tria_ccrd_h_p.push_back(circumradius/h_i);
                    }
                }
            }
        
            if(cvd_pt_update==false || Ncurrent!=Ntotal || ar_alpha_mean_rel_change>=0.002 || tria_iter_ct%5!=0){
            	edge_ratio=edge_ratio_tria(tax,tay,tbx,tby,tcx,tcy);
            	aspect_ratio=aspect_ratio_tria(tax,tay,tbx,tby,tcx,tcy);
            }
            alpha_mean_arsum_p+=sqrt(aspect_ratio);
            alpha_mean_ersum_p+=sqrt(edge_ratio);

            if(store_triangle_quality){
                if(tria_ct_p<tria_ar_p_tct){
                    tria_ar_p[tria_ct_p]=aspect_ratio;
                }
                else{
                    tria_ar_p.push_back(aspect_ratio);
                }

                if(tria_ct_p<tria_er_p_tct){
                    tria_er_p[tria_ct_p]=edge_ratio;
                }
                else{
                    tria_er_p.push_back(edge_ratio);
                }
            }
            if(store_max_tria_ar){
            	if(aspect_ratio>max_tria_ar_){max_tria_ar_=aspect_ratio;}
            }

            if(do_tria_ct){
                tria_ct_p++;
            }


        }
    }
    return tria_is_valid;
}



/**
 * @brief Checks if the Voronoi neighbors of a particle satisfy a certain requirement: 
 * It must have two consecutive neighbors that can form a potential triangle abc. 
 *
 * @param[out] recompute Set to true if the Voronoi neighbors need to be recomputed, false otherwise.
 * @param[in] neigh Vector of neighboring particles.
 * @param[in] x X-coordinate of the particle.
 * @param[in] y Y-coordinate of the particle.
 */
void mesh_alg_2d::check_voro_nei_particle_requirement(bool &recompute,std::vector<int> neigh,double x,double y)
{
	int i=0;
    int neigh_size=(signed int) neigh.size();
    while(recompute==true && i<neigh_size){
    	if (neigh[i]>=0){ //a valid particle neighbor

        	//the index of the next consecutive neighbor c, which construct a potential triangle abc
            int nei2n;    
            //check counterclockwise side
            if(i!=neigh_size-1){nei2n=i+1;}
            else {nei2n=0;}
            if(neigh[nei2n]>=0){
            	recompute=false;
            }

            if(recompute==true){
            	//check clockwise side if there is a valid particle neighbor
                if(i!=0){nei2n=i-1;}
                else {nei2n=neigh_size-1;}
                if(neigh[nei2n]>=0){
                	recompute=false;
                }
            }
            
        }
    	i++;
    }
     //if not satisfy, increase the wall_2d wis local cut off length, and recompute the Voronoi cell
    if(recompute==true){
    	//find the geometry grid the pt lies in
	    int gi=int((x-ax)*inv_gdx); 
	    int gj=int((y-ay)*inv_gdy); 
	    int gij=gnx*gj+gi;
        #pragma omp critical
        {
		    //update wall_2d to be larger
		    wis.update_tria_length_cri(gij,2.0);
		}
    }
}

/**
 * Computes and stores various information related to Voronoi cells and triangles.
 *
 * @param store_tria_bar Flag indicating whether to store triangle bar coordinates.
 * @param store_triangle_quality Flag indicating whether to store triangle quality information.
 * @param store_tria_vertex Flag indicating whether to store triangle vertex information.
 * @param store_tria_centroid Flag indicating whether to store triangle centroid information.
 * @param store_ccrd_h_ratio Flag indicating whether to store circumradius-to-characteristic-edge ratio for triangles.
 * @param store_max_tria_ar Flag indicating whether to store maximum triangle aspect ratio.
 * @param cvd_pt_update Flag indicating whether to update point positions using centroidal Voronoi diagram (CVD) meshing.
 * @param xy_id_store_array Pointer to the array for storing point x-y coordinates and IDs.
 */
void mesh_alg_2d::voro_compute_and_store_info(bool store_tria_bar,bool store_triangle_quality,
	bool store_tria_vertex,
	bool store_tria_centroid,
    bool store_ccrd_h_ratio, bool store_max_tria_ar,
    bool cvd_pt_update,
    bool store_xy_id, double *&xy_id_store_array){

	if(store_max_tria_ar){
		max_tria_ar_prev=max_tria_ar;
		max_tria_ar=-(lx+ly); 
	}

	//clear container and update points in container.
    pm2d->con->clear();
    pm2d->con->add_parallel(pm2d->xy_id, Ncurrent, num_t);
    pm2d->con->put_reconcile_overflow();

	//Voronoi computation
	//triangle extraction: 
		//for each triangle, calculate its circumradius. 
		//for each triangle, calculate its circumradius/hi ratio, store in array tria_ccrd and tria_ccrd_h
	bool do_tria_ct=true;

	container_2d::iterator cli;

	#pragma omp parallel num_threads(num_t)
	{   
	//thread-private variables: avoid memory rushing issue
	#ifdef _OPENMP
	    const int ithread = omp_get_thread_num();
	#else
	    const int ithread=0;
	#endif
	    voronoicell_neighbor_2d c(*(pm2d->con));

	    unsigned int seed_i=omp_get_thread_num()+100;

	    alpha_mean_arsum_private[ithread]=0.0;
	    alpha_mean_ersum_private[ithread]=0.0;

	    if(store_max_tria_ar){max_tria_ar_private[ithread]=-(lx+ly);}

	    int barid_p_barct=(signed int) barid_private[ithread].size()/2;
	    int tria_vertex_p_tct=(signed int) tria_vertex_private[ithread].size()/3;
	    int tria_centroid_p_tct=(signed int) tria_centroid_private[ithread].size()/2;
	    int tria_ccrd_h_p_tct=(signed int) tria_ccrd_h_private[ithread].size();
	    int tria_ar_p_tct=(signed int) tria_ar_private[ithread].size();
	    int tria_er_p_tct=(signed int) tria_er_private[ithread].size();

	    tria_ct_private[ithread]=0;
	    bar_ct_private[ithread]=0;
	    
	    #pragma omp for schedule(dynamic) reduction(+:inner_pt_over_stop_Continue_mvmt_thres_ct)
	    for(cli=pm2d->con->begin(); cli<pm2d->con->end(); cli++){

	    	int ij=(*cli).ijk;
			int q=(*cli).q;
	    	int id=pm2d->con->id[ij][q];
	    	double *pp=pm2d->con->p[ij]+pm2d->con->ps*q;
	    	double x=*(pp++);double y=*pp;

	    	if(store_xy_id==true){
	    		int id2=2*id;
				xy_id_store_array[id2]=x;
				xy_id_store_array[id2+1]=y;
	    	}

	    	//neighbors requirement: the Voronoi cell should have at least two consecutive particle neighbors. 
	    	//If not, increase the wall_2d wis local cut off length, and recompute the Voronoi cell
	    	//Until the requirement is satisfied. 
	    	bool recompute=true;
	    	while(recompute){
		        if(pm2d->con->compute_cell(c, cli)) {
			        std::vector<int> neigh;   //neighbors of a node

			        //check if Voronoi cell satisfy neighbors requirement
			        c.neighbors_sorted(neigh); 
			        check_voro_nei_particle_requirement(recompute,neigh,x,y);

			        //Voronoi cell satisfy neighbor requirement
			        if(recompute==false){

			        	bool use_random_point=false;
			        	double random_point_x=0.0;
			        	double random_point_y=0.0;

				        //Collect triangles info if needed; And calculate alpha_mean
				        for(int i=0;i<(signed int) neigh.size();i++){
				            if (neigh[i]>=0 && neigh[i]>id){ //store bars uniquely: barid2>barid1

				            	bool barid_store=false;

				            	//the index of the second neighbor c, which construct the triangle abc
				                int nei2n;    
		                        //check counterclockwise side
		                        if(i!=(signed int)neigh.size()-1){nei2n=i+1;}
		                        else {nei2n=0;}

		                        if(neigh[nei2n]>=0){
		                        	if(extract_tria_info(id,neigh[i],neigh[nei2n],

				            		store_tria_vertex,do_tria_ct,
									store_tria_centroid,store_ccrd_h_ratio,
									store_triangle_quality,store_max_tria_ar,

									tria_vertex_private[ithread],
									tria_centroid_private[ithread],tria_ccrd_h_private[ithread],
									tria_ar_private[ithread],tria_er_private[ithread],

									tria_ct_private[ithread],

									alpha_mean_arsum_private[ithread],alpha_mean_ersum_private[ithread],

									tria_vertex_p_tct,
    								tria_centroid_p_tct,
    								tria_ccrd_h_p_tct,
    								tria_ar_p_tct,
    								tria_er_p_tct, 

    								max_tria_ar_private[ithread],

    								cvd_pt_update, use_random_point,random_point_x,random_point_y, seed_i

    								)){barid_store=true;}
		                        }
				            	
		                        //check clockwise side if is a valid triangle
		                        if(i!=0){nei2n=i-1;}
		                        else {nei2n=(signed int)neigh.size()-1;}

		                        if(neigh[nei2n]>=0){
			                        if(extract_tria_info(id,neigh[i],neigh[nei2n],

					            		store_tria_vertex,do_tria_ct,
										store_tria_centroid,store_ccrd_h_ratio,
										store_triangle_quality,store_max_tria_ar,

										tria_vertex_private[ithread],
										tria_centroid_private[ithread],tria_ccrd_h_private[ithread],
										tria_ar_private[ithread],tria_er_private[ithread],

										tria_ct_private[ithread],

										alpha_mean_arsum_private[ithread],alpha_mean_ersum_private[ithread],
										
										tria_vertex_p_tct,
    									tria_centroid_p_tct,
    									tria_ccrd_h_p_tct,
    									tria_ar_p_tct,
    									tria_er_p_tct, 

    									max_tria_ar_private[ithread],

    									cvd_pt_update, use_random_point,random_point_x,random_point_y, seed_i

    									)){barid_store=true;}
			                    }
			                
		                        //bar is a side of a valid triangle, store to barid, increase bar_ct
								if(store_tria_bar){
									if(barid_store==true){

										if(bar_ct_private[ithread]<barid_p_barct){
											int store_ind_2=2*bar_ct_private[ithread];
											barid_private[ithread][store_ind_2]=id;
											barid_private[ithread][store_ind_2+1]=neigh[i];
										}
										else{
											barid_private[ithread].push_back(id);
											barid_private[ithread].push_back(neigh[i]);
										}
										
										bar_ct_private[ithread]++;

									}
								}
		                    //exiting if (neigh[i]>=0 && neigh[i]>id)
		                    }
		                //exiting loop through pt's neighbors
		                }

		                if(cvd_pt_update==true){
				        	//CVD: Compute centroid and new point position here
				        	//and store the new pt in xy_id
				        	update_pt_position_cvd(id,c,x,y,inner_pt_over_stop_Continue_mvmt_thres_ct,
				        		use_random_point,random_point_x,random_point_y);
				        }

		            //exit if(recompute==false), collect Voro and tria info
			        } 
			    //exit if Voronoi cell is true
			    }
			    else{ //No Voronoi cell for the particle. 
			    	recompute=false;
			    }
			//exit while(recompute)
			}
        //exit looping over particles
	    }
	}

	//Half-edge data structure needs to be re-calculated
	HE_exist=false;

	//obtain full set of triangle ct, triangle vertex, triangle circumradius, triangle circumradius/h_i	
	bar_ct=0;
	tria_ct=0;
	//alpha_mean_ar and alpha_mean_er computation
	double alpha_mean_arsum=0.0; double alpha_mean_ersum=0.0;
    for(int i=0;i<num_t;i++){
        tria_ct+=tria_ct_private[i];
        bar_ct+=bar_ct_private[i];
        alpha_mean_arsum+=alpha_mean_arsum_private[i];
        alpha_mean_ersum+=alpha_mean_ersum_private[i];
        if(store_max_tria_ar){
        	if(max_tria_ar<max_tria_ar_private[i]){max_tria_ar=max_tria_ar_private[i];}
        }
    }

    double mean_ar_er_fac=1.0/(1.0*tria_ct);
    alpha_mean_tria_ar=alpha_mean_arsum*mean_ar_er_fac;
	alpha_mean_tria_er=alpha_mean_ersum*mean_ar_er_fac;
	alpha_mean_tria_ar=alpha_mean_tria_ar*alpha_mean_tria_ar;
	alpha_mean_tria_er=alpha_mean_tria_er*alpha_mean_tria_er;

	if(store_tria_vertex){
		if(tria_ct>tria_vertex_tct){
			delete [] tria_vertex;
    		tria_vertex=new int[3*tria_ct];
    		tria_vertex_tct=tria_ct;
		}
    	
    }
    if(store_tria_centroid){
    	if(tria_ct>tria_centroid_tct){
	    	delete [] tria_centroid;
	    	tria_centroid=new double[2*tria_ct];
	    	tria_centroid_tct=tria_ct;
	    }
    }
    if(store_ccrd_h_ratio){
    	if(tria_ct>tria_ccrd_h_tct){
	    	delete [] tria_ccrd_h;
	    	tria_ccrd_h=new double[tria_ct];
	    	tria_ccrd_h_tct=tria_ct;
	    }
    }

    if(store_triangle_quality){
    	if(tria_ct>tria_ar_tct){
    		delete [] tria_ar;
    		tria_ar=new double[tria_ct];
    		tria_ar_tct=tria_ct;
    	}
    	if(tria_ct>tria_er_tct){
    		delete [] tria_er;
    		tria_er=new double[tria_ct];
    		tria_er_tct=tria_ct;
    	}
    }
	
    
    if(store_tria_bar){
    	if(bar_ct>barid_barct){
    		delete [] barid;
			barid=new int[2*bar_ct];
			barid_barct=bar_ct;
    	}
    }


    #pragma omp parallel num_threads(num_t)
    {
    //thread-private variables: avoid memory rushing issue
	#ifdef _OPENMP
	    const int i = omp_get_thread_num();
	#else
	    const int i=0;
	#endif

	    //triangles: collect vertex id from each thread to tri_vertex
        if(store_tria_vertex || store_ccrd_h_ratio || store_triangle_quality){
        	int ind_offset_tria=0;
	        if(i!=0){
	            for(int j=0;j<=i-1;j++){
	                ind_offset_tria+=tria_ct_private[j];
	            }
	        }
	        for(int ii=0;ii<tria_ct_private[i];ii++){
	        	if(store_tria_vertex){
	        		int indp_tri=3*ind_offset_tria;
	            	int ii3=3*ii;
		            tria_vertex[ii3+indp_tri]=tria_vertex_private[i][ii3];
		            tria_vertex[ii3+indp_tri+1]=tria_vertex_private[i][ii3+1];
		            tria_vertex[ii3+indp_tri+2]=tria_vertex_private[i][ii3+2];
	        	}
	        	if(store_ccrd_h_ratio){
	        		tria_ccrd_h[ii+ind_offset_tria]=tria_ccrd_h_private[i][ii];
	        	}
	        	if(store_triangle_quality){
	        		tria_ar[ii+ind_offset_tria]=tria_ar_private[i][ii];
	        		tria_er[ii+ind_offset_tria]=tria_er_private[i][ii];
	        	}
	        	if(store_tria_centroid){
	        		int indp_twi=2*ind_offset_tria;
	        		int ii2=2*ii;
	            	tria_centroid[ii2+indp_twi]=tria_centroid_private[i][ii2];
	            	tria_centroid[ii2+indp_twi+1]=tria_centroid_private[i][ii2+1];
	        	}
	        }
	        if(store_tria_vertex){tria_order_ccw=false;}
	    }

        if(store_tria_bar){
	        //edge barid: 
			int ind_offset_bar=0;
		    if(i!=0){
		        for(int j=0;j<=i-1;j++){
		            ind_offset_bar+=bar_ct_private[j];
		        }
		    }
		    int indp_bar=2*(ind_offset_bar);
		    for(int ii=0;ii<bar_ct_private[i];ii++){
		        int ii2=2*ii;
		        barid[ii2+indp_bar]=barid_private[i][ii2];
		        barid[ii2+indp_bar+1]=barid_private[i][ii2+1];
		    }
		}
    }

	//triangle quality statistics
	if(store_triangle_quality){
		max_tria_ar=-(lx+ly); 
   		max_tria_er=-(lx+ly); 
   		double arsum=0.0; double arsum2=0.0; //aspect ratio
    	double ersum=0.0; double ersum2=0.0;  //edge ratio
		for(int i=0;i<tria_ct;i++){
			arsum+=tria_ar[i]; arsum2+=sqr(tria_ar[i]);
            ersum+=tria_er[i]; ersum2+=sqr(tria_er[i]);
            if(tria_ar[i]>max_tria_ar){max_tria_ar=tria_ar[i];}
            if(tria_er[i]>max_tria_er){max_tria_er=tria_er[i];}
		}

		std::sort(tria_ar, tria_ar + tria_ct); 
    	std::sort(tria_er, tria_er + tria_ct); 
		
        median_tria_ar=tria_ar[tria_ct/2];
        median_tria_er=tria_er[tria_ct/2];

        mean_tria_ar=arsum/tria_ct;
        std_tria_ar=sqrt(arsum2/tria_ct-sqr(mean_tria_ar));
        mean_tria_er=ersum/tria_ct;
        std_tria_er=sqrt(ersum2/tria_ct-sqr(mean_tria_er));
    }
}

/**
 * Determines whether to continue triangulation based on triangle quality metrics.
 */
void mesh_alg_2d::determine_Continue_tria_quality(){

	double er_alpha_mean_rel_change=fabs(alpha_mean_tria_er-previous_alpha_mean_tria_er)/previous_alpha_mean_tria_er;
    double ar_alpha_mean_rel_change=fabs(alpha_mean_tria_ar-previous_alpha_mean_tria_ar)/previous_alpha_mean_tria_ar;

    if(ar_alpha_mean_rel_change<stop_Continue_ar_mean_rel_change 
    	&& 
    	er_alpha_mean_rel_change<stop_Continue_er_mean_rel_change){
    	Continue=false; 
    }
    hybrid_determine_method_switch(ar_alpha_mean_rel_change);
    
    previous_alpha_mean_tria_er=alpha_mean_tria_er; 
    previous_alpha_mean_tria_ar=alpha_mean_tria_ar; 

    if(max_tria_ar_prev>0){
        double max_tria_ar_rel_change=fabs(max_tria_ar-max_tria_ar_prev)/max_tria_ar_prev;
        if(max_tria_ar_rel_change>stop_Continue_max_ar_rel_change && max_tria_ar>3.0){
            Continue=true;
        }
    }
    if(Continue==false){printf("\nquality stop \n");}
}

        
/**
 * Adds new points to the mesh.
 */
void mesh_alg_2d::add_new_pts(){
printf("start add pt \n");
	int Nadd=add_pt_fac*Ncurrent;
	if(Nadd>Nremain){
		Nadd=Nremain;
	}
	if(Nadd>tria_ct){Nadd=tria_ct;}

//---------------------------------------------------------------------

	//Voronoi computation
		//triangle extraction: 
			//for each triangle, calculate its circumradius. 
			//for each triangle, calculate its circumradius/hi ratio, store in array tria_ccrd and tria_ccrd_h
	bool store_barid=false;
	bool store_triangle_quality=false;
	bool store_tria_vertex=false;
    bool store_tria_centroid=true;
    bool store_ccrd_h_ratio=true;
    bool store_max_tria_ar=false;
    bool cvd_pt_update=false;
    bool store_xy_id=false;
    double *dummy_array=new double[1];
	voro_compute_and_store_info(store_barid,store_triangle_quality,
		store_tria_vertex,
		store_tria_centroid,
		store_ccrd_h_ratio, store_max_tria_ar,
		cvd_pt_update,
		store_xy_id,dummy_array);
	delete [] dummy_array;

//---------------------------------------------------------------------
	//sort triangles index based on ccrd/h, in descending order
    std::vector<int> ccrd_h_sorted_ind(tria_ct);
    #pragma omp parallel for num_threads(num_t) 
    for(int i=0;i<tria_ct;i++){
    	ccrd_h_sorted_ind[i]=i;
    }
    std::sort(  std::begin(ccrd_h_sorted_ind), 
                std::end(ccrd_h_sorted_ind),
                [&](int i1, int i2) { return tria_ccrd_h[i1] > tria_ccrd_h[i2]; } );
	//add points: Add to the centroids of the triangles with largest ccrd/h_i ratios
double t0=omp_get_wtime();
    #pragma omp parallel for num_threads(num_t) schedule(guided) 
    for(int i=0;i<Nadd;i++){
    	int ind=ccrd_h_sorted_ind[i];
    	int ind2=2*ind;
    	double new_x=tria_centroid[ind2];
    	double new_y=tria_centroid[ind2+1];

    	int m=Ncurrent+i;
    	int m2=2*m;
    	pm2d->xy_id[m2]=new_x;
    	pm2d->xy_id[m2+1]=new_y;

    }
t_add_pt_centroid+=omp_get_wtime()-t0;


	Nremain-=Nadd;
	int Ncurrent_old=Ncurrent;
	Ncurrent+=Nadd;

	pm2d->Ncurrent=Ncurrent;
	pm2d->Nremain=Nremain;
//---------------------------------------------------------------------
    //update adaptive geometric arrays
    double adaptive_scale_fac=1.0/(sqrt(1.0*(Ncurrent/Ncurrent_old)));
    pm2d->chrtrt_len_h_avg*=adaptive_scale_fac;
	#pragma omp parallel num_threads(num_t)
	{
		#pragma omp for
		for(int gi=0;gi<geo_bgrid_ct;gi++){
			int ij=pm2d->geo_bgrid_ij(gi);
			pm2d->chrtrt_len_h[ij]*=adaptive_scale_fac;
			pm2d->bgrid_geps[gi]*=adaptive_scale_fac;
			pm2d->bgrid_deps[gi]*=adaptive_scale_fac;
			update_retria_movement_thres(ij);
			update_stop_Continue_mvmt_thres(ij);
			update_pt_mvmt_dis_thres(ij);
		}
		#pragma omp for
		for(int gi=0;gi<geo_igrid_ct;gi++){
			int ij=pm2d->geo_igrid_ij(gi);
			pm2d->chrtrt_len_h[ij]*=adaptive_scale_fac;
			update_retria_movement_thres(ij);
			update_stop_Continue_mvmt_thres(ij);
			update_pt_mvmt_dis_thres(ij);
		}
	}
	wis.update_tria_length_cri();

//---------------------------------------------------------------------
    //loop through points against new geps, 
    //re-ctgr them to pt_ctgr
    //and project any outside points back to the new geps
    inner_pt_ct_temp=0;   
    #pragma omp parallel for num_threads(num_t) schedule(dynamic) reduction(+: inner_pt_ct_temp)
	for(int i=0;i<Ncurrent;i++){
		int i2=2*i;
		double x=pm2d->xy_id[i2];
		double y=pm2d->xy_id[i2+1];
		int g_ind0=geo_grid(x,y);
    	if(g_ind0<0){ //inner grid
    		pm2d->pt_ctgr[i]=-1; //-1: inner pt
    		inner_pt_ct_temp+=1;
    	}
    	else if(g_ind0>0 && g_ind0<=gnxy){ 
    		double pt_adf=pm2d->geo_bgrid_adf(x,y);
    		double pt_bgrid_geps=get_bgrid_geps(x,y);
    		if(abs(pt_adf)<=pt_bgrid_geps){ //bdry pt
    			pm2d->pt_ctgr[i]=1;
    		}
    		else if(pt_adf<-pt_bgrid_geps){
				pm2d->pt_ctgr[i]=-1; //-1: inner pt
    			inner_pt_ct_temp+=1;
			}
			else{ //pt_adf>pt_bgrid_geps: outside of geometry with new geps
				if(i<Nfixed){
					printf("ERROR: Add new pt procedure: fixed pt (%g %g) in geo bdry grid but not in/on geometry.\n",x,y);
			    	throw std::exception();
				}

				double new_x, new_y;
				if(pt_projection(new_x,new_y,x,y,-1,-1)==true){ //(new_x,new_y) on geo bdry/inner grid now
					pm2d->xy_id[i2]=new_x;
					pm2d->xy_id[i2+1]=new_y;

					//decide if new pt is bdry pt or inner pt
			    	int g_ind=geo_grid(new_x,new_y);
			    	if(g_ind>0 && g_ind<=gnxy){ //bdry grid
			    		double new_xy_adf=pm2d->geo_bgrid_adf(new_x,new_y);
			    		double new_xy_bgrid_geps=get_bgrid_geps(new_x,new_y);
			    		if(abs(new_xy_adf)<=new_xy_bgrid_geps){ //bdry pt
			    			pm2d->pt_ctgr[i]=1;
			    		}
			    		else if(new_xy_adf<-new_xy_bgrid_geps){
			    			pm2d->pt_ctgr[i]=-1; //-1: inner pt
			    			inner_pt_ct_temp+=1;
			    		}
			    		else{
			    			printf("ERROR: add_pt: projected bdry pt now (%g %g) in geo bdry grid but not in/on geometry.\n",new_x,new_y);
			    			throw std::exception();
			    		}
			    	}
			    	else if(g_ind<0){ //inner grid
			    		pm2d->pt_ctgr[i]=-1; //-1: inner pt
			    		inner_pt_ct_temp+=1;
			    	}
			    	else{
			    		printf("ERROR: add_pt: projected bdry pt now (%g %g) not in geometry inner/bdry grid.\n",new_x,new_y);
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
					pm2d->xy_id[i2]=random_point_x;
					pm2d->xy_id[i2+1]=random_point_y;
					pm2d->pt_ctgr[i]=-1; //tag inner pt
                	inner_pt_ct_temp+=1;

				}
			}
		}
		else{
			if(i<Nfixed){
				printf("ERROR: Add new point procedure: Fixed pt (%g %g) not in inner nor geo bdry grid.\n",x,y);
    			throw std::exception();
			}
			else{
				printf("ERROR: Add new point procedure: prev bdry pt (%g %g) not in inner nor geo bdry grid.\n",x,y);
				throw std::exception();
			}
    		
    	}
		
	}
	pm2d->inner_pt_ct=inner_pt_ct_temp;

printf("finished added points \n");
//---------------------------------------------------------------------

}

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
void mesh_alg_2d::do_sphere_tracing(double &final_x,double &final_y,double px_old,double py_old,double px_new,double py_new){
	double dis=pm2d->sdf_adf(px_old,py_old);
    double adis=fabs(dis);
    double v_len=d_points(px_old,py_old,px_new,py_new);
    double v_unit_x=(px_new-px_old)/v_len;
    double v_unit_y=(py_new-py_old)/v_len;
    double tempx=px_old;
    double tempy=py_old;
    double geps_temp=get_bgrid_geps(px_new,py_new);
    double old_tempx=tempx; 
    double old_tempy=tempy;
    bool binary=false; bool binary_full=false;
    while(dis<-geps_temp){ //point from interior going to edge 
        old_tempx=tempx; 
        old_tempy=tempy;
        tempx=tempx+adis*v_unit_x;
        tempy=tempy+adis*v_unit_y;
        //if (tempx,tempy) now in bdry grid, update geps_temp
        if(geo_grid(tempx,tempy)>0 && geo_grid(tempx,tempy)<=gnxy){
        	geps_temp=get_bgrid_geps(tempx,tempy);
        }
        dis=pm2d->sdf_adf(tempx,tempy);
        adis=fabs(dis);
        if(dis>geps_temp){ //for the case when SDF is not so accurate
            //as this shouldn't happen is SDF is accuate, now we 
            //switch to using binary search for the root
            binary=true; 
        }
    }
    double dis_xy_tempxy=d_points(px_old,py_old,tempx,tempy); //should be smaller than v_len
    if(dis_xy_tempxy>v_len){
        binary=true;
        binary_full=true;
    }
    if(binary==true){
        if(binary_full==true){
            //-: x,y
            //+: cx, cy
            old_tempx=px_old;
            old_tempy=py_old;
            tempx=px_new;
            tempy=py_new;
            geps_temp=get_bgrid_geps(px_new,py_new);
        }

        //-: old_tempx, old_tempy
        //+: tempx, tempy
        double mid_tempx=0.5*(old_tempx+tempx);
        double mid_tempy=0.5*(old_tempy+tempy);
        //if (mid_tempx,mid_tempy) now in bdry grid, update geps_temp
        if(geo_grid(mid_tempx,mid_tempy)>0 && geo_grid(mid_tempx,mid_tempy)<=gnxy){
        	geps_temp=get_bgrid_geps(mid_tempx,mid_tempy);
        }
        dis=pm2d->sdf_adf(mid_tempx,mid_tempy);
        adis=fabs(dis);
        while(adis>geps_temp){
            if(dis<-geps_temp){
                old_tempx=mid_tempx;
                old_tempy=mid_tempy;
                mid_tempx=0.5*(old_tempx+tempx);
                mid_tempy=0.5*(old_tempy+tempy);
                //if (mid_tempx,mid_tempy) now in bdry grid, update geps_temp
		        if(geo_grid(mid_tempx,mid_tempy)>0 && geo_grid(mid_tempx,mid_tempy)<=gnxy){
		        	geps_temp=get_bgrid_geps(mid_tempx,mid_tempy);
		        }
        		dis=pm2d->sdf_adf(mid_tempx,mid_tempy);
        		adis=fabs(dis);
            }
            else{
                tempx=mid_tempx;
                tempy=mid_tempy;
                mid_tempx=0.5*(old_tempx+tempx);
                mid_tempy=0.5*(old_tempy+tempy);
                //if (mid_tempx,mid_tempy) now in bdry grid, update geps_temp
		        if(geo_grid(mid_tempx,mid_tempy)>0 && geo_grid(mid_tempx,mid_tempy)<=gnxy){
		        	geps_temp=get_bgrid_geps(mid_tempx,mid_tempy);
		        }
        		dis=pm2d->sdf_adf(mid_tempx,mid_tempy);
        		adis=fabs(dis);
            }
        }
        tempx=mid_tempx;
        tempy=mid_tempy;
    }
    final_x=tempx;
    final_y=tempy;
}


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
void mesh_alg_2d::pt_mvmt_treatment_projection(int pid,double px_old,double py_old,double px_new,double py_new, 
	double &final_x, double &final_y){
	if(pt_ctgr(pid)==0){ 
		printf("ERROR: pt (px_old,py_old) is neither inner nor bdry point.\n");
    	throw std::exception();
	}

	final_x=px_old;
	final_y=py_old;
	double mvmt_dis=d_points(px_old,py_old,px_new,py_new);
	int i_old=(px_old-ax)*inv_gdx;
    int j_old=(py_old-ay)*inv_gdy;
    int ij_old=j_old*gnx+i_old;
	
	
	//A. check pt mvmt_dis if too large: 
	//if pt mvmt_dis>pt_mvmt_dis_thres[ij_old], move pt in p_new direction for pt_mvmt_dis_thres[ij_old], and get new p_new.
    if(mvmt_dis>pt_mvmt_dis_thres[ij_old]){

    	double new_dis_fac=pt_mvmt_dis_thres[ij_old]/mvmt_dis;
    	double new_dis_dx=(px_new-px_old)*new_dis_fac;
    	double new_dis_dy=(py_new-py_old)*new_dis_fac;

    	px_new=px_old+new_dis_dx;
    	py_new=py_old+new_dis_dy;
    	mvmt_dis=d_points(px_old,py_old,px_new,py_new);
    }

    //B. check if p_new lies in inner/bdry/outer grid
    //if p_new in outer grid, further shorten its mvmt distance, until it lies in inner/bdry grid.
    int geo_grid_ij_p_new=geo_grid(px_new,py_new);
    if(geo_grid_ij_p_new>gnxy){ //p_new at outer grid, move pt for d=0.9d in p_new direction instead. Repeat until p_new at bdry/inner grid.
    	while(geo_grid_ij_p_new>gnxy){
    		double new_dis_dx=(px_new-px_old)*0.9;
    		double new_dis_dy=(py_new-py_old)*0.9;

    		px_new=px_old+new_dis_dx;
    		py_new=py_old+new_dis_dy;
    		geo_grid_ij_p_new=geo_grid(px_new,py_new);

    	}
    }

    //C. now, p_new is in inner/bdry grid, pt treatment & projection
    //New point categories
    if(geo_grid_ij_p_new<0){ //p_new at inner grid: move p_old to p_new 
    	final_x=px_new;
    	final_y=py_new;
    	if(pt_ctgr(pid)==-1){ //pt (px_old,py_old) was inner point
    		//pt p_new is inner pt too. No need to update pt_ctgr or inner_pt_ct
		}
		else{ //pt_ctgr(pid)==1, pt (px_old,py_old) was boundary point. Need to update pt_ctgt and inner_pt_ct
			pm2d->pt_ctgr[pid]=-1;
			#pragma omp atomic
			pm2d->inner_pt_ct++;
		}
    }
    else if(geo_grid_ij_p_new>0 && geo_grid_ij_p_new<=gnxy){ //p_new at bdry grid
    	//a. test if p_new SDF outside geometry
    	//if p_new inside/on bdry of geometry, move p_old to p_new. 
    	double sdf_val_p_new=pm2d->sdf_adf(px_new,py_new);
    	double p_new_bgrid_geps=get_bgrid_geps(px_new,py_new);
    	if(sdf_val_p_new<=p_new_bgrid_geps){  //inner/bdry pt
    		final_x=px_new;
    		final_y=py_new;
    		if(fabs(sdf_val_p_new)<=p_new_bgrid_geps){ //p_new bdry pt
    			if(pt_ctgr(pid)==-1){ //pt (px_old,py_old) was inner point
    				pm2d->pt_ctgr[pid]=1;
    				#pragma omp atomic
					pm2d->inner_pt_ct--;
				}
				else{ //pt_ctgr(pid)==1, pt (px_old,py_old) was boundary point
				}
    		}
    		else if(sdf_val_p_new<-p_new_bgrid_geps){ //p_new inner pt
    			if(pt_ctgr(pid)==-1){ //pt (px_old,py_old) was inner point
				}
				else{ //pt_ctgr(pid)==1, pt (px_old,py_old) was boundary point
					pm2d->pt_ctgr[pid]=-1;
					#pragma omp atomic
					pm2d->inner_pt_ct++;
				}
    		}
    		else{
    			printf("ERROR: pt (px_new,py_new) is neither inner nor bdry point. Should not happen.\n");
    			throw std::exception();
    		}
    	}
    	else{ //p_new outside of geometry. Need to move it back to geo bdry.
    		if(pt_ctgr(pid)==-1){ //pt (px_old,py_old) was inner point
    			//sphere tracing to find bdry pt, update p_new.
    			do_sphere_tracing(final_x,final_y,px_old,py_old,px_new,py_new);
    			//Tag pt bdry pt, inner_pt_ct--
    			pm2d->pt_ctgr[pid]=1;
    			#pragma omp atomic
				pm2d->inner_pt_ct--;
			}
			else{ //pt_ctgr(pid)==1, pt (px_old,py_old) was boundary point
				 double geps_prime=get_bgrid_geps(px_new,py_new);
				 double deps_prime=get_bgrid_deps(px_new,py_new);
				 bool store=true;

				//projection of p_new. 
				if(pt_projection(final_x,final_y,px_new,py_new,deps_prime,geps_prime)){ //projection successful. find pt on bdry/inner grid
					
					//test1: dis(p_new,p_proj)<=dis(p_new,p_old) && dis(p_proj,p_old)<=dis(p_new,p_old), then accept. Otherwise, pt stays at p_old
					double dis_p_new_p_old=d_points(px_old,py_old,px_new,py_new);
					double dis_p_new_p_proj=d_points(px_new,py_new,final_x,final_y);
					double dis_p_old_p_proj=d_points(px_old,py_old,final_x,final_y);
					if(dis_p_new_p_proj>dis_p_new_p_old || dis_p_old_p_proj>dis_p_new_p_old){
						store=false;
					}

					if(store==true){
						//test2:
						//Check if pt lie on bdry grid:
						//and has sdf<geps_i
						double final_xy_sdf=pm2d->sdf_adf(final_x,final_y);
						int final_xy_g_ind=geo_grid(final_x,final_y);
						if(final_xy_g_ind>0 && final_xy_g_ind<=gnxy){ //bdry grid
							double final_xy_bgrid_geps=get_bgrid_geps(final_x,final_y);
							if(fabs(final_xy_sdf)<=final_xy_bgrid_geps){ //bdry pt
								store=true; //as bdry pt
							}
							else if(final_xy_sdf<-final_xy_bgrid_geps){ //inner pt
								store=true; //as inner pt
								pm2d->pt_ctgr[pid]=-1;
								#pragma omp atomic
								pm2d->inner_pt_ct++;
							}
							else{ //outer pt in bdry grid
								store=false;
							}
						}
						else if(final_xy_g_ind<0){ //inner grid
							store=true; //as inner pt
							pm2d->pt_ctgr[pid]=-1;
							#pragma omp atomic
							pm2d->inner_pt_ct++;

						}
						else{ //outer grid: should not happen, since projection() return false if pt in outer grid
							printf("ERROR: pt after projection return true lie in outer grid. Should not happen.\n");
    						throw std::exception();
						}
					}
					

					
				}
				else{ //projection diverges. Keep pt at p_old
					store=false;
				}

				if(store==false){
					final_x=px_old;
					final_y=py_old;
				}

				//bdry pt->bdry pt, no need to update pt ctgr nor inner_pt_ct
			}
    	}
    }
    else{
    	printf("ERROR: pt (px_new,py_new) is neither inner nor bdry point. Should not happen.\n");
    	throw std::exception();
    }
}


/**
 * Check if the movement of an inner point exceeds the stop_Continue_mvmt_thres threshold.
 * This function compares the movement distance of an inner point with the threshold.
 *
 * @param[in] pi The index of the inner point.
 * @param[in] px_old The x-coordinate of the old point.
 * @param[in] py_old The y-coordinate of the old point.
 * @param[in] px_new The x-coordinate of the new point.
 * @param[in] py_new The y-coordinate of the new point.
 * @return True if the inner point movement exceeds the threshold, False otherwise.
 */
bool mesh_alg_2d::inner_pt_movement_over_stop_Continue_mvmt_thres(int pi, double px_old, double py_old, double px_new, double py_new){
	//test if inner point pi movement larger than stop_Continue_mvmt_thres
    if(pt_ctgr(pi)==-1){ //-1: inner pt
        double pt_movement=d_points(px_old,py_old,px_new,py_new);
        int i_old=(px_old-ax)*inv_gdx;
        int j_old=(py_old-ay)*inv_gdy;
        int ij_old=j_old*gnx+i_old;
        if(pt_movement>stop_Continue_mvmt_thres[ij_old]){
            return true;
        }
    }
    
    return false;
}

/**
 * Check if the movement of an inner point exceeds the stop_Continue_mvmt_thres threshold.
 * This function compares the movement distance of an inner point with the threshold.
 *
 * @param[in] pi The index of the inner point.
 * @param[in] pt_movement The movement distance of the inner point.
 * @param[in] stop_Continue_mvmt_thres_ij The threshold for the inner point movement.
 * @return True if the inner point movement exceeds the threshold, False otherwise.
 */
bool mesh_alg_2d::inner_pt_movement_over_stop_Continue_mvmt_thres(int pi,double pt_movement,double stop_Continue_mvmt_thres_ij){
	//test if inner point pi movement larger than stop_Continue_mvmt_thres
    if(pt_ctgr(pi)==-1){ //-1: inner pt
        if(pt_movement>stop_Continue_mvmt_thres_ij){
            return true;
        }
    }
    
    return false;
}

/**
 * Determine whether to continue or stop meshing based on the current state of the algorithm.
 * If the number of current points is equal to the total number of points and there are no inner points that exceed the stop_Continue_mvmt_thres threshold, then meshing will be stopped.
 */
void mesh_alg_2d::determine_Continue_pt_movement(){
	if(Ncurrent==Ntotal){
		if(inner_pt_over_stop_Continue_mvmt_thres_ct==0){
			Continue=false;
		}
	}
}

/**
 * @brief Determine whether to add new points and retriangulate points in the next iteration based on triangle quality.
 *
 * This function evaluates the relative change in median aspect ratio and median element quality to decide
 * whether new points should be added and points should be retriangulated in the next iteration. If the relative changes in 
 * mean aspect ratio (`ar_alpha_mean_rel_change`) and mean edge ratio (`er_alpha_mean_rel_change`) 
 * are both below a certain threshold (`determine_addPt_relative_change_mean_ar` and 
 * `determine_addPt_relative_change_mean_er` respectively), the function sets the `addPt` and `reTria` flags
 * to true, indicating that new points should be added and points should be retriangulated in the next iteration.
 *
 * Additionally, the function may call the `hybrid_determine_method_switch` function to perform additional
 * processing based on the calculated relative change in mean aspect ratio.
 *
 * After processing, the function updates the previous values of the mean aspect ratio and mean edge ratio 
 * for future comparisons.
 */
void mesh_alg_2d::determine_addPt_and_reTria_tria_quality(){
	//Relative criteria: Triangle quality
	//if relative change in median aspect ratio < etermine_addPt_relative_change_mean_ar=0.005
	//add pt	   

	double er_alpha_mean_rel_change=fabs(alpha_mean_tria_er-previous_alpha_mean_tria_er)/previous_alpha_mean_tria_er;
    double ar_alpha_mean_rel_change=fabs(alpha_mean_tria_ar-previous_alpha_mean_tria_ar)/previous_alpha_mean_tria_ar;

    if(ar_alpha_mean_rel_change<determine_addPt_relative_change_mean_ar
    	&&
    	er_alpha_mean_rel_change<determine_addPt_relative_change_mean_er)
    {
    	addPt=true;
		reTria=true;
    }

    hybrid_determine_method_switch(ar_alpha_mean_rel_change);

    previous_alpha_mean_tria_ar=alpha_mean_tria_ar; 
    previous_alpha_mean_tria_er=alpha_mean_tria_er;
}


/**
 * @brief Determine whether to add new points and retriangulate points in the next iteration based on movement and iteration count.
 *
 * This function evaluates the absolute criteria of point movement and iteration count to decide whether new points
 * should be added and points should be retriangulated in the next iteration.
 *
 * If the iteration count (`inbtw_tria_iter_ct`) exceeds a certain threshold (`inbtw_tria_iter_ct_thres`), or if the
 * proportion of inner points with movement over a specified threshold (`proportion_inner_pt_mvmt_over_thres`) falls
 * below a criteria value (`proportion_inner_pt_mvmt_over_thres_criteria`), the function sets the `addPt` and `reTria`
 * flags to true, indicating that new points should be added and points should be retriangulated in the next iteration.
 */
void mesh_alg_2d::determine_addPt_and_reTria_movement(){
	//Absolute criteria: Pt movement, iterations ct
	//if a large proportion of inner points have small movement 
	//or if inbtw_tria_iter_ct too large
	//then add new points in
	double proportion_inner_pt_mvmt_over_thres=1.0*inner_pt_over_stop_Continue_mvmt_thres_ct/pm2d->inner_pt_ct;
	if(inbtw_tria_iter_ct>inbtw_tria_iter_ct_thres||
		proportion_inner_pt_mvmt_over_thres<proportion_inner_pt_mvmt_over_thres_criteria)
	{
		addPt=true;
		reTria=true;
	}
	
}

/**
 * @brief Update the stop_Continue_mvmt_thres threshold for a specific grid point.
 *
 * This function updates the stop_Continue_mvmt_thres threshold for a given grid point index.
 * If the corresponding geometric grid value is less than or equal to gnxy, which represents
 * an inner or boundary grid, the stop_Continue_mvmt_thres value for that index is updated
 * based on the stop_Continue_fac and chrtrt_len_h factors.
 *
 * @param ij The index of the grid point for which the threshold needs to be updated.
 *
 */
void mesh_alg_2d::update_stop_Continue_mvmt_thres(int ij){
    if(geo_grid(ij)<=gnxy){ //inner or bdry grid 
        stop_Continue_mvmt_thres[ij]=stop_Continue_fac*chrtrt_len_h(ij);
    }
}

/**
 * @brief Update the pt_mvmt_dis_thres threshold for a specific grid point.
 *
 * This function updates the pt_mvmt_dis_thres threshold for a given grid point index.
 * If the corresponding geometric grid value is less than or equal to gnxy, which represents
 * an inner or boundary grid, the pt_mvmt_dis_thres value for that index is updated
 * based on the pt_mvmt_dis_thres_fac and chrtrt_len_h factors.
 *
 * @param ij The index of the grid point for which the threshold needs to be updated.
 *
 */
void mesh_alg_2d::update_pt_mvmt_dis_thres(int ij){
    if(geo_grid(ij)<=gnxy){ //inner or bdry grid 
        pt_mvmt_dis_thres[ij]=pt_mvmt_dis_thres_fac*chrtrt_len_h(ij);
    }
}


/**
 * @brief Initialize the meshing algorithm.
 *
 * This function initializes the meshing algorithm by performing the following steps:
 * - Initialize the stop_Continue_mvmt_thres and pt_mvmt_dis_thres arrays based on the geometry sizing grid.
 * - Update the stop_Continue_mvmt_thres and pt_mvmt_dis_thres values for each grid point.
 * - Add walls to Voronoi cells.
 * - Generate initial points for meshing.
 */
void mesh_alg_2d::meshing_init(){
t_add_pt_centroid=0.0;
double t0=omp_get_wtime();
	//stop_Continue_mvmt_thres threshold init: adaptive to geometry sizing grid
	//-1:dummy; if all inner points have movement <thres, Continue=false
    delete [] stop_Continue_mvmt_thres;
    stop_Continue_mvmt_thres=new double[gnxy];
    delete [] pt_mvmt_dis_thres;
    pt_mvmt_dis_thres=new double[gnxy];
	
    #pragma omp parallel for num_threads(num_t)
    for(int ij=0;ij<gnxy;ij++){
    	mesh_alg_init(ij);
        stop_Continue_mvmt_thres[ij]=-1;
        pt_mvmt_dis_thres[ij]=-1;
        update_stop_Continue_mvmt_thres(ij);
        update_pt_mvmt_dis_thres(ij);
    }

    pm2d->con->add_wall(wis);
	printf("Initial %d points\n", Ncurrent);
double t1=omp_get_wtime();
t_meshing_init=t1-t0;
t_algorithm_all_non_voro=t_meshing_init;
}


/**
 * @brief Perform meshing iterations until the stopping condition is met.
 *
 * This function performs meshing iterations until the stopping condition is met. The iterations involve the following steps:
 * - Enter a while loop that continues as long as the Continue flag is true.
 * - Increment the all_iter_ct counter and reset the inner_pt_over_stop_Continue_mvmt_thres_ct counter.
 * - Determine whether retriangulation is necessary by calling the determine_reTria function.
 * - If retriangulation is required, perform the following actions:
 *     - Increment the tria_iter_ct counter and reset the inbtw_tria_iter_ct counter.
 *     - If new points need to be added, call the add_new_pts function.
 *     - Reset the addPt flag and print the current number of points.
 *     - Perform Voronoi computation and extraction of triangle information by calling the voro_compute_retria_extract_info function.
 *     - Reset the reTria flag.
 *     - If the current number of points is equal to the total number of points, call the determine_Continue_tria_quality function to determine whether meshing should continue based on triangle quality.
 *     - If the current number of points is less than the total number of points, call the determine_addPt_and_reTria_tria_quality function to determine whether new points should be added and whether retriangulation should occur based on triangle quality.
 * - If the Continue flag is true and addPt is false, perform the following actions:
 *     - Increment the inbtw_tria_iter_ct counter.
 *     - Update the positions of the points and determine whether retriangulation is necessary by calling the update_pt_position_and_determine_reTria function.
 *     - Test the stopping criteria by calling the determine_Continue_pt_movement function if the current number of points is equal to the total number of points, or the determine_addPt_and_reTria_movement function if the current number of points is less than the total number of points.
 * - Exit the while loop when the Continue flag is false.
 */
void mesh_alg_2d::meshing_iter(){
double t0=omp_get_wtime();
t_add_pt=0.0;
t_voro_computation=0.0;
t_update_pt_position=0.0;

t_dm_getBarinfo=0.0;
t_dm_applyBarForce=0.0;
t_dm_updatePtPos=0.0;


	while(Continue){
		all_iter_ct++;
		inner_pt_over_stop_Continue_mvmt_thres_ct=0;

		//determine retriangulation or not, update bool reTria
		determine_reTria();
		//Retriangulation
		if(reTria==true){

			tria_iter_ct++;
			inbtw_tria_iter_ct=0;
			if(addPt==true){
double t1=omp_get_wtime();
				add_new_pts();
t_add_pt+=omp_get_wtime()-t1;
				addPt=false;
				printf("added pts. now: %d points; all_iter_ct %d tria_iter_ct %d \n", Ncurrent, all_iter_ct,tria_iter_ct);
			}
			//Voronoi computation
				//triangle extraction, triangle ct
				//tria quality calculation
				//*if meshing algorithm is DM: need to store triangle bars: barid
				//*if meshing algorithm is CVD: need to calcuate pt new position inside Voro computation, and store in xy_id
double t2=omp_get_wtime();
			voro_compute_retria_extract_info();
t_voro_computation+=omp_get_wtime()-t2;
//print_tria_quality_stat_overall();
			reTria=false;

			//test triangle quality stopping criteria
			if(Ncurrent==Ntotal){
				//if all triangles reach triangle quality threshold: Continue=false
				//if triangle quality not improving much anymore: Continue=false
				//In hybrid, this also determines to switch from DM to CVD or not for final refinement of mesh quality
				determine_Continue_tria_quality();
			}
			else{
				determine_addPt_and_reTria_tria_quality();
			}
		}

		//triangles have not reach the quality stopping criteria, continue meshing iterations
		if(Continue==true && addPt==false){
			inbtw_tria_iter_ct++;
			//use the chosen meshing algorithm (dm/cvd) to obtain the new point position, xy_id_new
			//point treatment: 
			//update xy_id to new position, based on rules & projection() if outside of geometry
			//determine reTria: true for CVD, false for dm
double t3=omp_get_wtime();
			update_pt_position_and_determine_reTria();
t_update_pt_position+=omp_get_wtime()-t3;

//printf("finished update points, new check \n");
//pm2d->check_pt_ctgr();

			//test stopping criteria: if all points movement small in this meshing iteration, stop
			if(Ncurrent==Ntotal){
				determine_Continue_pt_movement();
			}
			else{ //Ncurrent<Ntotal
				//test add point criteria: 
					//if a large portion of pts have small movement, 
					//or in-between-triangulation-iterations > threshold, add points
				determine_addPt_and_reTria_movement();
			}
		}	
	//exiting while(Continue) loop	
	}
double t4=omp_get_wtime();
t_algorithm_remain=t4-t0
					-t_add_pt
					-t_voro_computation
					-t_update_pt_position;
t_algorithm_all_non_voro+=(t4-t0-t_voro_computation);
}


/**
 * @brief Perform the final meshing iteration.
 *
 * This function performs the final meshing iteration. It involves the following steps:
 * - Increment the tria_iter_ct counter.
 * - Call the determine_printOutputs function to determine if the outputs need to be printed.
 * - Perform the last round of retriangulation and store triangle information if necessary.
 * - Call the voro_compute_and_store_info function to perform Voronoi computation and store information.
 * - Output the results if necessary.
 */
void mesh_alg_2d::meshing_final(){
double t0=omp_get_wtime();
	tria_iter_ct++;
	determine_printOutputs();

	//last round of retriangulation, and store triangle info for print out later
	bool store_barid=false;
	bool store_triangle_quality=false;
	bool store_tria_vertex=false;
    bool store_tria_centroid=false;
    bool store_ccrd_h_ratio=false;
    bool store_max_tria_ar=false;
    bool cvd_pt_update=false;
    bool store_xy_id=false;
    double *dummy_array=new double[1];
    if(printOutputs==true){
		//store_barid=true;
		store_triangle_quality=true;
		store_tria_vertex=true;
	}
double t1=omp_get_wtime();
	voro_compute_and_store_info(store_barid,store_triangle_quality,
		store_tria_vertex,
		store_tria_centroid,
		store_ccrd_h_ratio, store_max_tria_ar,
		cvd_pt_update,
		store_xy_id,dummy_array);
double t2=omp_get_wtime();
t_voro_computation+=t2-t1;

	delete [] dummy_array;

	//output
	if(printOutputs==true){
		construct_HE();
		do_print_final_outputs();
	}
double t3=omp_get_wtime();
//print_tria_quality_stat_overall();
t_algorithm_all_non_voro+=(t3-t0-(t2-t1));
}

/**
 * @brief Construct the Half-edge data structure for the triangulation.
 * 		  The boundary HE are those with NULL face pointer.
 * 		  Requires tria_vertex to have been stored in the voro-tria calculation.
 */
void mesh_alg_2d::construct_HE(){

printf("constructing HE\n");
	//clear barid storage to avoid redundancy in storage space
	delete [] barid;
	barid=new int[2*1];
	barid_barct=1;
	
	//Initialize variables
	Edges.clear();
	Bdry_Edges_Start.clear();
	delete [] Faces;
	delete [] Vertices;
	Faces=new HE_face[tria_ct];
	Vertices=new HE_vert[Ncurrent];
	

	//Make sure triangle vertices are CCW oriented
	if(tria_order_ccw==false){sort_tria_vertex_ids_ccw();}
printf("A\n");
	//Construct HE data structure
	//Loop through the triangles
	for(int triai=0;triai<tria_ct;triai++){
		int Aid=tria_vertex[3*triai];
		int Bid=tria_vertex[3*triai+1];
		int Cid=tria_vertex[3*triai+2];

		std::vector<std::pair<unsigned int, unsigned int>> tedges;
		tedges.push_back(std::make_pair(Aid, Bid));
		tedges.push_back(std::make_pair(Bid, Cid));
		tedges.push_back(std::make_pair(Cid, Aid));

		//Construct all the HE (u,v) and (v,u) related to the triangle
		for(int i=0;i<3;i++){
			std::pair<unsigned int, unsigned int> tedges_i=tedges[i];
			std::pair<unsigned int, unsigned int> tedges_i_pair=std::make_pair(tedges[i].second, tedges[i].first);
			if(Edges.find(KKey(tedges_i))==Edges.end()){
				Edges[KKey(tedges_i)]=new HE_edge();
				Edges[KKey(tedges_i)]->eu=tedges_i.first;
				Edges[KKey(tedges_i)]->ev=tedges_i.second;
			}
			if(Edges.find(KKey(tedges_i_pair))==Edges.end()){
				Edges[KKey(tedges_i_pair)]=new HE_edge();
				Edges[KKey(tedges_i_pair)]->eu=tedges_i_pair.first;
				Edges[KKey(tedges_i_pair)]->ev=tedges_i_pair.second;
			}
			if(Edges[KKey(tedges_i)]->pair == NULL){
				Edges[KKey(tedges_i)]->pair = Edges[KKey(tedges_i_pair)];
			}
			if(Edges[KKey(tedges_i_pair)]->pair == NULL){
				Edges[KKey(tedges_i_pair)]->pair = Edges[KKey(tedges_i)];
			}
			Edges[KKey(tedges_i)]->face=&Faces[triai];

			if(Edges[KKey(tedges_i)]->vert == NULL){
				Edges[KKey(tedges_i)]->vert = &Vertices[tedges_i.second];
			}
			if(Edges[KKey(tedges_i_pair)]->vert == NULL){
				Edges[KKey(tedges_i_pair)]->vert = &Vertices[tedges_i.first];
			}
			//HE vertice u
			if(Vertices[tedges_i.first].edge==NULL){
				Vertices[tedges_i.first].edge=Edges[KKey(tedges_i)];
				Vertices[tedges_i.first].HE_vert_id=tedges_i.first;
			}

			//HE triangle Face triai
			if(i==0){
				Faces[triai].edge=Edges[KKey(tedges_i)];
				Faces[triai].HE_face_id=triai;
			}
		}
		//Fill the HE data structure for the edges, vertices, and faces adjacent to the triangle
		for(int i=0;i<3;i++){
			std::pair<unsigned int, unsigned int> tedges_i=tedges[i];
			std::pair<unsigned int, unsigned int> tedges_i_pair=std::make_pair(tedges_i.second, tedges_i.first);

			//HE for (u,v)
			if(i==2){Edges[KKey(tedges_i)]->next=Edges[KKey(tedges[0])];} else{Edges[KKey(tedges_i)]->next=Edges[KKey(tedges[i+1])];}
			if(i==0){Edges[KKey(tedges_i)]->prev=Edges[KKey(tedges[2])];} else{Edges[KKey(tedges_i)]->prev=Edges[KKey(tedges[i-1])];}
		}
	}

printf("B\n");

	//Fix single triangle holes:
	//First, find vertices v with >1 (with 2) outgoing bdry edges (i.e. connecting to 4 bdry HE)
	//Obtain boundary edges connectivity: bdry HE prev and next
	//Loop through all edges and obtain the bdry HE's
	std::vector<std::vector<int>> eu_bdry_connect; //Vector: [index] ev1, (ev2)
	std::unordered_map<int,int> eu_bdry_ind; //Mapping: key: eu ---- Value: index
	std::vector<int> eu_problematic; //problematic eu's: Vector: eu1,eu2,eu3...
	int bdry_ind=0;
	for(std::unordered_map<size_t, HE_edge*>::iterator itr_e=Edges.begin();itr_e!=Edges.end();itr_e++){
		if(Edges[itr_e->first]->face==NULL){
			int eu=itr_e->second->eu;
			int ev=itr_e->second->ev;
			if(eu_bdry_ind.find(eu)==eu_bdry_ind.end()){
				eu_bdry_ind[eu]=bdry_ind++;
				std::vector<int> bdry_ind_temp; bdry_ind_temp.push_back(ev);
				eu_bdry_connect.push_back(bdry_ind_temp);
			}
			else{
				eu_bdry_connect[eu_bdry_ind[eu]].push_back(ev);
				eu_problematic.push_back(eu);
			}
		}
	}


	//loop through the problematic vertices
	for(int i=0;i<(signed int)eu_problematic.size();i++){
printf("HERE O\n");
		//Loop through the connecting incoming edges and the corresponding triangle
		bool have_deleted=false; //only delete one of the connecting triangles
		int eu=eu_problematic[i];
		if((signed int)eu_bdry_connect[eu_bdry_ind[eu]].size()>2){
			printf("ERROR: more than 4 bdry half edge connecting to a single bdry vertex %d (%g %g)!\n %d \n nei: %d (%g %g), %d(%g %g)\n", 
				eu, current_x(eu),current_y(eu),
				eu_bdry_connect[eu_bdry_ind[eu]].size(), 
				eu_bdry_connect[eu_bdry_ind[eu]][0], current_x(eu_bdry_connect[eu_bdry_ind[eu]][0]),current_y(eu_bdry_connect[eu_bdry_ind[eu]][0]),
				eu_bdry_connect[eu_bdry_ind[eu]][1], current_x(eu_bdry_connect[eu_bdry_ind[eu]][1]),current_y(eu_bdry_connect[eu_bdry_ind[eu]][1])

				);
			throw std::exception();
		}
		//delete a single triangle who two edges are bdry edges and connecting to the vertex
		for(int ej=0;ej<(signed int)eu_bdry_connect[eu_bdry_ind[eu]].size();ej++){
			if(have_deleted==false){
printf("HERE A\n");
				int ev=eu_bdry_connect[eu_bdry_ind[eu]][ej];
				//check if the pair edge's Next is the pair edge of another (a)bdry edge and (b)connecting to vertec eu
				if(Edges[KKey(eu,ev)]->pair->next->pair->face==NULL){
					//Delete the triangle associated with the Pair HE of eu,ev.
					//And move the last triangle to the slot
					int problem_tria_id=Edges[KKey(eu,ev)]->pair->face->HE_face_id;
					int problem_tria_id3=3*problem_tria_id;
					int last_tria_id=tria_ct-1;
					int last_tria_id3=3*last_tria_id;
					tria_vertex[problem_tria_id3]=tria_vertex[last_tria_id3];
					tria_vertex[problem_tria_id3+1]=tria_vertex[last_tria_id3+1];
					tria_vertex[problem_tria_id3+2]=tria_vertex[last_tria_id3+2];

					//Update triangle ct
					tria_ct--;

					//Update the HE data strcuture: Faces, Edges, Vertices
					Faces[problem_tria_id].HE_face_id=last_tria_id;
					Faces[problem_tria_id].edge=Faces[last_tria_id].edge;

					//Delete HE (eu,ev) (ev,eu) (eu,ev2) (ev2,eu)
					int ev2=Edges[KKey(eu,ev)]->pair->next->ev;
					//Make new bdry HE
					Edges[KKey(eu,ev)]->pair->next->next->face=NULL;
					Edges[KKey(eu,ev)]->pair->next->next->prev=NULL;
					Edges[KKey(eu,ev)]->pair->next->next->next=NULL;
					//Erase deleted triangle edges
					Edges.erase(KKey(eu,ev));
					Edges.erase(KKey(ev,eu));
					Edges.erase(KKey(eu,ev2));
					Edges.erase(KKey(ev2,eu));

					int ejj=ej+1; if(ejj>=(signed int)eu_bdry_connect[eu_bdry_ind[eu]].size()){ejj=0;}
					Vertices[eu].edge=Edges[KKey(eu,eu_bdry_connect[eu_bdry_ind[eu]][ejj])];
					Vertices[ev].edge=Edges[KKey(ev,ev2)];
					Vertices[ev2].edge=Edges[KKey(ev2,ev)];

					//Change the flag to signal a connecting triangle has been deleted.
					have_deleted=true;


				}
			}
		}
		//No single triangle found; 
		//Simply delete all triangles bounded by two bdry edges on one side connecting to the vertex
		//PS. choose the side with less number of triangles
		if(have_deleted==false){

printf("HERE B\n");
			int ej_chose=0;
			int e_tria_ct_min=100;
			for(int ejjj=0;ejjj<(signed int)eu_bdry_connect[eu_bdry_ind[eu]].size();ejjj++){
				int evvv=eu_bdry_connect[eu_bdry_ind[eu]][ejjj];
				//count the number of triangles on this side:
				int e_tria_ct=1;
				HE_edge *checkE=Edges[KKey(eu,evvv)]->pair->next->pair;
				while(checkE->face!=NULL){
					checkE=checkE->next->pair;
					e_tria_ct++; printf("!!!\n");
				}
				if(e_tria_ct<e_tria_ct_min){e_tria_ct_min=e_tria_ct;ej_chose=ejjj;}
			}
			printf("triangles to delete near bdry: %d \n",e_tria_ct_min);

			int evvv=eu_bdry_connect[eu_bdry_ind[eu]][ej_chose];
			//collect HE edges to delete:
			std::vector<std::vector<int>> edges_to_delete;
			std::vector<std::vector<int>> new_bdry_edges;
			std::vector<int> tria_to_delete;
			std::vector<std::pair<int, HE_edge*>> affected_vertices;


			HE_edge *checkE=Edges[KKey(eu,evvv)]->pair->next;

			std::vector<int> edge_ii; edge_ii.push_back(eu); edge_ii.push_back(evvv);
			std::vector<int> edge_iip; edge_iip.push_back(evvv); edge_iip.push_back(eu);
			edges_to_delete.push_back(edge_ii);
			edges_to_delete.push_back(edge_iip);
			edge_ii.clear(); edge_iip.clear();


			edge_ii.push_back(checkE->eu); edge_ii.push_back(checkE->ev);
			edge_iip.push_back(checkE->pair->eu); edge_iip.push_back(checkE->pair->ev);
			edges_to_delete.push_back(edge_ii);
			edges_to_delete.push_back(edge_iip);
			edge_ii.clear(); edge_iip.clear();

			edge_ii.push_back(checkE->next->eu); edge_ii.push_back(checkE->next->ev);
			new_bdry_edges.push_back(edge_ii);
			affected_vertices.push_back(std::make_pair(edge_ii[1],Edges[KKey(edge_ii[1],edge_ii[0])]));
			edge_ii.clear();

			edge_ii.push_back(checkE->next->pair->next->eu); edge_ii.push_back(checkE->next->pair->next->ev);
			affected_vertices.push_back(std::make_pair(edge_ii[0],Edges[KKey(edge_ii[0],edge_ii[1])]));
			edge_ii.clear();

			tria_to_delete.push_back(checkE->face->HE_face_id);



			while(checkE->pair->face!=NULL){
				checkE=checkE->pair->next;

				//Collect HE edges to delete
				edge_ii.push_back(checkE->eu); edge_ii.push_back(checkE->ev);
				edge_iip.push_back(checkE->pair->eu); edge_iip.push_back(checkE->pair->ev);
				edges_to_delete.push_back(edge_ii);
				edges_to_delete.push_back(edge_iip);
				edge_ii.clear(); edge_iip.clear();

				//Collect new Bdry HE's
				edge_ii.push_back(checkE->next->eu); edge_ii.push_back(checkE->next->ev);
				new_bdry_edges.push_back(edge_ii);
				edge_ii.clear();

				edge_ii.push_back(checkE->next->pair->next->eu); edge_ii.push_back(checkE->next->pair->next->ev);
				affected_vertices.push_back(std::make_pair(edge_ii[0],Edges[KKey(edge_ii[0],edge_ii[1])]));
				edge_ii.clear();

				//Collect triangles to delete
				tria_to_delete.push_back(checkE->face->HE_face_id);
				
			}

			//Re-calculate the HE structure after deleting the triangles\
			//Delete the triangle associated with the Pair HE of eu,ev.
			//And move the last triangle to the slot
			for(int triaiii=0;triaiii<e_tria_ct_min;triaiii++){
				int problem_tria_id=tria_to_delete[triaiii];
				int problem_tria_id3=3*problem_tria_id;
				int last_tria_id=tria_ct-1;
				int last_tria_id3=3*last_tria_id;
				tria_vertex[problem_tria_id3]=tria_vertex[last_tria_id3];
				tria_vertex[problem_tria_id3+1]=tria_vertex[last_tria_id3+1];
				tria_vertex[problem_tria_id3+2]=tria_vertex[last_tria_id3+2];

				//Update triangle ct
				tria_ct--;

				//Update the HE data strcuture: Faces, Edges, Vertices
				Faces[problem_tria_id].HE_face_id=last_tria_id;
				Faces[problem_tria_id].edge=Faces[last_tria_id].edge;
			}

			//Erase deleted triangle edges
			for(int edgeiii=0;edgeiii<edges_to_delete.size();edgeiii++){
				int deu=edges_to_delete[edgeiii][0];
				int dev=edges_to_delete[edgeiii][1];
				
				Edges.erase(KKey(deu,dev));
			}
			
			//Make new bdry HE
			for(int edgeiii=0;edgeiii<new_bdry_edges.size();edgeiii++){
				int aeu=new_bdry_edges[edgeiii][0];
				int aev=new_bdry_edges[edgeiii][1];
				
				Edges[KKey(aeu,aev)]->face=NULL;
				Edges[KKey(aeu,aev)]->prev=NULL;
				Edges[KKey(aeu,aev)]->next=NULL;
			}

			//Update the Vertices HE structure: Update one of the edge emanating from the vertex
			int ejj=ej_chose+1; if(ejj>=(signed int)eu_bdry_connect[eu_bdry_ind[eu]].size()){ejj=0;}
			Vertices[eu].edge=Edges[KKey(eu,eu_bdry_connect[eu_bdry_ind[eu]][ejj])];

			for(std::vector<std::pair<int, HE_edge*>>::iterator viii=affected_vertices.begin();viii!=affected_vertices.end();viii++){
				Vertices[viii->first].edge=viii->second;
			}

			have_deleted=true;

		}
		if(have_deleted==false){
			printf("Failed to resolve bdry triangles through HE\n");
			throw std::exception();
		}
		
	}
	
printf("C\n");


	//Obtain boundary edges connectivity: bdry HE prev and next
	//Loop through all edges and obtain the bdry HE's
	for(std::unordered_map<size_t, HE_edge*>::iterator itr_e=Edges.begin();itr_e!=Edges.end();itr_e++){
		if(Edges[itr_e->first]->face==NULL){
			int eu=itr_e->second->eu;
			int ev=itr_e->second->ev;
			Bdry_Edges[eu]=std::make_pair(ev,0); //0 signal not connected yet
		}
	}

	//clear memory of intermediate variables
	eu_bdry_connect.clear();
	eu_bdry_ind.clear();
	eu_problematic.clear(); 

printf("D\n");

	//Loop through the boundary half edges, and find the conenectivity
	for(std::unordered_map<unsigned int, std::pair<unsigned int, unsigned int>>::iterator itr_e=Bdry_Edges.begin();itr_e!=Bdry_Edges.end();itr_e++){
		int connected=itr_e->second.second;
		if(connected==0){
			int eu=itr_e->first;
			int ev=itr_e->second.first;
			Bdry_Edges_Start.push_back(std::make_pair(eu,ev));
			bool Continue_bdry=true;
			int eu0=eu; int ev0=ev;
			while(Continue_bdry){
				int nextu=ev;
				int nextv=Bdry_Edges[ev].first;
				Edges[KKey(eu,ev)]->next=Edges[KKey(nextu,nextv)];
				Edges[KKey(nextu,nextv)]->prev=Edges[KKey(eu,ev)];
				Bdry_Edges[eu].second=1;
				Bdry_Edges[nextu].second=1;
				eu=nextu;
				ev=nextv;
				if(eu==eu0){
					Continue_bdry=false;
				}
			}
		}
	}

	Bdry_Edges.clear();
	HE_exist=true;
printf("E\n");
}




/**
 * @brief Constructs a new wall_is_2d object.
 *
 * @param pm2d_ Pointer to an instance of the parallel_meshing_2d class.
 */
wall_is_2d::wall_is_2d(parallel_meshing_2d *pm2d_) : 
pm2d(pm2d_), tria_length_cri(new double[pm2d_->gnxy]), gnxy(pm2d_->gnxy),
gnx(pm2d_->gnx),gny(pm2d_->gny),
gdx(pm2d_->gdx),gdy(pm2d_->gdy),
inv_gdx(pm2d_->inv_gdx),inv_gdy(pm2d_->inv_gdy),

ax(pm2d_->ax),ay(pm2d_->ay),bx(pm2d_->bx),by(pm2d_->by),
cut_cell_base_octE_fac((1.0+1.0/(1.0+sqrt(2.0))))
{
	update_tria_length_cri();
}

/**
 * @brief Destructor for the wall_is_2d object.
 */
wall_is_2d::~wall_is_2d(){
	delete [] tria_length_cri;
}

/**
 * @brief Updates the triangle length criterion based on the mesh characteristic edge length.
 */
void wall_is_2d::update_tria_length_cri(){
	#pragma omp parallel num_threads(pm2d->num_t)
	for(int ij=0;ij<gnxy;ij++){
		tria_length_cri[ij]=voro_wall_fac * pm2d->chrtrt_len_h[ij];
	}
}

/**
 * @brief Updates the triangle length criterion for a specific grid cell.
 *
 * @param ij The index of the grid cell.
 * @param fac The factor to multiply the length criterion by.
 */
void wall_is_2d::update_tria_length_cri(int ij, double fac){
	tria_length_cri[ij]*=fac;
	int i=ij%gnx;
    int j=ij/gnx;
	//also, update wall of all grids within the new wall bound to be larger
	double half_new_length_cri=tria_length_cri[ij]*0.5;
    int gil=i-ceil(half_new_length_cri*inv_gdx); 
    int gih=i+ceil(half_new_length_cri*inv_gdx); 
    int gjl=j-ceil(half_new_length_cri*inv_gdy); 
    int gjh=j+ceil(half_new_length_cri*inv_gdy);  
    for(int gii=gil;gii<=gih;gii++){
    	for(int gjj=gjl;gjj<=gjh;gjj++){
    		if(gii>=0 && gii<gnx && gjj>=0 && gjj<gny){
    			int giijj=gnx*gjj+gii;
    			if(giijj!=ij){
    				tria_length_cri[giijj]=tria_length_cri[ij];
    			}
    		}
    	}
    }
}	

/**
 * @brief Cuts a cell with surrounding walls.
 *
 * @tparam v_cell_2d The type of the cell object.
 * @param c The cell object to be cut.
 * @param x The x-coordinate of the cell.
 * @param y The y-coordinate of the cell.
 * @return True if the cell is successfully cut, false otherwise.
 */
template<class v_cell_2d>
bool wall_is_2d::cut_cell_base(v_cell_2d &c,double x,double y) {

	double x_ax=x-ax; double y_ay=y-ay;
	double bx_x=bx-x; double by_y=by-y;
	int i=(x_ax)*inv_gdx;
	int j=(y_ay)*inv_gdy;
	int ij=gnx*j+i;
	double tria_length_cri_local=tria_length_cri[ij];
    double octR=tria_length_cri_local*0.5;
    double cell_ax=-octR; double cell_bx=octR;
    double cell_ay=-octR; double cell_by=octR; 
    if(x_ax<octR){cell_ax=-x_ax;}
    if(y_ay<octR){cell_ay=-y_ay;}
    if(bx_x<octR){cell_bx=bx_x;}
    if(by_y<octR){cell_by=by_y;}

    c.init(cell_ax,cell_bx,cell_ay,cell_by);
    double octE=cut_cell_base_octE_fac*octR;
    //cut the cell into an octagon shape
    c.nplane(-octE,-octE,-1);   //numbering the plane wall to be id -1
    c.nplane(-octE,octE,-1);
    c.nplane(octE,-octE,-1);
    c.nplane(octE,octE,-1);
    return true;

}




