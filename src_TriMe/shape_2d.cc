#include "shape_2d.hh"

using namespace voro;

/**
 * @brief Constructor for the shape_2d class.
 *
 * This constructor sets up the geometry grid based on a factor of the dimensions of the container grid.
 *
 * @param con_ The container_2d object.
 * @param num_t_ Number of parallel threads.
 */
shape_2d::shape_2d(container_2d &con_, int num_t_)
	: num_t(num_t_), 
	ax(con_.ax), bx(con_.bx), ay(con_.ay), by(con_.by), 
    is_custom_shape_contour(false), shape_scaling(false),

	geo_grid(new int[1]), 
    geo_igrid_ij(new int[1]), 
    geo_bgrid_ij(new int[1]),
    geo_ogrid_ij(new int[1]),
	geo_bgrid_ct(0), geo_igrid_ct(0), geo_ogrid_ct(0)
	{
        gnx=geo_grid_nxy_fac_con_grid*con_.nx;
        gny=geo_grid_nxy_fac_con_grid*con_.ny;
        gnxy=gnx*gny;
        gdx=(bx-ax)/gnx;
        gdy=(by-ay)/gny;
        inv_gdx=1.0/gdx;
        inv_gdy=1.0/gdy;
        diag_gdxy=sqrt(gdx*gdx+gdy*gdy);

    }

/**
 * @brief Destructor for the shape_2d class.
 *
 * This destructor frees the memory allocated for the geometry grid arrays.
 */
shape_2d::~shape_2d(){
	delete [] geo_grid;
	delete [] geo_bgrid_ij;
	delete [] geo_igrid_ij;
    delete [] geo_ogrid_ij;
}


/**
 * @brief Generates the geometry grid, which categorize each grid cell as inner/boundary/outer grid:  -1,-2,-3: inside; 1,2,3: boundary; gnxy+1,gnxy+2,...: outside.
 *
 * This function sets up the geometry grid by looping through the cells of the grid,
 * checking the mid-point signed distance values, and determining if the cell is
 * an inner cell, boundary cell, or outer cell within a bandwidth of geometry.
 * It assigns unique identifiers to the cells based on their type and stores them
 * in the corresponding geometry grid arrays.
 */
void shape_2d::get_geometryGrid(){
    
    delete [] geo_grid;
    geo_grid=new int[gnxy];
    geo_bgrid_ct=0;
    geo_igrid_ct=0;  
    geo_ogrid_ct=0;  

    //loop through the cells, check the mid-point fd values, and determine if the cell is inner cell or boundary within a band grid of geometry
    #pragma omp parallel for num_threads(num_t) 
    for(int ij=0; ij<gnxy; ij++){
        
        int i=ij%gnx;
        int j=ij/gnx;
        
        //mid point inside grid, find its signed distance
        double x=ax+(0.5+i)*gdx; 
        double y=ay+(0.5+j)*gdy;
        double fd_ij=sdf(x,y); //Use sdf(). construct adf() later on bdry grids.

        //>0, ind+1: bdry grid: sdf of midpoint |d|<=diag_gdxy
        //<0, -(ind+1): inner grid
        //outer grid: gnxy+1,gnxy+2...
        if(abs(fd_ij)<=diag_gdxy){
        	#pragma omp atomic capture
            geo_grid[ij]=++geo_bgrid_ct;  //1,2,3,...bdry grid
        }
        else if(fd_ij<-diag_gdxy){
        	#pragma omp atomic capture
        	geo_grid[ij]=++geo_igrid_ct; //1,2,3,... inner grid

            geo_grid[ij]=-(geo_grid[ij]); //-1,-2,-3,... 
        }
        else{ //outer grid: gnxy+1,gnxy+2,...
            #pragma omp atomic capture
            geo_grid[ij]=++geo_ogrid_ct;

            geo_grid[ij]=geo_grid[ij]+gnxy;
        }
    }

    if(geo_bgrid_ct==0 || geo_igrid_ct==0){
    	printf("ERROR: geometry bdry and/or inner grids ct is 0.\n");
    	throw std::exception();
    }

    delete [] geo_bgrid_ij;
    geo_bgrid_ij=new int[geo_bgrid_ct];

    delete [] geo_igrid_ij;
    geo_igrid_ij=new int[geo_igrid_ct];

    delete [] geo_ogrid_ij;
    geo_ogrid_ij=new int[geo_ogrid_ct];

    //get the inner and bdry grid ij's arrays
    #pragma omp parallel for num_threads(num_t)
    for(int ij=0; ij<gnxy; ij++){

        if(geo_grid[ij]>gnxy){ //outer grid
            int ind=geo_grid[ij]-gnxy-1;
            geo_ogrid_ij[ind]=ij;
        }
        
        if(geo_grid[ij]>0 && geo_grid[ij]<=gnxy){ //bdry grid
        	int ind=geo_grid[ij]-1;
        	geo_bgrid_ij[ind]=ij;
        }

        if(geo_grid[ij]<0){ //inner grid
        	int ind=abs(geo_grid[ij])-1;
        	geo_igrid_ij[ind]=ij;
        }
    }

}


/**
 * @brief Determines the category of the geometry grid that a point lies in.
 *
 * This function calculates the indices of the grid cell containing the point and
 * retrieves the corresponding geometry grid category. The categories are as follows:
 *   - Inner Grid: grid_ctgr < 0
 *   - Boundary Grid: grid_ctgr > 0 && grid_ctgr <= gnxy
 *   - Outer Grid: grid_ctgr > gnxy
 *
 * @param x The x-coordinate of the point.
 * @param y The y-coordinate of the point.
 * @return The category of the geometry grid that the point lies in.
 */
int shape_2d::pt_geo_grid_val(double x, double y){
	//find the geometry grid the pt lies in
    int i=int((x-ax)*inv_gdx); if(fabs(x-bx)<1e-14){i=gnx-1;} //if point at near upper boundary, count as the last grid
    int j=int((y-ay)*inv_gdy); if(fabs(y-by)<1e-14){j=gny-1;}
    if(i<0 || i>=gnx || j<0 || j>=gny){
        //(x,y) outside range. 
        printf("ERROR: (%g,%g) outside range, not in any geometry grid. Count as outer grid.\n", x,y);
        return 3*gnxy;
    } 
    int ij=gnx*j+i;
    int grid_ctgr=geo_grid[ij];
    return grid_ctgr;
}



/**
 * @brief Prints the geometry grid and signed distance field (SDF) values to separate files.
 *
 * This function prints the geometry grid values and corresponding SDF values for each grid cell
 * to separate text files. The geometry grid values are written to a file named "geo_grid.txt"
 * and the SDF values are written to a file named "sdf.txt".
 *
 * @param case_name The name of the case or directory where the files will be saved.
 */
void shape_2d::print_geo_grid_sdf_to_file(const char *case_name){
    char fg[256];
    sprintf(fg,"%s/geo_grid.txt",case_name);
    FILE *fgout=fopen(fg,"a");

    char fsd[256];
    sprintf(fsd,"%s/sdf.txt",case_name);
    FILE *fsdout=fopen(fsd,"a");
    
    for(int i=0;i<gnx;i++){
        for(int j=0;j<gny;j++){

            double x=ax+(0.5+i)*gdx; 
            double y=ay+(0.5+j)*gdy;
            int ij=gnx*j+i;

            fprintf(fgout,"%d ",geo_grid[ij]);
            fprintf(fsdout,"%g ",sdf(x,y));
        }
        fprintf(fgout,"\n");
        fprintf(fsdout,"\n");
    }
    fclose(fgout);
    fclose(fsdout);
}



