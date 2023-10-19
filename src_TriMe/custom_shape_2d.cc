#include "custom_shape_2d.hh"

using namespace voro;

/**
 * @brief Constructor for the custom_shape_2d class.
 *
 * This constructor initializes a custom_shape_2d object with the provided parameters.
 *
 * @param boundaries_ The boundaries of the shape specified as a vector of vector of doubles.
 *                    Each inner vector represents a boundary, and the outer vector holds all the boundaries of the shape.
 * @param ax_ The minimum x-coordinate of the domain container.
 * @param bx_ The maximum x-coordinate of the domain container.
 * @param ay_ The minimum y-coordinate of the domain container.
 * @param by_ The maximum y-coordinate of the domain container.
 * @param normalize_model_ Flag indicating whether to normalize the model within the container domain.
 *                          If true, the shape will be centered and scaled to fit within
 *                          [ax + 0.1*(bx-ax), ax + 0.9*(bx-ax)] x [ay + 0.1*(by-ay), ay + 0.9*(by-ay)].
 * @param num_t_ The number of parallel threads to be used.
 */
custom_shape_2d::custom_shape_2d(std::vector<std::vector<double>> boundaries_,
    double ax_, double bx_, double ay_, double by_, 
    bool normalize_model_, int num_t_):ax(ax_), bx(bx_), 
    ay(ay_), by(by_), lx(bx-ax), ly(by-ay), 
    normalize_model(normalize_model_), num_t(num_t_)
    {
        //if normalize_model==true:
        //center shape boundary input and scale so that the shape
        //lies in [ax+0.1*(bx-ax), ax+0.9*(bx-ax)]x[ay+0.1*(by-ay),ay+0.9*(by-ay)]
        create_b_pts(boundaries_);
        get_cell_lineSeg();
    }

/**
 * @brief Destructor for the custom_shape_2d class.
 *
 * This destructor releases the dynamically allocated memory used by the custom_shape_2d object.
 */
custom_shape_2d::~custom_shape_2d(){
    for(int i=0; i<b_ct; i++){
        delete [] b_pts[i];
    }
    delete [] b_pts;
    delete [] seg_ct;
}


/**
 * @brief Read and store the boundary points of the shape.
 *
 * This function reads the boundary points provided as input and stores them in the `b_pts` array.
 * It also calculates the number of line segments for each boundary and initializes the `seg_ct` array.
 * If the `normalize_model` flag is set to true, it normalizes the shape by centering it and scaling it
 * to fit within the specified domain range.
 *
 * @param boundaries_ The boundary points of the shape.
 */
void custom_shape_2d::create_b_pts(std::vector<std::vector<double>> boundaries_){

    //read boundary points and store them in b_pts array
    b_ct=(signed int)boundaries_.size();
    b_pts=new double*[b_ct];
    seg_ct = new int[b_ct];
double t0=omp_get_wtime();
    for(int i=0; i<b_ct; i++){
        int isize=(signed int)boundaries_[i].size();
        b_pts[i]=new double[isize];

        #pragma omp parallel for num_threads(num_t)
        for(int j=0; j<isize; j++){
            b_pts[i][j]=boundaries_[i][j];
        }
        seg_ct[i]=isize/2-1;
    }
double t1=omp_get_wtime();
t_c1=t1-t0;

    if(normalize_model){
        //normalize:
        //center shape boundary input and scale so that the shape
        //lies in [ax+0.1*(bx-ax), ax+0.9*(bx-ax)]x[ay+0.1*(by-ay),ay+0.9*(by-ay)]
        double xmin=lx+ly;
        double xmax=-(lx+ly);
        double ymin=lx+ly;
        double ymax=-(lx+ly);
        
t0=omp_get_wtime();
        //loop through all boundaries
        for(int i=0;i<b_ct;i++){
            //loop through all line segments of the current boundary
            #pragma omp parallel for num_threads(num_t) reduction(max: xmax,ymax) reduction(min: xmin,ymin)
            for(int j=0;j<seg_ct[i];j++){
                double x0=b_pts[i][2*j];
                double y0=b_pts[i][2*j+1];
                double x1=b_pts[i][2*(j+1)];
                double y1=b_pts[i][2*(j+1)+1];

                if(x0<xmin){xmin=x0;}
                if(x1<xmin){xmin=x1;}

                if(y0<ymin){ymin=y0;}
                if(y1<ymin){ymin=y1;}

                if(x0>xmax){xmax=x0;}
                if(x1>xmax){xmax=x1;}

                if(y0>ymax){ymax=y0;}
                if(y1>ymax){ymax=y1;}
            }
        }
t1=omp_get_wtime();
t_c2=t1-t0;

        scale_max_range=xmax-xmin;
        if(ymax-ymin>scale_max_range){scale_max_range=ymax-ymin;}
        scale_xmid=0.5*(xmin+xmax);
        scale_ymid=0.5*(ymin+ymax);
        
        scale_min_domain_range=std::min(bx-ax,by-ay);
        
t0=omp_get_wtime();
        //normalize shape so that the model fix into domain box, biggest dimension of moedel is normalize_model_length_fac of the size of the minimum length side of the box
        //centered at the box midpoint.
        //loop through all boundaries
        for(int i=0;i<b_ct;i++){
            //loop through all points of the current boundary; 
            //for boundary i, pt_ct[i]=seg_ct[i]+1
            #pragma omp parallel for num_threads(num_t)
            for(int j=0;j<seg_ct[i]+1;j++){
                //scale and center x component of point
                b_pts[i][2*j]=scale_min_domain_range*normalize_model_length_fac*(b_pts[i][2*j]-scale_xmid)/scale_max_range+0.5*(bx-ax);
                //scale and center y component of point
                b_pts[i][2*j+1]=scale_min_domain_range*normalize_model_length_fac*(b_pts[i][2*j+1]-scale_ymid)/scale_max_range+0.5*(by-ay);
            }
        }
t1=omp_get_wtime();
t_c3=t1-t0;

    }


//---------------------print out the (normalized) geometry contour line segments---------------------
    /*
    char fb[256];
    sprintf(fb,"geo_bdry_line_segs.txt");
    FILE *fbout=fopen(fb,"a");
    
    for(int i=0;i<b_ct;i++){
            //loop through all points of the current boundary; 
            //for boundary i, pt_ct[i]=seg_ct[i]+1
            for(int j=0;j<seg_ct[i]+1;j++){
                //scale and center x component of point
                double linex=b_pts[i][2*j];
                //scale and center y component of point
                double liney=b_pts[i][2*j+1];
                fprintf(fbout,"%g %g \n", linex,liney);
            }
        }
    fclose(fbout);
    */
//----------------------------------------------------------------------------------------------------


    double lineSeg_len_avg=0.0;
    int lineSeg_ct=0;
t0=omp_get_wtime();
    for(int bi=0;bi<b_ct;bi++){
        #pragma omp parallel for num_threads(num_t) reduction(+: lineSeg_len_avg,lineSeg_ct)
        for(int sj=0;sj<seg_ct[bi];sj++){
            double x0=b_pts[bi][2*sj];
            double y0=b_pts[bi][2*sj+1];
            double x1=b_pts[bi][2*(sj+1)];
            double y1=b_pts[bi][2*(sj+1)+1];
            double lineSeg_len=d_points(x0,y0,x1,y1);
            lineSeg_len_avg+=lineSeg_len;
            lineSeg_ct+=1;
        }
    }
t1=omp_get_wtime();
t_c4=t1-t0;
    lineSeg_len_avg=lineSeg_len_avg/lineSeg_ct;
    nx=ceil(lx/lineSeg_len_avg);ny=ceil(ly/lineSeg_len_avg);
    dx=lx/nx;dy=ly/ny;
    inv_dx=1.0/dx;inv_dy=1.0/dy;
//printf("lineSeg ct: %d len avg: %g nx, ny: %d %d\n",lineSeg_ct, lineSeg_len_avg, nx, ny);
    
}


/**
 * @brief Calculate the list of boundary line segments intersect with each grid cell.
 *
 * This function populates the `cell_lineSeg` vector, which contains the list of boundary line segments intersect with each grid cell.
 * The line segments are represented by pairs of indices:
 * the index of the boundary (`bi`) and the index of the line segment within that boundary (`sj`).
 * The `cell_lineSeg` vector is indexed by `ij=nx*j+i` the grid index.
 */
void custom_shape_2d::get_cell_lineSeg(){
    cell_lineSeg.resize(nx*ny);
double t0=omp_get_wtime();
    //loop through each bdry line segment
    for(int bi=0;bi<b_ct;bi++){
        #pragma omp parallel for num_threads(num_t)
        for(int sj=0;sj<seg_ct[bi];sj++){
            double x0=b_pts[bi][2*sj];
            double y0=b_pts[bi][2*sj+1];
            double x1=b_pts[bi][2*(sj+1)];
            double y1=b_pts[bi][2*(sj+1)+1];

            //the grids that the end points lie in
            int ii0=(x0-ax)*inv_dx; int jj0=(y0-ay)*inv_dy;
            int ii1=(x1-ax)*inv_dx; int jj1=(y1-ay)*inv_dy;

            //find the minimum bounding box of the line segment
            double xmin=x0; double xmax=x1; if(x1<xmin){xmin=x1;xmax=x0;}
            double ymin=y0; double ymax=y1; if(y1<ymin){ymin=y1;ymax=y0;}
            int iil=(xmin-ax)*inv_dx-1; int iih=(xmax-ax)*inv_dx+1;  if(iil<0){iil=0;} if(iih>=nx){iih=nx-1;} 
            int jjl=(ymin-ay)*inv_dy-1; int jjh=(ymax-ay)*inv_dy+1;  if(jjl<0){jjl=0;} if(jjh>=ny){jjh=ny-1;} 
        
            //loop through all grids of the minimum bounding box,
            //and test if the grid intersect with the line segment
            for(int ii=iil;ii<=iih;ii++){
                for(int jj=jjl;jj<=jjh;jj++){
                    //bmin
                    double bminx=ax+dx*ii;
                    double bminy=ay+dy*jj;
                    
                    //bmax
                    double bmaxx=bminx+dx;
                    double bmaxy=bminy+dy;

                    bool intersection=true;

                    //test if the grid holds either endpoint
                    if(ii==ii0 && jj==jj0){}
                    else if(ii==ii1 && jj==jj1){}
                    else{
                        //test if bdry line segment intersect with either of the four line segments of the grid cell
                        if(line_seg_intersect(x0,y0, x1,y1, bminx,bminy, bminx,bmaxy)==false){
                            if(line_seg_intersect(x0,y0, x1,y1, bminx,bminy, bmaxx,bminy)==false){
                                if(line_seg_intersect(x0,y0, x1,y1, bmaxx,bminy, bmaxx,bmaxy)==false){
                                    if(line_seg_intersect(x0,y0, x1,y1, bmaxx,bmaxy, bminx,bmaxy)==false){
                                        intersection=false;
                                    }
                                }
                            }
                        }
                    }

                    if(intersection==true){
                        #pragma omp critical
                        {
                        cell_lineSeg[nx*jj+ii].push_back(std::make_pair(bi,sj));
                        }
                    }

                }
            }

        }
    }
double t1=omp_get_wtime();
t_c5=t1-t0;
}



/**
 * @brief Calculates the signed distance of a point (x, y) to the boundary of a custom shape in 2D.
 * 
 * @param x The x-coordinate of the point.
 * @param y The y-coordinate of the point.
 * @return The signed distance of the point to the shape boundary.
 */
double custom_shape_2d::f_b_pts(double x, double y){

    //the grid cell(pi,pj) that the point is in
    int pi=(x-ax)*inv_dx; if(x>=bx-1e-10&&x<=bx+1e-10){pi=nx-1;} 
    int pj=(y-ay)*inv_dy; if(y>=by-1e-10&&y<=by+1e-10){pj=ny-1;}
    int pip=pi; int pim=pi; int pjp=pj; int pjm=pj; 
    
    std::vector< std::pair <int,int> > lineSeg_list;
    std::vector< std::pair <int,int> > lineSeg_list_unique;
    bool region_lineSeg=false;
    int cell_lineSeg_ind=nx*pj+pi;
    //if the cell the point is in has lineSegs
    if(cell_lineSeg[cell_lineSeg_ind].size()!=0){
        region_lineSeg=true;
        //store the lineSeg index
        for(int lineSegi=0;lineSegi<cell_lineSeg[cell_lineSeg_ind].size(); lineSegi++){
            lineSeg_list.push_back(cell_lineSeg[cell_lineSeg_ind][lineSegi]);
        }
    }
    int pim_prev=pim; int pip_prev=pip; int pjm_prev=pjm; int pjp_prev=pjp;
    //the increment layer to test to see if there are lineSeg
    int incm=0;
    while(region_lineSeg==false){
        incm++;
        pip=pi+incm; pim=pi-incm; 
        pjp=pj+incm; pjm=pj-incm; 
        if(pip>=nx && pim<0 && pjp>=ny && pjm<0){
            printf("ERROR: custom_shape_2d: no bdry line segment found\n");
            throw std::exception();
        }
        if(pip>=nx){pip=nx-1;} if(pim<0){pim=0;}
        if(pjp>=ny){pjp=ny-1;} if(pjm<0){pjm=0;}
        for(int ii=pim;ii<=pip;ii++){
            for(int jj=pjm;jj<=pjp;jj++){
                
                //only need to consider the outer new layer
                if(ii<pim_prev || ii>pip_prev || jj<pjm_prev || jj>pjp_prev){
                    cell_lineSeg_ind=nx*jj+ii;
                    if(cell_lineSeg[cell_lineSeg_ind].size()!=0){
                        region_lineSeg=true;
                        //store the lineSeg index
                        for(int lineSegi=0;lineSegi<cell_lineSeg[cell_lineSeg_ind].size();lineSegi++){
                            lineSeg_list.push_back(cell_lineSeg[cell_lineSeg_ind][lineSegi]);
                        }
                    }
                }
            }
            
        }
        pim_prev=pim; pip_prev=pip; pjm_prev=pjm; pjp_prev=pjp;
    }

    //Now region_lineSeg==true
    //get lineSegs of the bounding box of the previous region square
    double fac0=incm+1.0;
    
    double aa=fac0*dx; double bb=fac0*dy;
    
    //double a=std::max(aa,bb);
    //double sphere_r=sqrt(2)*a;
    
    double sphere_r=sqrt(aa*aa+bb*bb);
    
    double cellx0=ax+pi*dx; double cellx1=cellx0+dx;
    double celly0=ay+pj*dy; double celly1=celly0+dy;
    
    int ih=(cellx1+sphere_r-ax)*inv_dx; if(ih>=nx){ih=nx-1;}
    int il=(cellx0-sphere_r-ax)*inv_dx; if(il<0){il=0;}
    int jh=(celly1+sphere_r-ay)*inv_dy; if(jh>=ny){jh=ny-1;}
    int jl=(celly0-sphere_r-ay)*inv_dy; if(jl<0){jl=0;}
    
    for(int ii=il;ii<=ih;ii++){
        for(int jj=jl;jj<=jh;jj++){

            //only need to consider the outer new layer
            if(ii<pim || ii>pip || jj<pjm ||jj>pjp ){
                cell_lineSeg_ind=nx*jj+ii;
                if(cell_lineSeg[cell_lineSeg_ind].size()!=0){
                    //store the lineSeg index
                    for(int lineSegi=0;lineSegi<cell_lineSeg[cell_lineSeg_ind].size();lineSegi++){
                        lineSeg_list.push_back(cell_lineSeg[cell_lineSeg_ind][lineSegi]);
                    }
                }
            }
        
        }
    }

    //sort the lineSeg_list and delete the duplicated ones, and get a unique list of lineSegs
    //default sort into an ascending array of lineSeg bi-sj's,
    //sorting first by bi, then by sj
    std::sort(lineSeg_list.begin(), lineSeg_list.end());
    //delete duplicated ones 
    lineSeg_list_unique.push_back(lineSeg_list[0]);
    for(int i=0;i<lineSeg_list.size();i++){
        if(lineSeg_list_unique.back()!=lineSeg_list[i]){
            lineSeg_list_unique.push_back(lineSeg_list[i]);
        }
    }
    lineSeg_list.clear();

    double d=lx+ly;
    int d_status;
    int d_seg[2]; //track which line segment is the one with shortest distance
    double p_lx; 
    double p_ly; 

    double closestx,closesty; int status;
    //loop through the unique list of lineSegs 
    for(int li=0; li<lineSeg_list_unique.size(); li++){
        int i=lineSeg_list_unique[li].first;
        int j=lineSeg_list_unique[li].second;

        double x0=b_pts[i][2*j];
        double y0=b_pts[i][2*j+1];
        double x1=b_pts[i][2*(j+1)];
        double y1=b_pts[i][2*(j+1)+1];
        double d_temp=f_line_seg(x,y,x0,y0,x1,y1,closestx,closesty,status);

        if(d_temp<d){
            d_seg[0]=i;  //which boundary
            d_seg[1]=j;  //which line segment of the boundary
            d_status=status; //end points or in between
            d=d_temp; //shortest distance
            p_lx=closestx; //point on line segment closest to pt, x component
            p_ly=closesty; //point on line segment closest to pt, y component
        }

    }




/*
    double closestx,closesty; int status;
    //loop through all boundaries
    for(int i=0;i<b_ct;i++){
        //loop through all line segments of the current boundary
        for(int j=0;j<seg_ct[i];j++){
            double x0=b_pts[i][2*j];
            double y0=b_pts[i][2*j+1];
            double x1=b_pts[i][2*(j+1)];
            double y1=b_pts[i][2*(j+1)+1];
            double d_temp=f_line_seg(x,y,x0,y0,x1,y1,closestx,closesty,status);

            if(d_temp<d){
                d_seg[0]=i;  //which boundary
                d_seg[1]=j;  //which line segment of the boundary
                d_status=status; //end points or in between
                d=d_temp; //shortest distance
                p_lx=closestx; //point on line segment closest to pt, x component
                p_ly=closesty; //point on line segment closest to pt, y component
            }
            
        }
    }
*/
    
    //decide the sign of the shortest distance
    //the line segment looking at
    double x0=b_pts[d_seg[0]][2*d_seg[1]];
    double y0=b_pts[d_seg[0]][2*d_seg[1]+1];
    double x1=b_pts[d_seg[0]][2*(d_seg[1]+1)];
    double y1=b_pts[d_seg[0]][2*(d_seg[1]+1)+1];
    
    //the normal of the line segment p0p1 based on the right hand rule
    double nx=y1-y0;
    double ny=x0-x1;
    //unitize the normal vector
    double nlen=sqrt(nx*nx+ny*ny);
    nx=nx/nlen;
    ny=ny/nlen;
    
    if(d_status==1){ //closest point is in between the line segment
        //p_l to p
        double plpx=x-p_lx;
        double plpy=y-p_ly;
        //plp dot n: if positive, the same side, inside geometry, -; 
        //if negative, opposite side, outside geometry, +
        double dot_temp=plpx*nx+plpy*ny;
        if (dot_temp>0){
            d=-d;
        }
    }
    else if(d_status==0){ //closest point is p0
        //get the other connecting bar of p0: p00p0
        double x00;
        double y00;
        if(d_seg[1]==0){
            x00=b_pts[d_seg[0]][2*(seg_ct[d_seg[0]]-1)];
            y00=b_pts[d_seg[0]][2*(seg_ct[d_seg[0]]-1)+1];
        }
        else{
            x00=b_pts[d_seg[0]][2*(d_seg[1]-1)];
            y00=b_pts[d_seg[0]][2*(d_seg[1]-1)+1];
        }
        //the normal of the line segment p00p0 based on the right hand rule
        double nx0=y0-y00;
        double ny0=x00-x0;
        //unitize the normal vector
        double nlen0=sqrt(nx0*nx0+ny0*ny0);
        nx0=nx0/nlen0;
        ny0=ny0/nlen0;
        //average of the two normal of the two connecting line segment
        //used as the "normal" at point p0
        double nx_p0=0.5*(nx+nx0);
        double ny_p0=0.5*(ny+ny0);
        //p_0 to p
        double p0px=x-x0;
        double p0py=y-y0;
        //p0p dot n: if positive, the same side, inside geometry, -; 
        //if negative, opposite side, outside geometry, +
        double dot_temp=p0px*nx_p0+p0py*ny_p0;
        if (dot_temp>0){
            d=-d;
        }
    }
    else { //d_status==2: closest point is p1
        //get the other connecting bar of p1: p1p11
        double x11;
        double y11; 
        if(d_seg[1]==seg_ct[d_seg[0]]-1){
            x11=b_pts[d_seg[0]][2];
            y11=b_pts[d_seg[0]][3];
        }
        else{
            x11=b_pts[d_seg[0]][2*(d_seg[1]+2)];
            y11=b_pts[d_seg[0]][2*(d_seg[1]+2)+1];
        }
        //the normal of the line segment p1p11 based on the right hand rule
        double nx1=y11-y1;
        double ny1=x1-x11;
        //unitize the normal vector
        double nlen1=sqrt(nx1*nx1+ny1*ny1);
        nx1=nx1/nlen1;
        ny1=ny1/nlen1;
        //average of the two normal of the two connecting line segment
        //used as the "normal" at point p0
        double nx_p1=0.5*(nx+nx1);
        double ny_p1=0.5*(ny+ny1);
        //p_1 to p
        double p1px=x-x1;
        double p1py=y-y1;
        //p0p dot n: if positive, the same side, inside geometry, -; 
        //if negative, opposite side, outside geometry, +
        double dot_temp=p1px*nx_p1+p1py*ny_p1;
        if (dot_temp>0){
            d=-d;
        }
    }

    return d;
}






