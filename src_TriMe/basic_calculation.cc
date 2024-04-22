#include "basic_calculation.hh"

using namespace voro;


/**
 * @brief Calculates the unsigned distance from a point to a line segment.
 *
 * This function calculates the distance from a point (x, y) to a line segment defined by
 * two endpoints (x0, y0) and (x1, y1). It also determines the closest point on the line
 * segment to the given point and provides additional information about the relative position
 * of the closest point.
 *
 * @param x The x-coordinate of the point.
 * @param y The y-coordinate of the point.
 * @param x0 The x-coordinate of the first endpoint of the line segment.
 * @param y0 The y-coordinate of the first endpoint of the line segment.
 * @param x1 The x-coordinate of the second endpoint of the line segment.
 * @param y1 The y-coordinate of the second endpoint of the line segment.
 * @param [out] closestx The x-coordinate of the closest point on the line segment.
 * @param [out] closesty The y-coordinate of the closest point on the line segment.
 * @param [out] status A reference to an integer variable that indicates the position of the
 *                     closest point: 0 for the first endpoint, 1 for between the endpoints,
 *                     and 2 for the second endpoint.
 * @return The unsigned distance from the point to the line segment.
 */
double basic_calculation_2d::f_line_seg(double x, double y, 
	double x0, double y0, double x1, double y1, 
	double &closestx, double &closesty, int &status){

    double d; 
    
    //the line segment vector v=p1-p0
    double vx=x1-x0;
    double vy=y1-y0;
    //the Point p to p0 vector w=p-p0
    double wx=x-x0;
    double wy=y-y0;


    
    //c1=w.dot(v)
    double c1=wx*vx+wy*vy;
    if(c1<=0){ //before p0
        d=d_points(x, y, x0, y0);
        closestx=x0;
        closesty=y0;
        status=0;
        return d;
    }
    else{
        //c2=v.dot(v)
        double c2=vx*vx+vy*vy;
        
        if(c2<=c1){ //after p1
            d=d_points(x, y, x1, y1);
            closestx=x1;
            closesty=y1;
            status=2;
            return d;
        }
        else{
            double b=c1/c2;
            //closest (perpendicular) point on line segment to point p
            double p_lx=x0+b*vx;
            double p_ly=y0+b*vy;
            d=d_points(x, y, p_lx, p_ly);
            closestx=p_lx;
            closesty=p_ly;
            status=1;
            return d;
        }
    }
}

/**
 * @brief Calculates the unsigned distance from a point to a line segment.
 *
 * This function calculates the distance from a point (x, y) to a line segment defined by
 * two endpoints (x0, y0) and (x1, y1).
 *
 * @param x The x-coordinate of the point.
 * @param y The y-coordinate of the point.
 * @param x0 The x-coordinate of the first endpoint of the line segment.
 * @param y0 The y-coordinate of the first endpoint of the line segment.
 * @param x1 The x-coordinate of the second endpoint of the line segment.
 * @param y1 The y-coordinate of the second endpoint of the line segment.
 * @return The unsigned distance from the point to the line segment.
 */
double basic_calculation_2d::f_line_seg(double x, double y, 
    double x0, double y0, double x1, double y1){

    double d; 
    
    //the line segment vector v=p1-p0
    double vx=x1-x0;
    double vy=y1-y0;
    //the Point p to p0 vector w=p-p0
    double wx=x-x0;
    double wy=y-y0;

    //c1=w.dot(v)
    double c1=wx*vx+wy*vy;
    if(c1<=0){ //before p0
        d=d_points(x, y, x0, y0);
        return d;
    }
    else{
        //c2=v.dot(v)
        double c2=vx*vx+vy*vy;
        
        if(c2<=c1){ //after p1
            d=d_points(x, y, x1, y1);
            return d;
        }
        else{
            double b=c1/c2;
            //closest (perpendicular) point on line segment to point p
            double p_lx=x0+b*vx;
            double p_ly=y0+b*vy;
            d=d_points(x, y, p_lx, p_ly);
            return d;
        }
    }
}


/**
 * @brief Projects a point onto the geometry boundary.
 *
 * This function projects a point (xx1, yy1) onto the boundary of the geometry. The projection
 * is performed by iteratively moving the point towards the boundary using Newton's method.
 * The function takes into account the specified deps_prime and geps_prime values, and the
 * projection behavior can be controlled by the project_pt_outside_geo parameter.
 *
 * @param[in,out] new_x The x-coordinate of the projected point. Upon successful projection,
 *                      it will be updated with the x-coordinate of the projected point.
 * @param[in,out] new_y The y-coordinate of the projected point. Upon successful projection,
 *                      it will be updated with the y-coordinate of the projected point.
 * @param xx1 The x-coordinate of the point to be projected.
 * @param yy1 The y-coordinate of the point to be projected.
 * @param deps_prime The value of deps_prime, which is the default deps value to use if no valid deps
 *                   is available. It is used to determine the step size for Newton's method.
 * @param geps_prime The value of geps_prime, which is the default geps value to use if no valid geps
 *                   is available. It represents the geometry boundary tolerance used in the stop criteria of
 *                   Newton's method. If the distance to the boundary is less than geps_prime, the
 *                   projection is considered successful.
 * @param project_pt_outside_geo Specifies whether to project points outside the geometry boundary
 *                               or not. If set to true, only points outside the boundary will be
 *                               projected. If set to false, all points (inside and outside) not on
 *                               the boundary will be projected.
 * @return True if the projection is successful and the point lies on the boundary, false otherwise.
 *         If the point lies in the outer grid, the function immediately returns false.
 */
bool basic_calculation_2d::pt_projection(double &new_x, double &new_y,
     double xx1, double yy1, 
     double deps_prime, double geps_prime,
     bool project_pt_outside_geo)
{
    new_x=xx1;
    new_y=yy1;
    double deps=deps_func(xx1,yy1,deps_prime);
    double geps=geps_func(xx1,yy1,geps_prime);
    double disk=sdf_func(xx1,yy1);
    bool do_projection=false;
    if(project_pt_outside_geo==true){
        if(disk>geps){ //points outside   
            do_projection=true;
        }
    }
    else{ //project any pt (inner, outside) not on bdry onto bdry
        if(abs(disk)>geps){
            do_projection=true;
        }
    }
    
    if(do_projection){ 
        int dnewtonCt =0;
        double xxk=xx1; 
        double yyk=yy1; 
        double mul=1.0/(deps*deps);
        double mul2=1.0/(2.0*deps);
        //Stop criteria: abs(sdf_func(xxk,yyk))<geps, i.e. when new point is on geo bdry
        while (abs(disk) > geps && dnewtonCt < max_projection_step){
            dnewtonCt ++;
            double fd_pdelx_y=sdf_func(xxk+deps,yyk);
            double fd_x_pdely=sdf_func(xxk,yyk+deps);
            double fd_mdelx_y=sdf_func(xxk-deps,yyk);
            double fd_x_mdely=sdf_func(xxk,yyk-deps);

            //Numerical gradient: Center difference approximation; 
            double gradxx=(fd_pdelx_y-fd_mdelx_y)*mul2; //fx
            double gradyy=(fd_x_pdely-fd_x_mdely)*mul2; //fy
            double gradxxxx=(fd_pdelx_y-2*disk+fd_mdelx_y)*mul; //fxx
            double gradyyyy=(fd_x_pdely-2*disk+fd_x_mdely)*mul; //fyy
            double gradxxyy=(sdf_func(xxk+deps,yyk+deps)-sdf_func(xxk+deps,yyk-deps)-sdf_func(xxk-deps,yyk+deps)+sdf_func(xxk-deps,yyk-deps))*mul*(1.0/4.0); //fxy

            //Note: la=disk;ja=gradxx;jb=gradyy;
            double lb=(xxk-xx1)*gradyy-(yyk-yy1)*gradxx;
            double jc=gradyy+(xxk-xx1)*gradxxyy-(yyk-yy1)*gradxxxx;
            double jd=-gradxx-(yyk-yy1)*gradxxyy+(xxk-xx1)*gradyyyy;
            double detJ=gradxx*jd-gradyy*jc;

            double mul5=-alpha*1.0/detJ;
            double cx=mul5*(jd*disk-gradyy*lb);
            double cy=mul5*(-jc*disk+gradxx*lb);
            xxk=xxk+cx;
            yyk=yyk+cy;

            if(pt_in_outer_grid(xxk,yyk)==true){ //outer grid, return false immediately and count as diverge.
                return false;
            }

            disk=sdf_func(xxk, yyk); 
            deps=deps_func(xxk,yyk,deps);
            geps=geps_func(xxk,yyk,geps);
        }
        new_x=xxk;
        new_y=yyk;
    }
    if(abs(disk) > geps){ //divergence, fail to find pt on bdry
        return false;
    }
    else{
        return true; //new pt lies on bdry
    }
}





/**
 * @brief Calculates the distance and closest point information between a point and a line segment.
 *
 * This function computes the distance (dis_temp) between a point P0(x0, y0, z0) and a line segment
 * defined by two endpoints (lx1, ly1, lz1) and (lx2, ly2, lz2). It also determines the closest
 * point (x_temp, y_temp, z_temp) on the line segment to the given point. The status_temp variable
 * indicates the relative position of the closest point on the line segment:
 *   - 0: The closest point is the first endpoint (lx1, ly1, lz1).
 *   - 1: The closest point is between the two endpoints.
 *   - 2: The closest point is the second endpoint (lx2, ly2, lz2).
 *
 * @param x0 The x-coordinate of the point P0.
 * @param y0 The y-coordinate of the point P0.
 * @param z0 The z-coordinate of the point P0.
 * @param lx1 The x-coordinate of the first endpoint of the line segment.
 * @param ly1 The y-coordinate of the first endpoint of the line segment.
 * @param lz1 The z-coordinate of the first endpoint of the line segment.
 * @param lx2 The x-coordinate of the second endpoint of the line segment.
 * @param ly2 The y-coordinate of the second endpoint of the line segment.
 * @param lz2 The z-coordinate of the second endpoint of the line segment.
 * @param[out] x_temp The x-coordinate of the closest point on the line segment.
 * @param[out] y_temp The y-coordinate of the closest point on the line segment.
 * @param[out] z_temp The z-coordinate of the closest point on the line segment.
 * @param[out] status_temp The status indicating the relative position of the closest point.
 * @param[out] dis_temp The distance between the point and the line segment.
 */
void basic_calculation_3d::p_lineSeg_dis(double x0, double y0, double z0, 
	double lx1, double ly1, double lz1, double lx2, double ly2, double lz2, 
    double &x_temp, double &y_temp, double &z_temp, int &status_temp, double &dis_temp){
    
    double d;
    //the line segment vector v=lp2-lp1
    double vx=lx2-lx1;
    double vy=ly2-ly1;
    double vz=lz2-lz1;
    //the Point p to p0 vector w=p0-lp1
    double wx=x0-lx1;
    double wy=y0-ly1;
    double wz=z0-lz1;
    
    //c1=w.dot(v)
    double c1=wx*vx+wy*vy+wz*vz;
    if(c1<=0+1e-7){ //before lp1
        d=d_points(x0, y0, z0, lx1, ly1, lz1);
        x_temp=lx1;
        y_temp=ly1;
        z_temp=lz1;
        status_temp=0;
        dis_temp=d;
        return;
    }
    else{
        //c2=v.dot(v)
        double c2=vx*vx+vy*vy+vz*vz;
        
        if(c2<=c1+1e-7){ //after lp2
            d=d_points(x0, y0, z0, lx2, ly2, lz2);
            x_temp=lx2;
            y_temp=ly2;
            z_temp=lz2;
            status_temp=2;
            dis_temp=d;
            return;
        }
        else{
            double b=c1/c2;
            //closest (perpendicular) point on line segment to point p
            double p_lx=lx1+b*vx;
            double p_ly=ly1+b*vy;
            double p_lz=lz1+b*vz;
            d=d_points(x0, y0, z0, p_lx, p_ly, p_lz);
            x_temp=p_lx;
            y_temp=p_ly;
            z_temp=p_lz;
            status_temp=1;
            dis_temp=d;
            return;
        }
    }
}

/**
 * @brief Computes the unsigned distance from a point to a triangle and determines the closest point.
 *
 * This function calculates the unsigned distance (d) from a point P0(x0, y0, z0) to a triangle ABC,
 * with an outward normal (a, b, c). It also determines the closest point (closestx, closesty, closestz)
 * on the triangle to the given point and provides a status code (status) indicating the relative position
 * of the closest point:
 *   - 0: The closest point is on the line segment AB.
 *   - 1: The closest point is on the line segment BC.
 *   - 2: The closest point is on the line segment CA.
 *   - 3: The closest point is on vertex A.
 *   - 4: The closest point is on vertex B.
 *   - 5: The closest point is on vertex C.
 *   - 6: The closest point is inside the triangle.
 *
 * @param x0 The x-coordinate of the point P0.
 * @param y0 The y-coordinate of the point P0.
 * @param z0 The z-coordinate of the point P0.
 * @param xA The x-coordinate of vertex A of the triangle.
 * @param yA The y-coordinate of vertex A of the triangle.
 * @param zA The z-coordinate of vertex A of the triangle.
 * @param xB The x-coordinate of vertex B of the triangle.
 * @param yB The y-coordinate of vertex B of the triangle.
 * @param zB The z-coordinate of vertex B of the triangle.
 * @param xC The x-coordinate of vertex C of the triangle.
 * @param yC The y-coordinate of vertex C of the triangle.
 * @param zC The z-coordinate of vertex C of the triangle.
 * @param a The x-component of the outward normal of the triangle.
 * @param b The y-component of the outward normal of the triangle.
 * @param c The z-component of the outward normal of the triangle.
 * @param[out] closestx The x-coordinate of the closest point on the triangle.
 * @param[out] closesty The y-coordinate of the closest point on the triangle.
 * @param[out] closestz The z-coordinate of the closest point on the triangle.
 * @param[out] status The status indicating the relative position of the closest point.
 *
 * @return The unsigned distance between the point and the triangle.
 */
double basic_calculation_3d::f_tria_seg(double x0, double y0, double z0, 
        double xA, double yA, double zA, double xB, double yB, double zB, double xC, double yC, double zC,
        double a, double b, double c, double &closestx, double &closesty, double &closestz, int &status){
    
    double d; 
    double xx=x0-xA;
    double yy=y0-yA;
    double zz=z0-zA;
    d=a*xx+b*yy+c*zz;
    //projected point on plane, Pp, of P0
    double Ppx=x0-d*a; 
    double Ppy=y0-d*b; 
    double Ppz=z0-d*c; 
    //find Pp on plane's relative location to the triangle ABC
    //u=VB-VA
    double ux=xB-xA; double uy=yB-yA; double uz=zB-zA;
    //v=VC-VA
    double vx=xC-xA; double vy=yC-yA; double vz=zC-zA;
    //w=Pp-VA
    double wx=Ppx-xA; double wy=Ppy-yA; double wz=Ppz-zA;
    //uv=u.v; wu=w.u; wv=w.v; uu=u.u; vv=v.v
    double uv=ux*vx+uy*vy+uz*vz;
    double wu=wx*ux+wy*uy+wz*uz;
    double wv=wx*vx+wy*vy+wz*vz;
    double uu=ux*ux+uy*uy+uz*uz;
    double vv=vx*vx+vy*vy+vz*vz;
    //den=uv*uv-uu*vv
    double den=uv*uv-uu*vv;
    //s=(uv*wv-vv*wu)/den; t=(uv*wu-uu*wv)/den
    double s=(uv*wv-vv*wu)/den;
    double t=(uv*wu-uu*wv)/den;
    
    //in/on tria:
    if(s>=0-(1e-7) && t>=0-(1e-7) && s+t<=1+(1e-7)){
        //get distance
        d=abs(d);  //to get absolute distance, use abs
        //get status
        if(s>0-(1e-7) && t>0-(1e-7) && s+t<1+(1e-7)){ //interior
            status=6;
        }
        else if(s>=-(1e-7) && s<=(1e-7)){ //s=0
            if(t>=-(1e-7) && t<=(1e-7)) {status=0;} //t=0
            else if(t>=1-(1e-7) && t<=1+(1e-7)) {status=2;}
            else {status=5;}
        }
        else if(t>=-(1e-7) && t<=(1e-7)){ //t=0
            if(s>=1-(1e-7) && s<=1+(1e-7)) {status=1;}
            else {status=3;}
        }
        else {status=4;}
        //get closest pt info
        closestx=Ppx;
        closesty=Ppy;
        closestz=Ppz;
        return d;
    }
    //outside tria
    else{
        double x_temp[3]; double y_temp[3]; double z_temp[3];
        int status_temp[3]; //status_temp is 0,1,2; 0: closest point on line segment is lp1, 1: in between, 2: is lp2
        double dis_temp[3];
        //calculate distance and closest point to the three edges of the triangle
        p_lineSeg_dis(x0, y0, z0, xA, yA, zA, xB, yB, zB, x_temp[0], y_temp[0], z_temp[0], status_temp[0], dis_temp[0]);
        p_lineSeg_dis(x0, y0, z0, xB, yB, zB, xC, yC, zC, x_temp[1], y_temp[1], z_temp[1], status_temp[1], dis_temp[1]);
        p_lineSeg_dis(x0, y0, z0, xC, yC, zC, xA, yA, zA, x_temp[2], y_temp[2], z_temp[2], status_temp[2], dis_temp[2]);
        
        d=dis_temp[0]; closestx=x_temp[0]; closesty=y_temp[0]; closestz=z_temp[0]; 
        if(status_temp[0]==0){status=0;}
        else if(status_temp[0]==1){status=3;}
        else{status=1;}
        
        for(int k=1; k<3; k++){
            if(dis_temp[k]<d-1e-7){
                d=dis_temp[k]; closestx=x_temp[k]; closesty=y_temp[k]; closestz=z_temp[k]; 
                if(k==1){
                    if(status_temp[0]==0){status=1;}
                    else if(status_temp[0]==1){status=4;}
                    else{status=2;}
                }
                if(k==2){
                    if(status_temp[0]==0){status=2;}
                    else if(status_temp[0]==1){status=5;}
                    else{status=0;}
                }
            }
        }      
        return d;
    }
}