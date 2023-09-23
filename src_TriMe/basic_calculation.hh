/**
 * @file basic_calculation.hh
 * @brief Contains the basic_calculation namespace and classes for performing mathematical and geometric calculations.
 */

#ifndef BASIC_CALCULATION_HH
#define BASIC_CALCULATION_HH


#include "config.hh"

#include <stdlib.h>
#include <cstdlib>
#include <vector>
#include <limits>
#include <algorithm>
#include <fstream>
#include <cstdio>
#include <math.h> 
#include <stdexcept>
#include <string>
#include <cstring>
#include <iostream>
#include <map>
#include <ctime>
#include <omp.h>


/**
 * @namespace voro
 * @brief Namespace for the voro library.
 */

namespace voro {

	/**
	 * @class basic_calculation
	 * @brief Provides basic mathematical functions.
	 */
	class basic_calculation{
	public: 
		/**
		 * @brief Default constructor.
		 */
		basic_calculation(){};

		/**
		 * @brief Default destructor.
		 */
		~basic_calculation(){};

		/**
		 * @brief Calculates the square of a number.
		 * @param x The input number.
		 * @return The square of the input number.
		 */
	    inline double sqr(double x){return x*x;}

	    /**
		 * @brief Generates a random number between 0 and 1. Use in serial code.
		 * @return The random number.
		 */
	    inline double rnd() {return (1./RAND_MAX)*static_cast<double>(rand());}
	    
	    /**
		 * @brief Generates a random number between 0 and 1. Can be used in parallel code.
		 * @param seedptr Pointer to the seed value for the random number generator.
		 * @return The random number.
		 */
	    inline double rnd_r(unsigned int *seedptr){return (1./RAND_MAX)*static_cast<double>(rand_r(seedptr));}
	};

	/**
	 * @class basic_calculation_2d
	 * @brief Provides basic 2D geometric calculations.
	 */
	class basic_calculation_2d : public basic_calculation{
	public: 
		/**
		 * @brief Default constructor.
		 */
		basic_calculation_2d():basic_calculation(){};

		/**
		 * @brief Default destructor.
		 */
		~basic_calculation_2d(){};

		/**
		 * @brief Calculates the distance between two points in 2D space.
		 * @param x1 X-coordinate of the first point.
		 * @param y1 Y-coordinate of the first point.
		 * @param x2 X-coordinate of the second point.
		 * @param y2 Y-coordinate of the second point.
		 * @return The distance between the two points.
		 */
	    double d_points(double x1, double y1, double x2, double y2){
	        return sqrt(sqr(x1-x2)+sqr(y1-y2));
	    }

	    /**
		 * @brief Calculates the length of an edge in 2D space.
		 * @param x1 X-coordinate of the first endpoint of the edge.
		 * @param y1 Y-coordinate of the first endpoint of the edge.
		 * @param x2 X-coordinate of the second endpoint of the edge.
		 * @param y2 Y-coordinate of the second endpoint of the edge.
		 * @return The length of the edge.
		 */
	    double edge_len(double x1, double y1, double x2, double y2){
		    double d1=x1-x2; double d2=y1-y2;
		    double len=sqrt(d1*d1+d2*d2);
		    return len;
		}

		//triangle related
		/**
		 * @brief Calculates the semiperimeter of a triangle.
		 * @param a Length of the first edge.
		 * @param b Length of the second edge.
		 * @param c Length of the third edge.
		 * @return The semiperimeter of the triangle.
		 */
		double s_tria(double a, double b, double c){
		    return 0.5*(a+b+c);
		}

		/**
		 * @brief Calculates the area of a triangle using Heron's formula.
		 * @param a Length of the first edge.
		 * @param b Length of the second edge.
		 * @param c Length of the third edge.
		 * @param s Semiperimeter of the triangle.
		 * @return The area of the triangle.
		 */
		double area_tria(double a, double b, double c, double s){
		    return sqrt(s*(s-a)*(s-b)*(s-c));
		}

		/**
		 * @brief Calculates the aspect ratio of a triangle.
		 * @param a Length of the first edge.
		 * @param b Length of the second edge.
		 * @param c Length of the third edge.
		 * @param s Semiperimeter of the triangle.
		 * @return The aspect ratio of the triangle.
		 */
		double aspect_ratio_tria(double a, double b, double c, double s){
		    return a*b*c/(8*(s-a)*(s-b)*(s-c));
		}

		/**
		 * @brief Calculates the circumradius of a triangle using its side lengths and area.
		 * @param a Length of the first edge.
		 * @param b Length of the second edge.
		 * @param c Length of the third edge.
		 * @param area Area of the triangle.
		 * @return The circumradius of the triangle.
		 */
		double circumradius_tria(double a, double b, double c, double area){
			return 0.25*a*b*c/area;
		}

		/**
		 * @brief Calculates the edge ratio of a triangle using its side lengths.
		 * @param a Length of the first edge.
		 * @param b Length of the second edge.
		 * @param c Length of the third edge.
		 * @return The edge ratio of the triangle.
		 */
		double edge_ratio_tria(double a, double b, double c){
		    double lmin=a;
		    double lmax=a;
		    if(b<lmin){lmin=b;}
		    if(c<lmin){lmin=c;}
		    if(b>lmax){lmax=b;}
		    if(c>lmax){lmax=c;}
			if(lmax/lmin<1.0){printf("a,b,c: %g %g %g er: %g \n",a,b,c,lmax/lmin);};
			return lmax/lmin;
		}       

		/**
		 * @brief Calculates the circumradius of a triangle using its vertex coordinates.
		 * @param x0 X-coordinate of the first vertex.
		 * @param y0 Y-coordinate of the first vertex.
		 * @param x1 X-coordinate of the second vertex.
		 * @param y1 Y-coordinate of the second vertex.
		 * @param x2 X-coordinate of the third vertex.
		 * @param y2 Y-coordinate of the third vertex.
		 * @return The circumradius of the triangle.
		 */
		double circumradius_tria(double x0,double y0,double x1,double y1,double x2,double y2){
			double a=d_points(x0,y0,x1,y1);
			double b=d_points(x1,y1,x2,y2);
			double c=d_points(x0,y0,x2,y2);
			double s=s_tria(a,b,c);
			double area=area_tria(a,b,c,s);
			return circumradius_tria(a,b,c,area);
		}

		/**
		 * @brief Calculates the centroid of a triangle using its vertex coordinates.
		 * @param x0 X-coordinate of the first vertex.
		 * @param y0 Y-coordinate of the first vertex.
		 * @param x1 X-coordinate of the second vertex.
		 * @param y1 Y-coordinate of the second vertex.
		 * @param x2 X-coordinate of the third vertex.
		 * @param y2 Y-coordinate of the third vertex.
		 * @param cx Reference to the variable to store the X-coordinate of the centroid.
		 * @param cy Reference to the variable to store the Y-coordinate of the centroid.
		 */
		void centroid_tria(double x0,double y0,double x1,double y1,double x2,double y2,double &cx,double &cy){
			double one_third=1.0/3.0;
			cx=(x0+x1+x2)*one_third;
			cy=(y0+y1+y2)*one_third;
		}

		/**
		 * @brief Calculates the circumcenter of a triangle using its vertex coordinates.
		 * @param x1 X-coordinate of the first vertex.
		 * @param y1 Y-coordinate of the first vertex.
		 * @param x2 X-coordinate of the second vertex.
		 * @param y2 Y-coordinate of the second vertex.
		 * @param x3 X-coordinate of the third vertex.
		 * @param y3 Y-coordinate of the third vertex.
		 * @param ccx Reference to the variable to store the X-coordinate of the circumcenter.
		 * @param ccy Reference to the variable to store the Y-coordinate of the circumcenter.
		 */
		void circumcenter_tria(double x1,double y1,double x2,double y2,double x3,double y3,double &ccx,double &ccy){
			//reference: https://mathworld.wolfram.com/Circumcircle.html
			double a=x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2);
			double bx=x1*x1*(y3-y2) + x2*x2*(y1-y3) + x3*x3*(y2-y1) +
						y1*y1*(y3-y2) + y2*y2*(y1-y3) + y3*y3*(y2-y1);
			double by=y1*y1*(x2-x3) + y2*y2*(x3-x1) + y3*y3*(x1-x2) +
						x1*x1*(x2-x3) + x2*x2*(x3-x1) + x3*x3*(x1-x2);
			double fac_temp= -1.0/(2*a);
			ccx=bx*fac_temp;
			ccy=by*fac_temp;
		}

		/**
		 * @brief Calculates the aspect ratio of a triangle using its vertex coordinates.
		 * @param x0 X-coordinate of the first vertex.
		 * @param y0 Y-coordinate of the first vertex.
		 * @param x1 X-coordinate of the second vertex.
		 * @param y1 Y-coordinate of the second vertex.
		 * @param x2 X-coordinate of the third vertex.
		 * @param y2 Y-coordinate of the third vertex.
		 * @return The aspect ratio of the triangle.
		 */
		double aspect_ratio_tria(double x0,double y0,double x1,double y1,double x2,double y2){
		    double a=d_points(x0,y0,x1,y1);
			double b=d_points(x1,y1,x2,y2);
			double c=d_points(x0,y0,x2,y2);
			double s=s_tria(a,b,c);
		    return aspect_ratio_tria(a,b,c,s);
		}

		/**
		 * @brief Calculates the edge ratio of a triangle using its vertex coordinates.
		 * @param x0 X-coordinate of the first vertex.
		 * @param y0 Y-coordinate of the first vertex.
		 * @param x1 X-coordinate of the second vertex.
		 * @param y1 Y-coordinate of the second vertex.
		 * @param x2 X-coordinate of the third vertex.
		 * @param y2 Y-coordinate of the third vertex.
		 * @return The edge ratio of the triangle.
		 */
		double edge_ratio_tria(double x0,double y0,double x1,double y1,double x2,double y2){
		    double a=d_points(x0,y0,x1,y1);
			double b=d_points(x1,y1,x2,y2);
			double c=d_points(x0,y0,x2,y2);
		    return edge_ratio_tria(a,b,c);
		}       


		/**
		 * @brief Determines whether three points are in counterclockwise order.
		 * @param Ax X-coordinate of the first point.
		 * @param Ay Y-coordinate of the first point.
		 * @param Bx X-coordinate of the second point.
		 * @param By Y-coordinate of the second point.
		 * @param Cx X-coordinate of the third point.
		 * @param Cy Y-coordinate of the third point.
		 * @return True if the points are in counterclockwise order, False otherwise.
		 */
		bool ccw(double Ax, double Ay, double Bx, double By, double Cx, double Cy){
		    bool temp=false;
		    if((Cy-Ay) * (Bx-Ax) > (By-Ay) * (Cx-Ax)){
		        temp=true;
		    }
		     return temp;
		}

		/**
		 * @brief Checks if two line segments AB and CD intersect with each other.
		 * @param Ax X-coordinate of the first endpoint of line segment AB.
		 * @param Ay Y-coordinate of the first endpoint of line segment AB.
		 * @param Bx X-coordinate of the second endpoint of line segment AB.
		 * @param By Y-coordinate of the second endpoint of line segment AB.
		 * @param Cx X-coordinate of the first endpoint of line segment CD.
		 * @param Cy Y-coordinate of the first endpoint of line segment CD.
		 * @param Dx X-coordinate of the second endpoint of line segment CD.
		 * @param Dy Y-coordinate of the second endpoint of line segment CD.
		 * @return True if the line segments intersect, False otherwise.
		 */
		bool line_seg_intersect(double Ax, double Ay, double Bx, double By, 
		        double Cx, double Cy, double Dx, double Dy){
		    bool temp=false;
		    if((ccw(Ax,Ay,Cx,Cy,Dx,Dy) != ccw(Bx,By,Cx,Cy,Dx,Dy)) && (ccw(Ax,Ay,Bx,By,Cx,Cy) != ccw(Ax,Ay,Bx,By,Dx,Dy))){
		        temp=true;
		    }
		    return temp;
		}

		/**
		 * @brief Calculates the unsigned distance from a point (x, y) to a line segment defined by two points (x0, y0) and (x1, y1).
		 * @param x X-coordinate of the point.
		 * @param y Y-coordinate of the point.
		 * @param x0 X-coordinate of the first endpoint of the line segment.
		 * @param y0 Y-coordinate of the first endpoint of the line segment.
		 * @param x1 X-coordinate of the second endpoint of the line segment.
		 * @param y1 Y-coordinate of the second endpoint of the line segment.
		 * @param closestx Reference to the variable to store the X-coordinate of the closest point on the line segment.
		 * @param closesty Reference to the variable to store the Y-coordinate of the closest point on the line segment.
		 * @param status Reference to the variable to indicate the status of the point with respect to the line segment.
		 * @return The unsigned distance from the point to the line segment.
		 */
        double f_line_seg(double x, double y, double x0, double y0, double x1, double y1, double &closestx, double &closesty, int &status);
	
		
		/**
		 * @brief Projects a point (xx1, yy1) onto the geometry boundary using the provided functions.
		 * @param new_x Reference to the variable to store the projected X-coordinate.
		 * @param new_y Reference to the variable to store the projected Y-coordinate.
		 * @param xx1 X-coordinate of the point to be projected.
		 * @param yy1 Y-coordinate of the point to be projected.
		 * @param deps_prime The value of deps_prime, default deps value to use if no valid deps available.
		 * @param geps_prime The value of geps_prime, default geps value to use if no valid geps available.
		 * @param project_pt_outside_geo Flag indicating whether to project only point outside the geometry boundary.
		 * @return True if the point is successfully projected, False otherwise.
		 */
        bool pt_projection(double &new_x, double &new_y,
         double xx1, double yy1, 
         double deps_prime, double geps_prime,
         bool project_pt_outside_geo=true);

		/**
		 * @brief Virtual function to calculate the finite difference spatial stepsize (deps) value for a point (x, y) on the geometry boundary.
		 * @param x X-coordinate of the point.
		 * @param y Y-coordinate of the point.
		 * @param deps_prime The value of deps_prime.
		 * @return The calculated deps value.
		 */
		virtual double deps_func(double x, double y, double deps_prime){return 0.0;};

		/**
		 * @brief Virtual function to calculate the geometry boundary tolerence (geps) value for a point (x, y) on the geometry boundary.
		 * @param x X-coordinate of the point.
		 * @param y Y-coordinate of the point.
		 * @param deps_prime The value of deps_prime.
		 * @return The calculated geps value.
		 */
		virtual double geps_func(double x, double y, double deps_prime){return 0.0;};

		/**
		 * @brief Virtual function to calculate the signed distance function (SDF) for a point (x, y).
		 * @param x X-coordinate of the point.
		 * @param y Y-coordinate of the point.
		 * @return The calculated SDF value.
		 */
		virtual double sdf_func(double x, double y){return 0.0;};

		/**
		 * @brief Virtual function to check if a point (x, y) is inside the outer grid of the geometry.
		 * @param x X-coordinate of the point.
		 * @param y Y-coordinate of the point.
		 * @return True if the point is inside the outer grid, False otherwise.
		 */
		virtual bool pt_in_outer_grid(double x,double y){return false;};


	};

	/**
	 * @class basic_calculation_3d
	 * @brief Provides basic 3D geometric calculations.
	 */
	class basic_calculation_3d : public basic_calculation{
	public: 
		/**
		 * @brief Default constructor.
		 */
		basic_calculation_3d():basic_calculation(){};

		/**
		 * @brief Default destructor.
		 */
		~basic_calculation_3d(){};

		/**
		 * @brief Calculates the distance between two points in 3D space.
		 * @param x1 X-coordinate of the first point.
		 * @param y1 Y-coordinate of the first point.
		 * @param z1 Z-coordinate of the first point.
		 * @param x2 X-coordinate of the second point.
		 * @param y2 Y-coordinate of the second point.
		 * @param z2 Z-coordinate of the second point.
		 * @return The distance between the two points.
		 */
	    double d_points(double x1, double y1, double z1, double x2, double y2, double z2){
            return sqrt(sqr(x1-x2)+sqr(y1-y2)+sqr(z1-z2));
        }

        /**
		 * @brief Calculates the length of an edge in 3D space.
		 * @param x1 X-coordinate of the first endpoint of the edge.
		 * @param y1 Y-coordinate of the first endpoint of the edge.
		 * @param z1 Z-coordinate of the first endpoint of the edge.
		 * @param x2 X-coordinate of the second endpoint of the edge.
		 * @param y2 Y-coordinate of the second endpoint of the edge.
		 * @param z2 Z-coordinate of the second endpoint of the edge.
		 * @return The length of the edge.
		 */
	    double edge_len(double x1, double y1, double z1, double x2, double y2, double z2){
		    double d1=x1-x2; double d2=y1-y2; double d3=z1-z2;
		    double len=sqrt(d1*d1+d2*d2+d3*d3);
		    return len;
		}

		/**
         * @brief Calculates the distance and closest point information of a point to a line segment.
         * @param x0 X-coordinate of the point.
         * @param y0 Y-coordinate of the point.
         * @param z0 Z-coordinate of the point.
         * @param lx1 X-coordinate of the first endpoint of the line segment.
         * @param ly1 Y-coordinate of the first endpoint of the line segment.
         * @param lz1 Z-coordinate of the first endpoint of the line segment.
         * @param lx2 X-coordinate of the second endpoint of the line segment.
         * @param ly2 Y-coordinate of the second endpoint of the line segment.
         * @param lz2 Z-coordinate of the second endpoint of the line segment.
         * @param x_temp Reference to the variable to store the X-coordinate of the closest point.
         * @param y_temp Reference to the variable to store the Y-coordinate of the closest point.
         * @param z_temp Reference to the variable to store the Z-coordinate of the closest point.
         * @param status_temp Reference to the variable to indicate the status of the point with respect to the line segment.
         * @param dis_temp Reference to the variable to store the distance from the point to the line segment.
         */
		void p_lineSeg_dis(double x0, double y0, double z0, double lx1, double ly1, double lz1, double lx2, double ly2, double lz2, 
            double &x_temp, double &y_temp, double &z_temp, int &status_temp, double &dis_temp);

        /**
         * @brief Calculates the unsigned distance from a point (x0, y0, z0) to a triangle (VA, VB, VC).
         * @param x0 X-coordinate of the point.
         * @param y0 Y-coordinate of the point.
         * @param z0 Z-coordinate of the point.
         * @param xA X-coordinate of the first vertex of the triangle.
         * @param yA Y-coordinate of the first vertex of the triangle.
         * @param zA Z-coordinate of the first vertex of the triangle.
         * @param xB X-coordinate of the second vertex of the triangle.
         * @param yB Y-coordinate of the second vertex of the triangle.
         * @param zB Z-coordinate of the second vertex of the triangle.
         * @param xC X-coordinate of the third vertex of the triangle.
         * @param yC Y-coordinate of the third vertex of the triangle.
         * @param zC Z-coordinate of the third vertex of the triangle.
         * @param a Length of the triangle edge opposite to vertex A.
         * @param b Length of the triangle edge opposite to vertex B.
         * @param c Length of the triangle edge opposite to vertex C.
         * @param closestx Reference to the variable to store the X-coordinate of the closest point on the triangle.
         * @param closesty Reference to the variable to store the Y-coordinate of the closest point on the triangle.
         * @param closestz Reference to the variable to store the Z-coordinate of the closest point on the triangle.
         * @param status Reference to the variable to indicate the status of the point with respect to the triangle.
         * @return The unsigned distance from the point to the triangle.
         */
        double f_tria_seg(double x0, double y0, double z0, double xA, double yA, double zA, double xB, double yB, double zB, 
            double xC, double yC, double zC, double a, double b, double c, 
            double &closestx, double &closesty, double &closestz, int &status);


	};

}

#endif