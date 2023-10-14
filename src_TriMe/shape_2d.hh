#ifndef SHAPE_2D_HH
#define SHAPE_2D_HH

#include "custom_shape_2d.hh"
#include "basic_calculation.hh"

#include <algorithm>
#include <limits>
#include <vector>
#include "voro++.hh"
#include <fstream>

namespace voro {

    /**
     * @brief A 2D shape class that represents a geometric shape and provides various operations and calculations.
     */
    class shape_2d : public basic_calculation_2d {
    public:

        /**
         * @brief Constructor for the shape_2d class.
         * 
         * @param con_ A reference to a container_2d object.
         * @param num_t_ The number of parallel threads.
         */
        shape_2d(container_2d &con_, int num_t_);

        /**
         * @brief Destructor for the shape_2d class.
         */
        ~shape_2d();
        int num_t; /**< The number of parallel threads. */
        double ax, bx, ay, by; /**< The domain size. */

        /**
         * Signed distance function: outside +, inside -.
         * For the actual implementation of the custom shape sdf(x,y),
         * it should be constructed based on the geometry grid grometry_ij.
         * @param x X-coordinate of the point.
         * @param y Y-coordinate of the point.
         * @return The signed distance from the point to the shape.
         */
        virtual double sdf(double x, double y)=0;

        //Geometry grid
        int gnx, gny, gnxy; /**< Size of the geometry grid. */
        double gdx, gdy, diag_gdxy, inv_gdx, inv_gdy; /**< Geometry grid parameters. */
        int *geo_grid; /**< Underlying geometry grid. -1,-2,-3: inside; 1,2,3: boundary; gnxy+1,gnxy+2,...: outside. */
        int *geo_igrid_ij; /**< Inner grid ij's array. */
        int *geo_bgrid_ij; /**< Boundary grid ij's array. */
        int *geo_ogrid_ij; /**< Outer grid ij's array. */
        int geo_bgrid_ct; /**< Number of geometry boundary grids. */
        int geo_igrid_ct; /**< Number of geometry inner grids. */
        int geo_ogrid_ct; /**< Number of geometry outer grids. */

        //shape scaling parameters
        bool shape_scaling; /**< Flag to determine scaling shape input or not. Only available for custom_shape_2d input. */
        double scale_min_domain_range; /**< scaling parameter, std::min(bx-ax,by-ay). */
        double scale_xmid; /**< scaling parameter, midpoint x-coordinate of the input shape. */
        double scale_ymid; /**< scaling parameter, midpoint y-coordinate of the input shape. */
        double scale_max_range; /**< scaling parameter, max(height,width) of the input shape. */

        /**
         * @brief Gets the geometry grid of size gnx*gny categorizing inner/outer/boundary grid.
         * 
         * Useful for efficient outer point projection, fine feature resolution, and adf construction only on bdry grids.
         */
        void get_geometryGrid();
        
        /**
         * Helper function to get the value of the geometry grid at a given point (x, y).
         * @param x X-coordinate of the point.
         * @param y Y-coordinate of the point.
         * @return The value of the geometry grid at the given point.
         */
        int pt_geo_grid_val(double x, double y);
        
         /**
         * Print the fields (geo_grid, sdf) to a file.
         * @param case_name The name of the case.
         */
        void print_geo_grid_sdf_to_file(const char *case_name);

        /**
         * Helper variables and functions for using custom_shape_2d.
         */

        bool is_custom_shape_contour; /**< Helper flag for constructing shape from custom shape contour line segment inputs. */
        
        /**
         * Get the number of boundaries.
         * @return The number of boundaries.
         */
        virtual int b_ct(){return 0;} 

        /**
         * Get the number of line segments on a specific boundary.
         * @param bi Boundary index.
         * @return The number of line segments on the specified boundary.
         */
        virtual int seg_ct(int bi){return 0;} 

        /**
         * Get the x-coordinate of the starting point of a line segment on a specific boundary.
         * @param bi Boundary index.
         * @param li Line segment index.
         * @return The x-coordinate of the starting point of the specified line segment.
         */
        virtual double b_pts_x(int bi, int li){return 0.0;}
        
        /**
         * Get the y-coordinate of the starting point of a line segment on a specific boundary.
         * @param bi Boundary index.
         * @param li Line segment index.
         * @return The y-coordinate of the starting point of the specified line segment.
         */
        virtual double b_pts_y(int bi, int li){return 0.0;}
        
    };
    
    
//Derived Classes for Shape helper functions: 
//return signed distance of (x, y) to boundary

    
    /**
     * Class representing the union of two shapes, where the signed distance is defined as min(dA, dB).
     */
    class shape_2d_union : public shape_2d{
     public:
         shape_2d *shp1; /**< Pointer to the first shape. */
         shape_2d *shp2; /**< Pointer to the second shape. */
         
         /**
         * Constructor for shape_2d_union, a new shape that's the union of two shapes.
         * @param con_ The container_2d object.
         * @param num_t_ Number of parallel threads.
         * @param shp1_ Pointer to the first shape.
         * @param shp2_ Pointer to the second shape.
         */
         shape_2d_union(container_2d &con_, int num_t_, shape_2d *shp1_, shape_2d *shp2_)
             :shape_2d(con_, num_t_),
             shp1(shp1_), shp2(shp2_) 
             {
                get_geometryGrid();
             }

        /** Destructor for shape_2d_union. */
         ~shape_2d_union(){};
         
         
         /**
         * Signed distance function for the union of two shapes.
         * @param x The x-coordinate.
         * @param y The y-coordinate.
         * @return The signed distance: outside +, inside -
         */
         double sdf(double x, double y){
             double d1=shp1->sdf(x, y);
             double d2=shp2->sdf(x, y);
             return d1<d2?d1:d2;
         }
    };
    
    
    /**
     * Class representing the difference of two shapes, where the signed distance is defined as max(dA, -dB).
     */
    class shape_2d_difference : public shape_2d{
     public:
         shape_2d *shp1; /**< Pointer to the first shape. */
         shape_2d *shp2; /**< Pointer to the second shape. */

        /**
         * Constructor for shape_2d_difference, a new shape that's the difference of two shapes.
         *
         * @param con_ The container_2d object.
         * @param num_t_ Number of parallel threads.
         * @param shp1_ Pointer to the first shape.
         * @param shp2_ Pointer to the second shape.
         */
         shape_2d_difference(container_2d &con_, int num_t_, shape_2d *shp1_, shape_2d *shp2_)
             :shape_2d(con_,num_t_),
             shp1(shp1_), shp2(shp2_) 
             {
                get_geometryGrid();
             }

        /** Destructor for shape_2d_difference. */
         ~shape_2d_difference(){};

         /**
         * Virtual function to compute the signed distance from a point (x, y) to the difference of the two shapes.
         * The signed distance is defined as max(dA, -dB), where dA is the signed distance to the first shape
         * and dB is the signed distance to the second shape.
         *
         * @param x The x-coordinate of the point.
         * @param y The y-coordinate of the point.
         * @return The signed distance (outside +, inside -) from the point to the difference of the two shapes.
         */
         double sdf(double x, double y){
             double d1=shp1->sdf(x, y);
             double d2=shp2->sdf(x, y);
             return d1>(-d2)?d1:(-d2);
         }
    };
    
    
    /**
     * Class representing the intersection of two shapes, where the signed distance is defined as max(dA, dB).
     */
    class shape_2d_intersection : public shape_2d{
     public:
         shape_2d *shp1; /**< Pointer to the first shape. */
         shape_2d *shp2; /**< Pointer to the second shape. */

        /**
         * Constructor for shape_2d_intersection, a new shape that's the intersection of two shapes.
         *
         * @param con_ The container_2d object.
         * @param num_t_ Number of parallel threads.
         * @param shp1_ Pointer to the first shape.
         * @param shp2_ Pointer to the second shape.
         */
         shape_2d_intersection(container_2d &con_, int num_t_, shape_2d *shp1_, shape_2d *shp2_)
             :shape_2d(con_,num_t_),
             shp1(shp1_), shp2(shp2_) 
             {
                get_geometryGrid();
             }

        /** Destructor for shape_2d_intersection. */
         ~shape_2d_intersection(){};
         
         /**
         * Virtual function to compute the signed distance from a point (x, y) to the intersection of the two shapes.
         * The signed distance is defined as max(dA, dB), where dA is the signed distance to the first shape
         * and dB is the signed distance to the second shape.
         *
         * @param x The x-coordinate of the point.
         * @param y The y-coordinate of the point.
         * @return The signed distance from the point to the intersection of the two shapes.
         */
         double sdf(double x, double y){
             double d1=shp1->sdf(x, y);
             double d2=shp2->sdf(x, y);
             return d1>d2?d1:d2;
         }
    };
    

    /**
     * @brief The shape_2d_rectangle class represents a 2D rectangle shape.
     *
     * This class inherits from the base class shape_2d.
     */
    class shape_2d_rectangle : public shape_2d {
     public:
         /**
          * Rectangle with sides defined by [x0, x1], [y0, y1].
          */
          const double x0; /**< The x-coordinate of the left side of the rectangle. */
         const double x1; /**< The x-coordinate of the right side of the rectangle. */
         const double y0; /**< The y-coordinate of the bottom side of the rectangle. */
         const double y1; /**< The y-coordinate of the top side of the rectangle. */
         
             
         /**
          * Constructor for shape_2d_rectangle, representing a rectangle with sides defined by [x0, x1], [y0, y1].
          * @param con_ The container_2d object.
          * @param num_t_ Number of parallel threads.
          * @param x0_ The x-coordinate of the left side of the rectangle.
          * @param x1_ The x-coordinate of the right side of the rectangle.
          * @param y0_ The y-coordinate of the bottom side of the rectangle.
          * @param y1_ The y-coordinate of the top side of the rectangle.
          */
         shape_2d_rectangle(container_2d &con_, int num_t_, double x0_, double x1_, double y0_, double y1_)
            :shape_2d(con_,num_t_),
             x0(x0_), x1(x1_), y0(y0_), y1(y1_) 
             {
                get_geometryGrid();
             }

        /**
          * Destructor for shape_2d_rectangle.
          * Frees any resources allocated by the object.
          */
         ~shape_2d_rectangle(){};

         /**
          * Signed distance function for the rectangle.
          * @param x The x-coordinate.
          * @param y The y-coordinate.
          * @return The signed distance of the point (x, y) to the rectangle.
          *         Negative values indicate inside, positive values indicate outside.
          */
         double sdf(double x, double y){
             double d = 0.0;
            if (x>=x0 && x<=x1 && y>=y0 && y<=y1){
                d = -min(abs(x-x0), min(abs(x-x1), min(abs(y-y0), abs(y-y1))));
                return d;
            }
            else {
                if(x>=x0 && x<=x1) {
                    d = min(abs(y-y0), abs(y-y1));
                    return d;
                }
                else if(y>=y0 && y<=y1){
                    d = min(abs(x-x0), abs(x-x1));
                    return d;
                }
                else {
                    if(x<x0){
                        d = sqrt(min(sqr(x-x0)+sqr(y-y1), sqr(x-x0)+sqr(y-y0)));
                        return d;
                    }
                    else {
                        d = sqrt(min(sqr(x-x1)+sqr(y-y1), sqr(x-x1)+sqr(y-y0)));
                        return d;
                    }
                }
            }
         }
    };


    /**
     * @brief The shape_2d_circle class represents a 2D circle shape.
     *
     * This class inherits from the base class shape_2d and implements the signed distance function (SDF)
     * for the circle shape. The SDF returns the signed distance from a given point (x, y) to the circle,
     * with the sign indicating whether the point is inside or outside the circle.
     */
    class shape_2d_circle : public shape_2d {
     public:
         const double r; /**< The radius of the circle. */
         const double x0; /**< The x-coordinate of the center of the circle. */
         const double y0; /**< The y-coordinate of the center of the circle. */
         
         
         /**
          * @brief Constructor for the shape_2d_circle class.
          *
          * @param con_ The container_2d object.
          * @param num_t_ Number of parallel threads.
          * @param r_ The radius of the circle.
          * @param x0_ The x-coordinate of the center of the circle.
          * @param y0_ The y-coordinate of the center of the circle.
          */
         shape_2d_circle(container_2d &con_, int num_t_, double r_, double x0_, double y0_)
            :shape_2d(con_,num_t_),
             r(r_), x0(x0_), y0(y0_) 
             {
                get_geometryGrid();
             }

        /**
          * Destructor for shape_2d_circle.
          * Frees any resources allocated by the object.
          */
         ~shape_2d_circle(){};
         
         /**
          * @brief Calculates the signed distance from a point (x, y) to the circle.
          *
          * The signed distance is the distance from the point to the circle's boundary,
          * with the sign indicating whether the point is inside or outside the circle.
          *
          * @param x The x-coordinate of the point.
          * @param y The y-coordinate of the point.
          * @return The signed distance from the point to the circle.
          */
         double sdf(double x, double y){
            double d = sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0))-r;
            return d;
         }
    };
    
    
    
    /**
     * @brief The shape_2d_ellipse class represents a 2D ellipse shape.
     *
     * This class inherits from the base class shape_2d and implements the signed distance function (SDF)
     * for the ellipse shape. The SDF returns the signed distance from a given point (x, y) to the ellipse,
     * with the sign indicating whether the point is inside or outside the ellipse.
     */
    class shape_2d_ellipse : public shape_2d{
     public:
         const double x0; /**< The x-coordinate of the center of the ellipse. */
         const double y0; /**< The y-coordinate of the center of the ellipse. */
         const double a;  /**< The horizontal radius (semi-major axis) of the ellipse. */
         const double b;  /**< The vertical radius (semi-minor axis) of the ellipse. */

         
         /**
         * @brief Constructor for the shape_2d_ellipse class: (x-x0)^2/a^2 +(y-y0)^2/b^2 = 1.0.
         *
         * @param con_ The container_2d object.
         * @param num_t_ Number of parallel threads.
         * @param x0_ The x-coordinate of the center of the ellipse.
         * @param y0_ The y-coordinate of the center of the ellipse.
         * @param a_ The horizontal radius (semi-major axis) of the ellipse.
         * @param b_ The vertical radius (semi-minor axis) of the ellipse.
         */
         shape_2d_ellipse(container_2d &con_, int num_t_, double x0_, double y0_, double a_, double b_)
            :shape_2d(con_,num_t_),
             x0(x0_), y0(y0_), a(a_), b(b_) 
             {
                get_geometryGrid();
             }

        /**
          * Destructor for shape_2d_ellipse.
          * Frees any resources allocated by the object.
          */
         ~shape_2d_ellipse(){};

         /**
         * @brief Calculates the signed distance from a point (x, y) to the ellipse.
         *
         * The signed distance is the distance from the point to the ellipse's boundary,
         * with the sign indicating whether the point is inside or outside the ellipse.
         *
         * @param x The x-coordinate of the point.
         * @param y The y-coordinate of the point.
         * @return The signed distance from the point to the ellipse.
         */
         double sdf(double x, double y){
            double d = 1.0/sqr(a)*sqr(x-x0)+1.0/sqr(b)*sqr(y-y0)-1.0;
            return d;
         }
    };

    /**
     * @brief The shape_2d_superellipse class represents a 2D superellipse shape.
     *
     * This class inherits from the base class shape_2d and implements the signed distance function (SDF)
     * for the superellipse shape. The SDF returns the signed distance from a given point (x, y) to the superellipse,
     * with the sign indicating whether the point is inside or outside the superellipse.
     */
    class shape_2d_superellipse : public shape_2d{
     public: 
         const double r; /**< The radius of the superellipse. */
         const double x0; /**< The x-coordinate of the center of the superellipse. */
         const double y0; /**< The y-coordinate of the center of the superellipse. */

        /**
         * @brief Constructor for the shape_2d_superellipse class.
         *
         * @param con_ The container_2d object.
         * @param num_t_ Number of parallel threads.
         * @param r_ The radius of the superellipse.
         * @param x0_ The x-coordinate of the center of the superellipse.
         * @param y0_ The y-coordinate of the center of the superellipse.
         */
        shape_2d_superellipse(container_2d &con_, int num_t_, double r_, double x0_, double y0_)
        : shape_2d(con_,num_t_),
          r(r_),x0(x0_),y0(y0_) 
          {
            get_geometryGrid();
          }

        /**
          * Destructor for shape_2d_superellipse.
          * Frees any resources allocated by the object.
          */
        ~shape_2d_superellipse(){};

        /**
         * @brief Calculates the signed distance from a point (x, y) to the superellipse.
         *
         * The signed distance is the distance from the point to the superellipse's boundary,
         * with the sign indicating whether the point is inside or outside the superellipse.
         *
         * @param x The x-coordinate of the point.
         * @param y The y-coordinate of the point.
         * @return The signed distance from the point to the superellipse.
         */
         double sdf(double x, double y){
            double d = sqrt(sqrt(sqr(x-x0)*sqr(x-x0)+sqr(y-y0)*sqr(y-y0)))-r;
            return d;
         }
    };
    
    
    
    /**
     * Class representing a custom shape defined by user-input contour line segments, and generate an adaptive signed distance field based on it.
     */
    class shape_2d_contour_lines : public shape_2d {
     public:
            
         custom_shape_2d cshp; /**< The custom_shape_2d object representing the shape formed by contour line segments. */
         

         /**
         * Constructor for shape_2d_contour_lines, a custom shape defined by user-input contour line segments.
         *
         * @param con_ The container_2d object.
         * @param num_t_ Number of parallel threads.
         * @param num_t_cshp_ Number of parallel threads for the custom shape.
         * @param boundaries_ A vector of vectors representing the contour line segments of the shape.
         * @param normalize_model_ Flag indicating whether to normalize the model.
         */
         shape_2d_contour_lines(container_2d &con_, int num_t_, int num_t_cshp_, std::vector<std::vector<double>> boundaries_, bool normalize_model_)
            :shape_2d(con_,num_t_),
             cshp(custom_shape_2d(boundaries_, con_.ax, con_.bx, con_.ay, con_.by, normalize_model_,num_t_cshp_))
             {
                if(normalize_model_==true){
                    shape_scaling=true;
                    scale_min_domain_range=cshp.scale_min_domain_range;
                    scale_xmid=cshp.scale_xmid;
                    scale_ymid=cshp.scale_ymid;
                    scale_max_range=cshp.scale_max_range;
                }
                is_custom_shape_contour=true;
                get_geometryGrid();
             }
             
        /**
         * Destructor for shape_2d_contour_lines.
         */
         ~shape_2d_contour_lines(){};
         

         /**
         * Virtual function to compute the signed distance from a point (x, y) to the custom shape defined by the contour lines.
         *
         * @param x The x-coordinate of the point.
         * @param y The y-coordinate of the point.
         * @return The signed distance from the point to the custom shape.
         */
         double sdf(double x, double y){
            double d = cshp.f_b_pts(x,y);
            return d;
         }

         /**
         * Function to get the number of boundaries in the custom shape.
         *
         * @return The number of boundaries.
         */
        int b_ct(){return cshp.b_ct;} 

        /**
         * Function to get the number of line segments on a specific boundary of the custom shape.
         *
         * @param bi The index of the boundary.
         * @return The number of line segments on the boundary.
         */
        int seg_ct(int bi){return cshp.seg_ct[bi];} 

        /**
         * Function to get the x-component of the starting point of a line segment on a specific boundary of the custom shape.
         *
         * @param bi The index of the boundary.
         * @param li The index of the line segment on the boundary.
         * @return The x-component of the starting point of the line segment.
         */
        double b_pts_x(int bi, int li){return cshp.b_pts[bi][2*li];}

        /**
         * Function to get the y-component of the starting point of a line segment on a specific boundary of the custom shape.
         *
         * @param bi The index of the boundary.
         * @param li The index of the line segment on the boundary.
         * @return The y-component of the starting point of the line segment.
         */
        double b_pts_y(int bi, int li){return cshp.b_pts[bi][2*li+1];}
        


    };
    
    







}

#endif
