#include "parallel_meshing_2d.hh"

#include <sys/types.h>
#include <sys/stat.h>

#include <iostream>
#include <fstream>

using namespace voro;

/**
 * @brief A class representing a user-defined signed distance function in 2D. 
 * Here, we construct a SDF for a hexagram star shape. 
 * The code of the SDF is adapted from https://iquilezles.org/articles/distfunctions2d/. 
 *
 * This class derives from the base class `shape_2d` and provides user-defined implementations for
 * computing signed distance value at a given point.
 */
class shape_2d_mySDF : public shape_2d{
    public: 
        /**
         * Constructor for shape_2d_mySDF, a shape given by user-defined signed distance field.
         *
         * @param con_ The container_2d object.
         * @param num_t_ Number of parallel threads.
         * @param centerX_ Center of the hexagram shape, x-coordinate.
         * @param centerY_ Center of the hexagram shape, y-coordinate.
         * @param radius_ Radius of the circumscribed circle.
         */
         shape_2d_mySDF(container_2d &con_, int num_t_,double centerX_, double centerY_, double radius_)
            :shape_2d(con_,num_t_), centerX(centerX_), centerY(centerY_), radius(radius_)
             {
                get_geometryGrid();
             }
             
        /**
          * Destructor for shape_2d_mySDF.
          * Frees any resources allocated by the object.
          */
         ~shape_2d_mySDF(){};

        double centerX; //Center of the hexagram shape, x-coordinate.
        double centerY; //Center of the hexagram shape, y-coordinate.
        double radius; //Radius of the circumscribed circle.

        

         /**
          * Signed distance function calculation of a hexagram star shape. 
          * Code adapted from https://iquilezles.org/articles/distfunctions2d/. 
          * 
          * @param x The x-coordinate.
          * @param y The y-coordinate.
          * @return The signed distance of the point (x, y) to the shape.
          *         Negative values indicate inside, positive values indicate outside.
          */
         double sdf(double x, double y){
            double kx=-0.5;
            double ky=0.86602540378;
            double kz=0.57735026919;
            double kw=1.73205080757;
            //Translate (x,y) to be centered around the origin, since the hexagram SDF code is centered at the origin. 
            x-=centerX; 
            y-=centerY; 
            double px=abs(x);
            double py=abs(y);
            double fac1=min((kx*px+ky*py),0.0)*2.0;
            double fac2=min((ky*px+kx*py),0.0)*2.0;
            px-=fac1*kx;
            py-=fac1*ky;
            px-=fac2*ky;
            py-=fac2*kx;
            double fac3=min(max(px,radius*kz),radius*kw);
            px-=fac3;
            py-=radius;
            double d=sqrt(px*px+py*py);
            if(py<0.0){
                d=-d;
            }
            return d;

         };

};

int main() {
    
//------------------0.Set up parameters--------------------
    //Specify number of parallel threads
    int num_t=4;
    //Method index, 0 Distmesh; 1 CVD; 2 Hybrid
    int method_ind=2; 
    //Number of vertices in the mesh
    int Ntotal=5000; 
    //0:uniform sizing field
    double K=0.1; 
    //File output frequency. 0 no output; -1 output at termination; -2 output initial and final; 10, every 10 triangulations output
    int output_interval=-1; 

    //File output name prefix
    char case_name_base[256];
    sprintf(case_name_base,"mySDF_mesh_N_%d_K_%g",Ntotal,K);
    //Create directory to store file outputs
    mkdir(case_name_base,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

//------------------1.Create Voro++ container_2d----------------------
    printf("create container\n");
    int cnx=sqrt(Ntotal/3.3); int cny=cnx;
    container_2d con(0.0,1.0,0.0,1.0,cnx,cny,false,false,16,num_t);


//-------------------2.Create shape------------------------
    printf("create shape\n");
    //User-defined SDF shape
    shape_2d_mySDF myShp(con,num_t,0.5,0.5,0.2);
    //Circular shape
    shape_2d_circle cir(con,num_t,0.1,0.5,0.5);
    //Boolean difference of the two shapes
    shape_2d_difference shp(con,num_t,&myShp,&cir);


//--------------------3.Create adaptive sizing field------------------
    printf("create size_field\n");
    sizing_2d_automatic size_field(&shp,K);

//------------------4.Create parallel_meshing_2d object----------------------
    printf("create pm2d and pt init\n");
    //Optional: set seed for rand() for generating randome points in pt_init
    srand(10);
    //Create the pm2d object
    parallel_meshing_2d pm2d(&con, &shp, &size_field, num_t, output_interval, case_name_base);
    //Initialize meshing points
    pm2d.pt_init(Ntotal);

//======================Parallel meshing=====================
    printf("meshing\n");
    if(method_ind==0){ //DistMesh
        mesh_alg_2d_dm mesh_method(&pm2d);
        pm2d.meshing(&mesh_method); 
    }
    else if(method_ind==1){ //CVD
        mesh_alg_2d_cvd mesh_method(&pm2d);
        pm2d.meshing(&mesh_method); 
    }
    else if(method_ind==2){ //Hybrid
        mesh_alg_2d_hybrid mesh_method(&pm2d);
        pm2d.meshing(&mesh_method); 
    }

    return 0;   
}



