#include "parallel_meshing_2d.hh"

#include <sys/types.h>
#include <sys/stat.h>

#include <iostream>
#include <fstream>

using namespace voro;


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
    sprintf(case_name_base,"boolean_mesh_N_%d_K_%g",Ntotal,K);
    //Create directory to store file outputs
    mkdir(case_name_base,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

//------------------1.Create Voro++ container_2d----------------------
    printf("create container\n");
    int cnx=sqrt(Ntotal/3.3); int cny=cnx;
    container_2d con(0.0,1.0,0.0,1.0,cnx,cny,false,false,16,num_t);


//-------------------2.Create shape------------------------
    printf("create shape\n");

    //A. Intersection
    shape_2d_circle cdl(con,num_t,0.15,0.4,0.4);
    shape_2d_circle cdr(con,num_t,0.15,0.6,0.4);
    shape_2d_circle cul(con,num_t,0.15,0.4,0.6);
    shape_2d_circle cur(con,num_t,0.15,0.6,0.6);

    shape_2d_intersection leafd(con,num_t,&cdl,&cdr);
    shape_2d_intersection leafl(con,num_t,&cdl,&cul);
    shape_2d_intersection leafu(con,num_t,&cul,&cur);
    shape_2d_intersection leafr(con,num_t,&cdr,&cur);

    //B. Union
    shape_2d_circle cm(con,num_t,0.1,0.5,0.5);

    shape_2d_union m0(con,num_t,&leafd,&leafl);
    shape_2d_union m1(con,num_t,&m0,&leafu);
    shape_2d_union m2(con,num_t,&m1,&leafr);
    shape_2d_union m3(con,num_t,&m2,&cm);


    //C. Difference
    shape_2d_rectangle rec(con,num_t,0.1,0.9,0.1,0.9);

    shape_2d_difference shp_final(con,num_t,&rec,&m3);


//--------------------3.Create adaptive sizing field------------------
    printf("create size_field\n");
    sizing_2d_automatic size_field(&shp_final,K);

//------------------4.Create parallel_meshing_2d object----------------------
    printf("create pm2d and pt init\n");
    //Optional: set seed for rand() for generating randome points in pt_init
    srand(10);
    //Create the pm2d object
    parallel_meshing_2d pm2d(&con, &shp_final, &size_field, num_t, output_interval, case_name_base);
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



