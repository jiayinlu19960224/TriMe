#include "parallel_meshing_2d.hh"

#include <sys/types.h>
#include <sys/stat.h>

#include <iostream>
#include <fstream>

#include <sstream>

using namespace voro;

int main() {

//------------------Define a shape from custom contour line segments--------------------
    //Variables to store the line segments from multiple boundaries
    vector<vector<double>> boundaries;
    vector<vector<double>> boundaries_x;
    vector<vector<double>> boundaries_y;

    //Read in North America map boundaries from .txt file
    string fileName = "north_america_single_bdry_xy.txt"; 
    ifstream inFile(fileName);
    string line;
    int l_ct=0;
    while(getline(inFile,line)){
        l_ct++;
        if(!line.empty()){
            istringstream lstream(line);
        
            if(l_ct%2!=0){ //x
                vector<double> boundary_i_x;
                double x_temp;
                while(lstream >> x_temp){
                    boundary_i_x.push_back(x_temp);
                }
                boundaries_x.push_back(boundary_i_x);
            }
            else{
                vector<double> boundary_i_y;
                double y_temp;
                while(lstream >> y_temp){
                    boundary_i_y.push_back(y_temp);
                }
                boundaries_y.push_back(boundary_i_y);
            }

            
        }
    }

    //Store multiple boundaries to variable ``boundaries`` in the proper way
    int bdry_ct=0;
    for(int bi=0;bi<(signed int)boundaries_x.size();bi++){
        vector<double> boundary_i_x=boundaries_x[bi];
        vector<double> boundary_i_y=boundaries_y[bi];
        vector<double> boundary_i_xy;
        for(int bic=0;bic<(signed int)boundary_i_x.size();bic++){
            boundary_i_xy.push_back(boundary_i_x[bic]);
            boundary_i_xy.push_back(boundary_i_y[bic]);
            bdry_ct++;
        }
        
        boundaries.push_back(boundary_i_xy);
    }
    printf("Number of boundaries: %d\n",(signed int)boundaries.size());
//------------------0.Set up parameters--------------------
    int physical_core=4; //Specify number of physical cores in the computer. 
    int method_ind=2; //Method index, 0 DM; 1 CVD; 2 Hybrid
    int Ntotal=20000; //Number of vertices in the mesh
    double K=0.2; //0:uniform sizing field
    int output_interval=-1; //File output, int: 0, no output; -1, last final output; 10, every 10 triangulations output


    int num_t_meshing=2*physical_core;
    int num_t_setup=physical_core;
    //for custom_shape_2d: if the boundary contour line segments count is small, we use serial code nt=1 to avoid parallel overhead
    int num_t_setup_cshp=1; 
    if(bdry_ct>200){
        num_t_setup_cshp=num_t_setup;
    }

   char method_name[256];
   if(method_ind==0){
       sprintf(method_name,"dm");
   }else if(method_ind==1){
       sprintf(method_name,"cvd");
   }else{
       sprintf(method_name,"hybrid");
   }

   //File output name prefix
   char case_name_base[256];
   sprintf(case_name_base,"north_america_map_mesh_N_%d_K_%g",Ntotal,K);
   //Create directory to store file outputs
   mkdir(case_name_base,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

   //------------------1.Create container----------------------
    printf("create container\n");
    int cnx=sqrt(Ntotal/3.3); int cny=cnx;
    container_2d con(0.0,1.0,0.0,1.0,cnx,cny,false,false,16,num_t_setup);


   //-------------------2.Create shape------------------------
      printf("create shape\n");
      //final test shape, complicated poker shape
      bool normalize_model=true;
      shape_2d_contour_lines shp(con,num_t_setup,num_t_setup_cshp,boundaries,normalize_model);

    //-------------------3.Create sizing field------------------
      printf("create size_field\n");
      sizing_2d_automatic size_field(&shp,K);


   //------------------4.Create pm2d----------------------
     srand(10);//set seed for rand() so that pt_init are the same across different num_t_mesh
     printf("create pm2d and pt init\n");
     parallel_meshing_2d pm2d(&con, &shp, &size_field, num_t_setup, output_interval, case_name_base);
     pm2d.pt_init(Ntotal); //Ntotal=Nfix+Nmove: total number of fixed and unfixed points in the mesh

//======================Parallel meshing=====================
     
    printf("meshing\n");
   if(method_ind==0){
         mesh_alg_2d_dm mesh_method(&pm2d);
         mesh_method.change_number_thread(num_t_meshing);
         pm2d.meshing(&mesh_method); 
    }
    else if(method_ind==1){
         mesh_alg_2d_cvd mesh_method(&pm2d);
         mesh_method.change_number_thread(num_t_meshing);
         pm2d.meshing(&mesh_method); 
    }
    else if(method_ind==2){
         mesh_alg_2d_hybrid mesh_method(&pm2d);
         mesh_method.change_number_thread(num_t_meshing);
         pm2d.meshing(&mesh_method); 
    }

    return 0;   
}



