#include "parallel_meshing_2d.hh"

#include <sys/types.h>
#include <sys/stat.h>

#include <iostream>
#include <fstream>

using namespace voro;

int main() {

    //custom poker shape contour lines def start
     double boundaryx[117]=
    {
       0.585 , 0.7651, 0.7255, 0.685 , 0.6545, 0.623 , 0.5926, 0.5716,
       0.5511, 0.5386, 0.5629, 0.587 , 0.6165, 0.652 , 0.6805, 0.731 ,
       0.7545, 0.7815, 0.813 , 0.843 , 0.8686, 0.8901, 0.9056, 0.9161,
       0.9186, 0.9176, 0.9136, 0.9001, 0.8826, 0.8621, 0.8385, 0.81  ,
       0.7825, 0.748 , 0.724 , 0.69  , 0.661 , 0.634 , 0.604 , 0.578 ,
       0.5376, 0.5641, 0.5811, 0.6075, 0.6596, 0.6891, 0.7046, 0.7126,
       0.7161, 0.7156, 0.7071, 0.6921, 0.6716, 0.644 , 0.6075, 0.57  ,
       0.5335, 0.4965, 0.4615, 0.424 , 0.3905, 0.3594, 0.3384, 0.3209,
       0.3119, 0.3099, 0.3134, 0.3199, 0.3284, 0.3429, 0.3579, 0.3765,
       0.4015, 0.4315, 0.4504, 0.4779, 0.4989, 0.4721, 0.4505, 0.421 ,
       0.392 , 0.356 , 0.317 , 0.283 , 0.2525, 0.218 , 0.189 , 0.1664,
       0.1424, 0.1259, 0.1174, 0.1149, 0.1159, 0.1224, 0.1354, 0.1514,
       0.1709, 0.197 , 0.221 , 0.2505, 0.282 , 0.308 , 0.3465, 0.3735,
       0.409 , 0.437 , 0.4641, 0.4984, 0.4849, 0.4634, 0.4424, 0.4165,
       0.3865, 0.347 , 0.3075, 0.266 , 0.585
    };
    double boundaryy[117]=
    {
       0.0284, 0.031 , 0.0381, 0.0521, 0.0676, 0.0906, 0.121 , 0.1495,
       0.1885, 0.2205, 0.188 , 0.1639, 0.1429, 0.1259, 0.1184, 0.1174,
       0.1189, 0.1269, 0.1419, 0.1629, 0.1885, 0.2195, 0.2525, 0.2895,
       0.328 , 0.3645, 0.386 , 0.4235, 0.4555, 0.482 , 0.5041, 0.5231,
       0.5356, 0.5451, 0.5481, 0.5466, 0.5401, 0.5296, 0.5121, 0.4906,
       0.4405, 0.501 , 0.5275, 0.5574, 0.6055, 0.6445, 0.6785, 0.7085,
       0.742 , 0.7775, 0.823 , 0.859 , 0.8905, 0.9206, 0.9456, 0.9616,
       0.9691, 0.9696, 0.9631, 0.9481, 0.9261, 0.895 , 0.864 , 0.824 ,
       0.7855, 0.7405, 0.7115, 0.6855, 0.663 , 0.6355, 0.6145, 0.5944,
       0.5734, 0.5534, 0.5335, 0.487 , 0.4355, 0.472 , 0.4951, 0.5176,
       0.5326, 0.5436, 0.5476, 0.5446, 0.5366, 0.5206, 0.5001, 0.477 ,
       0.4425, 0.4075, 0.378 , 0.35  , 0.2995, 0.2665, 0.2335, 0.2065,
       0.182 , 0.1574, 0.1409, 0.1274, 0.1189, 0.1159, 0.1179, 0.1239,
       0.1384, 0.1564, 0.1815, 0.2295, 0.1905, 0.149 , 0.1195, 0.0931,
       0.0706, 0.0496, 0.0366, 0.0294, 0.0284
    };
    std::vector<std::vector<double>> boundaries;
    std::vector<double> boundary1;
    for(int i=116;i>=0;i--){
        boundaryx[i]=boundaryx[i];
        boundaryy[i]=boundaryy[i];
        boundary1.push_back(boundaryx[i]);
        boundary1.push_back(boundaryy[i]);
    }
    boundaries.push_back(boundary1);

//------------------0.Set up parameters--------------------
    int num_t_max=8;
    int physical_core=4;
    int num_t_setup=physical_core;
    int num_t_setup_cshp=4; //for custom_shape_2d: since the boundary line segments count is small, we use serial code nt=1 to avoid parallel overhead
    int method_ind_0=2; //Yes: 0, 1, 3


   int method_ind=method_ind_0; //0 DM; 1 CVD; 2 Hybrid
   int Ntotal=2000; 
   double K=0.05; //0.005; 0:uniform sizing field
   int output_interval=1; //int: 0, no output; -1, last final output; 10, every 10 triangulations output

   char method_name[256];
   if(method_ind==0){
       sprintf(method_name,"dm");
   }else if(method_ind==1){
       sprintf(method_name,"cvd");
   }else{
       sprintf(method_name,"hybrid");
   }

   char case_name_base[256];
   sprintf(case_name_base,"test_fixed_pts_N_%d_K_%g_pinitFac_%g_addPFac_%g",Ntotal,K,pt_init_frac,add_pt_fac);

   mkdir(case_name_base,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);



   //------------------1.Create container----------------------
    printf("create container\n");
    int cnx=sqrt(Ntotal/3.3); int cny=cnx;
    container_2d con(0.0,1.0,0.0,1.0,cnx,cny,false,false,16,num_t_setup);



   //-------------------2.Create shape------------------------

      printf("create shape\n");

      //final test shape, complicated poker shape
      bool normalize_model=true;
      shape_2d_contour_lines shp1(con,num_t_setup,num_t_setup_cshp,boundaries,normalize_model);
      shape_2d_circle shp2(con,num_t_setup,0.1,0.5,0.7);
      shape_2d_difference shp(con,num_t_setup,&shp1,&shp2);
      

    //-------------------3.Create sizing field------------------
      printf("create size_field\n");
      sizing_2d_automatic size_field(&shp,K);


   //------------------4.Create pm2d----------------------
     srand(10);//set seed for rand() so that pt_init are the same across different num_t_mesh
     printf("create pm2d and pt init\n");
     parallel_meshing_2d pm2d(&con, &shp, &size_field, num_t_setup, output_interval, case_name_base);

     int Nfixed=10;
     double *fixed_pt_list=new double[2*Nfixed];
     double fixed_pt_x[Nfixed]={0.2,0.3,0.4,0.2,0.2,0.5,0.5,0.7,0.7,0.7};
     double fixed_pt_y[Nfixed]={0.2,0.3,0.4,0.3,0.4,0.6,0.1,0.4,0.3,0.2};
     for(int i=0;i<Nfixed;i++){
        fixed_pt_list[2*i]=fixed_pt_x[i];
        fixed_pt_list[2*i+1]=fixed_pt_y[i];
     }
     pm2d.add_fixed_points_normailze(Nfixed, &shp1, fixed_pt_list);
     pm2d.pt_init(Ntotal); //Ntotal=Nfix+Nmove: total number of fixed and unfixed points in the mesh

//======================Parallel meshing=====================
     
   if(method_ind==0){
         mesh_alg_2d_dm mesh_method(&pm2d);
         mesh_method.change_number_thread(num_t_setup);
         pm2d.meshing(&mesh_method); 
    }
    else if(method_ind==1){
         mesh_alg_2d_cvd mesh_method(&pm2d);
         mesh_method.change_number_thread(num_t_setup);
         pm2d.meshing(&mesh_method); 
    }
    else if(method_ind==2){
         mesh_alg_2d_hybrid mesh_method(&pm2d);
         mesh_method.change_number_thread(num_t_setup);
         pm2d.meshing(&mesh_method); 
    }

    return 0;   

}



