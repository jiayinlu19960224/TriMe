#include "parallel_meshing_2d.hh"

#include <sys/types.h>
#include <sys/stat.h>

#include <iostream>
#include <fstream>

using namespace voro;

int main() {

    int population_id=130;
    int species_id=4598506;

    char bug[256];
    sprintf(bug,"/Users/jiayinlu/Desktop/Jiayin/projects/wing/mapping/data");


// Read in data files of shape contour boundary, bdry landmark pts, and inner landmark pts
    std::vector<double> boundaryx;
    std::vector<double> boundaryy;
    std::vector<double> bdry_lm_x;
    std::vector<double> bdry_lm_y;
    std::vector<double> in_lm_x;
    std::vector<double> in_lm_y;
    std::vector<double> lm_x;
    std::vector<double> lm_y;

    std::fstream myfile;
    char f1n[256]; sprintf(f1n,"%s/population_%d+FMNH_%d_contour_x.txt",bug,population_id,species_id);
    myfile.open(f1n, std::ios_base::in);
    double ra;
    while (myfile >> ra)
    {
        boundaryx.push_back(ra); 
    }
    myfile.close();
    char f2n[256]; sprintf(f2n,"%s/population_%d+FMNH_%d_contour_y.txt",bug,population_id,species_id);
    myfile.open(f2n, std::ios_base::in);
    while (myfile >> ra)
    {
        boundaryy.push_back(ra); 
    }
    myfile.close();
    char f3n[256]; sprintf(f3n,"%s/population_%d+FMNH_%d_bdry_lm_x.txt",bug,population_id,species_id);
    myfile.open(f3n, std::ios_base::in);
    while (myfile >> ra)
    {
        bdry_lm_x.push_back(ra); 
        lm_x.push_back(ra); 
    }
    myfile.close();
    char f4n[256]; sprintf(f4n,"%s/population_%d+FMNH_%d_bdry_lm_y.txt",bug,population_id,species_id);
    myfile.open(f4n, std::ios_base::in);
    while (myfile >> ra)
    {
        bdry_lm_y.push_back(ra);
        lm_y.push_back(ra);
    }
    myfile.close();
    char f5n[256]; sprintf(f5n,"%s/population_%d+FMNH_%d_inner_lm_x.txt",bug,population_id,species_id);
    myfile.open(f5n, std::ios_base::in);
    while (myfile >> ra)
    {
        in_lm_x.push_back(ra); 
        lm_x.push_back(ra); 
    }
    myfile.close();
    char f6n[256]; sprintf(f6n,"%s/population_%d+FMNH_%d_inner_lm_y.txt",bug,population_id,species_id);
    myfile.open(f6n, std::ios_base::in);
    while (myfile >> ra)
    {
        in_lm_y.push_back(ra); 
        lm_y.push_back(ra); 
    }
    myfile.close();

    // Obtain Meshing shape contour input
    std::vector<std::vector<double>> boundaries;
    std::vector<double> boundary1;
    for(int i=0;i<(signed int)boundaryx.size();i++){
        boundary1.push_back(boundaryx[i]);
        boundary1.push_back(boundaryy[i]);
    }
    boundaries.push_back(boundary1);

    //Obtain Meshing fixed point input
    int Nfixed=(signed int)lm_x.size();
    double *fixed_pt_list=new double[2*Nfixed];
    for(int i=0;i<Nfixed;i++){
        fixed_pt_list[2*i]=lm_x[i];
        fixed_pt_list[2*i+1]=lm_y[i];
    }
    

//------------------0.Set up parameters--------------------
    int num_t_max=8;
    int physical_core=4;
    int num_t_setup=physical_core;
    int num_t_setup_cshp=4; //for custom_shape_2d: since the boundary line segments count is small, we use serial code nt=1 to avoid parallel overhead
    int method_ind_0=2; //Yes: 0, 1, 2


   int method_ind=method_ind_0; //0 DM; 1 CVD; 2 Hybrid
   int Ntotal=13000; 
   double K=0.0; //0.005; 0:uniform sizing field
   int output_interval=-1; //int: 0, no output; -1, last final output; 10, every 10 triangulations output

   char method_name[256];
   if(method_ind==0){
       sprintf(method_name,"dm");
   }else if(method_ind==1){
       sprintf(method_name,"cvd");
   }else{
       sprintf(method_name,"hybrid");
   }

   char case_name_base[256];
   sprintf(case_name_base,"population_%d+FMNH_%d_N_%d",population_id,species_id,Ntotal);

   mkdir(case_name_base,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);



   //------------------1.Create container----------------------
    printf("create container\n");
    int cnx=sqrt(Ntotal/3.3); int cny=cnx;
    container_2d con(-100.0,5000.0,-100.0,5000.0,cnx,cny,false,false,16,num_t_setup);



   //-------------------2.Create shape------------------------

      printf("create shape\n");

      //final test shape, complicated poker shape
      bool normalize_model=false;
      shape_2d_contour_lines shp(con,num_t_setup,num_t_setup_cshp,boundaries,normalize_model);
      

    //-------------------3.Create sizing field------------------
      printf("create size_field\n");
      sizing_2d_automatic size_field(&shp,K);


   //------------------4.Create pm2d----------------------
     srand(10);//set seed for rand() so that pt_init are the same across different num_t_mesh
     printf("create pm2d and pt init\n");
     parallel_meshing_2d pm2d(&con, &shp, &size_field, num_t_setup, output_interval, case_name_base);

     pm2d.add_fixed_points_normailze(Nfixed, &shp, fixed_pt_list);
     pm2d.pt_init(Ntotal); //Ntotal=Nfix+Nmove: total number of fixed and unfixed points in the mesh

//======================Parallel meshing=====================
     
   if(method_ind==0){
         mesh_alg_2d_dm mesh_method(&pm2d);
         mesh_method.change_number_thread(num_t_setup);
         pm2d.meshing(&mesh_method); 

        //Output bdry and inner landmark IDs to file
        int bdry_lm_ct=(signed int)bdry_lm_x.size();
        char bug0[256];
        sprintf(bug0,"%s/bdry_lm_id.txt",case_name_base);
        FILE *outFile0 = fopen(bug0, "a");
        for(int pi=0;pi<bdry_lm_ct;pi++){
            fprintf(outFile0,"%d %d\n",pi,pm2d.pt_ctgr[pi]);
        }
        fclose(outFile0);

        int in_lm_ct=(signed int)in_lm_x.size();
        char bug00[256];
        sprintf(bug00,"%s/inner_lm_id.txt",case_name_base);
        FILE *outFile00 = fopen(bug00, "a");
        for(int pi=0;pi<in_lm_ct;pi++){
            fprintf(outFile00,"%d %d\n",pi+bdry_lm_ct,pm2d.pt_ctgr[pi+bdry_lm_ct]);
        }
        fclose(outFile00);
    }
    else if(method_ind==1){
         mesh_alg_2d_cvd mesh_method(&pm2d);
         mesh_method.change_number_thread(num_t_setup);
         pm2d.meshing(&mesh_method); 

         //Output bdry and inner landmark IDs to file
        int bdry_lm_ct=(signed int)bdry_lm_x.size();
        char bug0[256];
        sprintf(bug0,"%s/bdry_lm_id.txt",case_name_base);
        FILE *outFile0 = fopen(bug0, "a");
        for(int pi=0;pi<bdry_lm_ct;pi++){
            fprintf(outFile0,"%d %d\n",pi,pm2d.pt_ctgr[pi]);
        }
        fclose(outFile0);

        int in_lm_ct=(signed int)in_lm_x.size();
        char bug00[256];
        sprintf(bug00,"%s/inner_lm_id.txt",case_name_base);
        FILE *outFile00 = fopen(bug00, "a");
        for(int pi=0;pi<in_lm_ct;pi++){
            fprintf(outFile00,"%d %d\n",pi+bdry_lm_ct,pm2d.pt_ctgr[pi+bdry_lm_ct]);
        }
        fclose(outFile00);
    }
    else if(method_ind==2){
         mesh_alg_2d_hybrid mesh_method(&pm2d);
         mesh_method.change_number_thread(num_t_setup);
         pm2d.meshing(&mesh_method); 

         //Output bdry and inner landmark IDs to file
        int bdry_lm_ct=(signed int)bdry_lm_x.size();
        char bug0[256];
        sprintf(bug0,"%s/bdry_lm_id.txt",case_name_base);
        FILE *outFile0 = fopen(bug0, "a");
        for(int pi=0;pi<bdry_lm_ct;pi++){
            fprintf(outFile0,"%d %d\n",pi,pm2d.pt_ctgr[pi]);
        }
        fclose(outFile0);

        int in_lm_ct=(signed int)in_lm_x.size();
        char bug00[256];
        sprintf(bug00,"%s/inner_lm_id.txt",case_name_base);
        FILE *outFile00 = fopen(bug00, "a");
        for(int pi=0;pi<in_lm_ct;pi++){
            fprintf(outFile00,"%d %d\n",pi+bdry_lm_ct,pm2d.pt_ctgr[pi+bdry_lm_ct]);
        }
        fclose(outFile00);
    }

    return 0;   

}



