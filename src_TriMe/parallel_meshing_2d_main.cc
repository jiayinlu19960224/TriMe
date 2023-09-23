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
        boundaryx[i]=boundaryx[i]*0.9+0.04;
        boundaryy[i]=boundaryy[i]*0.9+0.04;
        boundary1.push_back(boundaryx[i]);
        boundary1.push_back(boundaryy[i]);
    }
    boundaries.push_back(boundary1);
    //custom poker shape contour lines def end



//------------------0.Set up parameters--------------------
    int num_t_max=8;
    int physical_core=4;
    int num_t_setup=physical_core;
    int num_t_setup_cshp=4; //for custom_shape_2d: since the boundary line segments count is small, we use serial code nt=1 to avoid parallel overhead
    int method_ind_0=2;

    double shape_duration_n1=-1;
    double sizing_duration_n1=-1;
    double pm2d_duration_n1=-1;
    double pm2d_generate_pt_in_grid_n1=-1;

    //for(int num_t_setup=1;num_t_setup<=num_t_max;num_t_setup++){

       int method_ind=method_ind_0; //0 DM; 1 CVD; 2 Hybrid
       int Ntotal=10000; 
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
       double t0,t1,duration,duration_n1,per_iter_time,per_iter_time_n1,parallel_efficiency;
       double container_duration, shape_duration, sizing_duration,pm2d_duration;
       int tria_iter_ct, all_iter_ct;

       char case_name_base[256];
       sprintf(case_name_base,"paper_plot_adaptiveCVD_poker_N_%d_K_%g_pinitFac_%g_addPFac_%g",Ntotal,K,pt_init_frac,add_pt_fac);

       mkdir(case_name_base,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);



   //------------------1.Create container----------------------
        printf("create container\n");
        t0=omp_get_wtime();
        int cnx=sqrt(Ntotal/3.3); int cny=cnx;
        container_2d con(0.0,1.0,0.0,1.0,cnx,cny,false,false,16,num_t_setup);
        t1=omp_get_wtime(); 
        container_duration=t1-t0;
        printf("finished creating container: duration %g\n",t1-t0);


       


   //-------------------2.Create shape------------------------

      printf("create shape\n");
      t0=omp_get_wtime();
      
      //fist test shape, circle
      //shape_2d_circle shp(con,num_t_setup,0.4,0.5,0.5);

      //rectangle shape test
      //shape_2d_rectangle shp(con,num_t_setup,0.1,0.9,0.1,0.9);

      //final test shape, complicated poker shape
      bool normalize_model=true;
      shape_2d_contour_lines shp1(con,num_t_setup,num_t_setup_cshp,boundaries,normalize_model);
      shp1.print_geo_grid_sdf_to_file(case_name_base);


      shape_2d_circle shp2(con,num_t_setup,0.1,0.5,0.7);
      shape_2d_difference shp(con,num_t_setup,&shp1,&shp2);
      
      shp.print_geo_grid_sdf_to_file(case_name_base);

      t1=omp_get_wtime(); 
      shape_duration=t1-t0;
      printf("finished creating shape: duration %g\n",t1-t0);
     

   //-------------------3.Create sizing field------------------
      printf("create size_field\n");
      t0=omp_get_wtime();
      sizing_2d_automatic size_field(&shp,K);
      t1=omp_get_wtime(); 
      sizing_duration=t1-t0;
      printf("finished size_field: duration %g\n",t1-t0);

      size_field.print_fields_to_file(case_name_base);
      size_field.print_pts_to_file(case_name_base);




   //------------------4.Create pm2d----------------------
     srand(10);//set seed for rand() so that pt_init are the same across different num_t_mesh
     printf("create pm2d and pt init\n");
     t0=omp_get_wtime();
     parallel_meshing_2d pm2d(&con, &shp, &size_field, num_t_setup, output_interval, case_name_base);
     pm2d.pt_init(Ntotal);
     t1=omp_get_wtime(); 
     pm2d_duration=t1-t0;
     printf("finished pm2d and pt init: duration %g\n",t1-t0);

     double pm2d_generate_pt_in_grid=pm2d.t_generate_pt_in_grid;
     pm2d.print_init_pts_to_file(case_name_base);


     //Store initial pm2d status before meshing, for parallelization reset use later:
     int Ncurrent_0=pm2d.Ncurrent;
     int Nremain_0=pm2d.Nremain;
     int inner_pt_ct_0=pm2d.inner_pt_ct;
     int chrtrt_len_h_avg_0=pm2d.chrtrt_len_h_avg;
     
     double *chrtrt_len_h_0=new double[pm2d.gnxy];
     #pragma omp parallel for num_threads(num_t_setup)
     for(int i=0;i<pm2d.gnxy;i++){
         chrtrt_len_h_0[i]=pm2d.chrtrt_len_h[i];
     }
     double *xy_id_0=new double[Ntotal*2];
     int *pt_ctgr_0=new int[Ntotal];
     #pragma omp parallel for num_threads(num_t_setup)
     for(int i=0;i<Ntotal;i++){
         xy_id_0[2*i]=pm2d.xy_id[2*i];
         xy_id_0[2*i+1]=pm2d.xy_id[2*i+1];
         pt_ctgr_0[i]=pm2d.pt_ctgr[i];
     }
     double *bgrid_geps_0=new double[pm2d.geo_bgrid_ct];
     double *bgrid_deps_0=new double[pm2d.geo_bgrid_ct];
     #pragma omp parallel for num_threads(num_t_setup)
     for(int i=0;i<pm2d.geo_bgrid_ct;i++){
         bgrid_geps_0[i]=pm2d.bgrid_geps[i];
         bgrid_deps_0[i]=pm2d.bgrid_deps[i];
     }


   //-----------------------Output set up timing statistics-------------------
     
/*
     if(num_t_setup==1){
         shape_duration_n1=shape_duration;
         sizing_duration_n1=sizing_duration;
         pm2d_duration_n1=pm2d_duration;
         pm2d_generate_pt_in_grid_n1=pm2d_generate_pt_in_grid;

     }
     double shape_duration_efficiency=shape_duration_n1/(shape_duration*num_t_setup);
     double sizing_duration_efficiency=sizing_duration_n1/(sizing_duration*num_t_setup);
     double pm2d_duration_efficiency=pm2d_duration_n1/(pm2d_duration*num_t_setup);
     double pm2d_generate_pt_in_grid_efficiency=pm2d_generate_pt_in_grid_n1/(pm2d_generate_pt_in_grid*num_t_setup);
      char bug00[256];
      sprintf(bug00,"%s_time_stat_setup.txt",case_name_base);
      FILE *outFile00 = fopen(bug00, "a");
      fprintf(outFile00,"%d %g %g %g %g %g %g %g %g %g \n",num_t_setup_temp, container_duration, 
         shape_duration, shape_duration_efficiency,
         sizing_duration, sizing_duration_efficiency,
         pm2d_duration, pm2d_duration_efficiency,
         pm2d_generate_pt_in_grid, pm2d_generate_pt_in_grid_efficiency);
      fclose(outFile00);
      
*/
//   }
        

//======================Parallel meshing=====================
     



   //warm up round: No need to store statistics for this
   printf("-----start warm up round method %d -----\n", method_ind);
   t0=omp_get_wtime();
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

    t1=omp_get_wtime();
    printf("-----finished warm up round: duration %g----\n",t1-t0);


   printf("-----start parallelization-----\n");
   for(method_ind=2;method_ind>=0;method_ind--){

     printf("create container\n");
     container_2d con_p(0.0,1.0,0.0,1.0,cnx,cny,false,false,16,num_t_setup);
       
     pm2d.con=&con_p;



     
       if(method_ind==0){
           sprintf(method_name,"dm");
       }else if(method_ind==1){
           sprintf(method_name,"cvd");
       }else{
           sprintf(method_name,"hybrid");
       }
       sprintf(case_name_base,"final_test_adaptiveCVD_poker_N_%d_K_%g_%s_pinitFac_%g_addPFac_%g",Ntotal,K,method_name,pt_init_frac,add_pt_fac);

       //Timing variables
        int hybrid_cvd_ct=0;
        int tria_iter_ct_n1;
        int all_iter_ct_n1;
        double t_meshing_init_n1;
        double t_add_pt_n1;
        double t_voro_computation_n1;
        double t_update_pt_position_n1;
        double t_algorithm_remain_n1;
        double t_algorithm_all_non_voro_n1;
        double t_add_pt_centroid_n1;

        double t_dm_getBarinfo_n1;
        double t_dm_applyBarForce_n1;
        double t_dm_updatePtPos_n1;

        double t_meshing_init;
        double t_add_pt;
        double t_voro_computation;
        double t_update_pt_position;
        double t_algorithm_remain;
        double t_algorithm_all_non_voro;
        double t_add_pt_centroid;

        double t_dm_getBarinfo;
        double t_dm_applyBarForce;
        double t_dm_updatePtPos;


      for(int num_t_mesh=num_t_max;num_t_mesh<=num_t_max+1;num_t_mesh++){
         if(num_t_mesh==num_t_max+1){pm2d.output_interval=1;output_interval=1;}
         else{pm2d.output_interval=0;output_interval=0;}


           //Create output directory
           char case_name[256];
           sprintf(case_name,"%s_nt_%d",case_name_base,num_t_mesh);
           //if there are outputs
           if(output_interval!=0){
               // Make the output directory if it doesn't already exist
               mkdir(case_name,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
           }

           //Update container and pm2d new num_t
            con_p.change_number_thread(num_t_mesh);
            pm2d.num_t=num_t_mesh;
            pm2d.file_name_prefix=case_name;

   //=====================if num_t_mesh>1: Reset pm2d to before mesh status==========
   //=====================And add initial points to container==========
   //reset: 
            //arrays: xy_id, pt_ctgr, chrtrt_len_h, bgrid_geps, bgrid_deps, 
            //variables: Ncurrent, Nremain, inner_pt_ct, chrtrt_len_h_avg
            //Reset pm2d
            pm2d.Ncurrent=Ncurrent_0;
            pm2d.Nremain=Nremain_0;
            pm2d.inner_pt_ct=inner_pt_ct_0;
            pm2d.chrtrt_len_h_avg=chrtrt_len_h_avg_0;

           #pragma omp parallel for num_threads(num_t_setup)
           for(int i=0;i<pm2d.gnxy;i++){
               pm2d.chrtrt_len_h[i]=chrtrt_len_h_0[i];
           }
           
           #pragma omp parallel for num_threads(num_t_setup)
           for(int i=0;i<Ntotal;i++){
               pm2d.xy_id[2*i]=xy_id_0[2*i];
               pm2d.xy_id[2*i+1]=xy_id_0[2*i+1];
               pm2d.pt_ctgr[i]=pt_ctgr_0[i];
           }
           
           #pragma omp parallel for num_threads(num_t_setup)
           for(int i=0;i<pm2d.geo_bgrid_ct;i++){
               pm2d.bgrid_geps[i]=bgrid_geps_0[i];
               pm2d.bgrid_deps[i]=bgrid_deps_0[i];
           }

           //Put initial points into container
            con_p.clear();
            con_p.add_parallel(pm2d.xy_id, pm2d.Ncurrent, num_t_mesh);
            con_p.put_reconcile_overflow();

          
   //------------------5.Start meshing--------------------
          printf("start meshing\n");
          if(method_ind==0){
               mesh_alg_2d_dm mesh_method(&pm2d);
               //mesh_method.change_number_thread(num_t_mesh);
               t0=omp_get_wtime();
               pm2d.meshing(&mesh_method); 
               t1=omp_get_wtime();
               tria_iter_ct=mesh_method.tria_iter_ct;
               all_iter_ct=mesh_method.all_iter_ct;

               t_meshing_init=mesh_method.t_meshing_init;
               t_add_pt=mesh_method.t_add_pt;
               t_voro_computation=mesh_method.t_voro_computation;
               t_update_pt_position=mesh_method.t_update_pt_position;
               t_algorithm_remain=mesh_method.t_algorithm_remain;
               t_algorithm_all_non_voro=mesh_method.t_algorithm_all_non_voro;
               t_add_pt_centroid=mesh_method.t_add_pt_centroid;

               t_dm_getBarinfo=mesh_method.t_dm_getBarinfo;
               t_dm_applyBarForce=mesh_method.t_dm_applyBarForce;
               t_dm_updatePtPos=mesh_method.t_dm_updatePtPos;
          }
          else if(method_ind==1){
               mesh_alg_2d_cvd mesh_method(&pm2d);
               //mesh_method.change_number_thread(num_t_mesh);
               t0=omp_get_wtime();
               pm2d.meshing(&mesh_method); 
               t1=omp_get_wtime();
               tria_iter_ct=mesh_method.tria_iter_ct;
               all_iter_ct=mesh_method.all_iter_ct;

               t_meshing_init=mesh_method.t_meshing_init;
               t_add_pt=mesh_method.t_add_pt;
               t_voro_computation=mesh_method.t_voro_computation;
               t_update_pt_position=mesh_method.t_update_pt_position;
               t_algorithm_remain=mesh_method.t_algorithm_remain;
               t_algorithm_all_non_voro=mesh_method.t_algorithm_all_non_voro;
               t_add_pt_centroid=mesh_method.t_add_pt_centroid;

               t_dm_getBarinfo=-1;
               t_dm_applyBarForce=-1;
               t_dm_updatePtPos=-1;
          }
          else if(method_ind==2){
               mesh_alg_2d_hybrid mesh_method(&pm2d);
               //mesh_method.change_number_thread(num_t_mesh);
               t0=omp_get_wtime();
               pm2d.meshing(&mesh_method); 
               t1=omp_get_wtime();
               tria_iter_ct=mesh_method.tria_iter_ct;
               all_iter_ct=mesh_method.all_iter_ct;

               t_meshing_init=mesh_method.t_meshing_init;
               t_add_pt=mesh_method.t_add_pt;
               t_voro_computation=mesh_method.t_voro_computation;
               t_update_pt_position=mesh_method.t_update_pt_position;
               t_algorithm_remain=mesh_method.t_algorithm_remain;
               t_algorithm_all_non_voro=mesh_method.t_algorithm_all_non_voro;
               t_add_pt_centroid=mesh_method.t_add_pt_centroid;

               t_dm_getBarinfo=mesh_method.t_dm_getBarinfo;
               t_dm_applyBarForce=mesh_method.t_dm_applyBarForce;
               t_dm_updatePtPos=mesh_method.t_dm_updatePtPos;
               hybrid_cvd_ct=mesh_method.cvd_iter_ct;
               
          }
          
          
          duration=t1-t0;
          per_iter_time=(duration-t_meshing_init)/tria_iter_ct;
          
          if(num_t_mesh==1){
            tria_iter_ct_n1=tria_iter_ct;
            all_iter_ct_n1=all_iter_ct;
            duration_n1=duration;
            per_iter_time_n1=per_iter_time;
            t_meshing_init_n1=t_meshing_init;
            t_add_pt_n1=t_add_pt;
            t_voro_computation_n1=t_voro_computation;
            t_update_pt_position_n1=t_update_pt_position;
            t_algorithm_remain_n1=t_algorithm_remain;
            t_algorithm_all_non_voro_n1=t_algorithm_all_non_voro;
            t_add_pt_centroid_n1=t_add_pt_centroid;

            t_dm_getBarinfo_n1=t_dm_getBarinfo;
            t_dm_applyBarForce_n1=t_dm_applyBarForce;
            t_dm_updatePtPos_n1=t_dm_updatePtPos;
         }

          parallel_efficiency=duration_n1
                  /
                  (
                     (per_iter_time*tria_iter_ct_n1+t_meshing_init)
                     *num_t_mesh
                  );
          double e_meshing_init=t_meshing_init_n1/(t_meshing_init*num_t_mesh);
          double e_add_pt=t_add_pt_n1/(t_add_pt*num_t_mesh);
          double e_voro_computation=(t_voro_computation_n1/tria_iter_ct_n1)/((t_voro_computation/tria_iter_ct)*num_t_mesh);
          double e_update_pt_position=(t_update_pt_position_n1/all_iter_ct_n1)/((t_update_pt_position/all_iter_ct)*num_t_mesh);
          double e_algorithm_remain=(t_algorithm_remain_n1/all_iter_ct_n1)/((t_algorithm_remain/all_iter_ct)*num_t_mesh);
          double e_algorithm_all_non_voro=
                  (t_algorithm_all_non_voro_n1)
                  /
                  (
                     ((t_algorithm_all_non_voro-t_meshing_init)/all_iter_ct*all_iter_ct_n1+t_meshing_init)
                     *num_t_mesh
                  );
         double e_add_pt_centroid=t_add_pt_centroid_n1/(t_add_pt_centroid*num_t_mesh);

         double e_dm_getBarinfo=t_dm_getBarinfo_n1/(t_dm_getBarinfo*num_t_mesh);
         double e_dm_applyBarForce=t_dm_applyBarForce_n1/(t_dm_applyBarForce*num_t_mesh);
         double e_dm_updatePtPos=t_dm_updatePtPos_n1/(t_dm_updatePtPos*num_t_mesh);
          printf("finished meshing: duration %g\n",duration);




   //------------------6.Output timing statistics-------------
          char bug0[256];
          sprintf(bug0,"%s_time_stat.txt",case_name_base);
          FILE *outFile0 = fopen(bug0, "a");
          fprintf(outFile0,"%d %d %d "
            "%g %g %g "

            "%g %g "
            "%g %g "
            "%g %g "
            "%g %g "
            "%g %g "
            "%g %g "

            "%g %g "
            "%g %g "
            "%g %g "

            "%g %g "
            "%d \n",
            num_t_mesh, tria_iter_ct, all_iter_ct,
            duration, per_iter_time, parallel_efficiency,

            t_meshing_init, e_meshing_init,
            t_add_pt, e_add_pt,
            t_voro_computation, e_voro_computation,
            t_update_pt_position, e_update_pt_position,
            t_algorithm_remain, e_algorithm_remain,
            t_algorithm_all_non_voro,e_algorithm_all_non_voro,

            t_dm_getBarinfo, e_dm_getBarinfo,
            t_dm_applyBarForce,e_dm_applyBarForce,
            t_dm_updatePtPos,e_dm_updatePtPos,

            t_add_pt_centroid,e_add_pt_centroid, 

            hybrid_cvd_ct
            );
          fclose(outFile0);

      }
   }


   delete [] chrtrt_len_h_0;
   delete [] xy_id_0;
   delete [] pt_ctgr_0;
   delete [] bgrid_geps_0;
   delete [] bgrid_deps_0;


    return 0;   

}



