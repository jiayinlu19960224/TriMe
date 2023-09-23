#include "custom_shape_3d.hh"

using namespace voro;

custom_shape_3d::custom_shape_3d(std::string file_name, double ax_, double bx_, double ay_, double by_, double az_, double bz_, 
        int nx_, int ny_, int nz_):ax(ax_), bx(bx_), ay(ay_), by(by_), az(az_), bz(bz_), nx(nx_),ny(ny_),nz(nz_),
        dx(1.0/nx*(bx-ax)), dy(1.0/ny*(by-ay)), dz(1.0/nz*(bz-az))
    {
        //1. get triangles and normals
        ReadSTL(file_name);
        //2. get vertices
        get_vertex();
        //3.get elements
        get_elements();
        //4.get edges, num_edges
        CreateEdges();
        //5.get ver_tria
        get_ver_tria();
        //6.get edge_tria
        get_edge_tria();
        //7.get cell_tria
        get_cell_tria();
    }


//Note: triangle v0,v1,v2 are all translated to be centered at (0,0,0)
//used in get_cell_tria, in testing if cube and triangle intersect
//SAT test 3
bool custom_shape_3d::tria_cube_test3(double ax, double ay, double az, 
        double hx, double hy, double hz, 
        double v0x, double v0y, double v0z, double v1x, double v1y, double v1z, double v2x, double v2y, double v2z){
    //1. find translated tria vertex on aij
    double p0=ax*v0x+ay*v0y+az*v0z;
    double p1=ax*v1x+ay*v1y+az*v1z;
    double p2=ax*v2x+ay*v2y+az*v2z;
    double minp=std::min(p0,std::min(p1,p2));
    double maxp=std::max(p0,std::max(p1,p2));
    
    //2. largest "radius" of box projection on aij: r
    double r=hx*abs(ax)+hy*abs(ay)+hz*abs(az);
    
    //3. compare and decide:
    if(minp>r+1e-7 || maxp<-r-1e-7){
        //no intersection; no further test needed; Do not put tria i in the cube cell
        return false;
    }
    else{
        return true;
    } 
}


//get cell_tria: list of triangles in each grid cell, stored in cell_tria, index by cell i, j, k: cell_tria[nx*ny*k+nx*j+i]
//Accurate list of triangles each cell, 
//Use SAT(Separating Axis Theorem) for triangle-cube intersection detection
void custom_shape_3d::get_cell_tria(){
    cell_tria.resize(nx*ny*nz);
    //loop through each triangle, 
    for(int i=0;i<triangles.size();i++){
        //triangle face normal
        double a=normals[i][0];
        double b=normals[i][1];
        double c=normals[i][2]; 
        double absa=abs(a);
        double absb=abs(b);
        double absc=abs(c);
        //get xmin, xmax, ymin, ymax, zmin, zmax
        double v1x=triangles[i][0][0];
        double v1y=triangles[i][0][1];
        double v1z=triangles[i][0][2];
        double v2x=triangles[i][1][0];
        double v2y=triangles[i][1][1];
        double v2z=triangles[i][1][2];
        double v3x=triangles[i][2][0];
        double v3y=triangles[i][2][1];
        double v3z=triangles[i][2][2];
        double vxmin=std::min(v1x,std::min(v2x,v3x));
        double vymin=std::min(v1y,std::min(v2y,v3y));
        double vzmin=std::min(v1z,std::min(v2z,v3z));
        double vxmax=std::max(v1x,std::max(v2x,v3x));
        double vymax=std::max(v1y,std::max(v2y,v3y));
        double vzmax=std::max(v1z,std::max(v2z,v3z));
        //tria plane d for a*x+b*y+c*z+d=0
        double d=-a*v1x-b*v1y-c*v1z;
        
        //Minimum bounding box of the triangle
        //will need to test cubes in this bounding box for triangle-cube intersection
        int iil=(vxmin-ax)/dx-1; int iih=(vxmax-ax)/dx+1;  if(iil<0){iil=0;} if(iih>=nx){iih=nx-1;} 
        int jjl=(vymin-ay)/dy-1; int jjh=(vymax-ay)/dy+1;  if(jjl<0){jjl=0;} if(jjh>=ny){jjh=ny-1;} 
        int kkl=(vzmin-az)/dz-1; int kkh=(vzmax-az)/dz+1;  if(kkl<0){kkl=0;} if(kkh>=nz){kkh=nz-1;} 
        
        //loop through the above bounidng box cells
        //and detect of triangle intersect with the cube.
        for(int ii=iil;ii<=iih;ii++){
            for(int jj=jjl;jj<=jjh;jj++){
                for(int kk=kkl;kk<=kkh;kk++){
                    //bmin
                    double bminx=ax+dx*ii;
                    double bminy=ay+dy*jj;
                    double bminz=az+dz*kk;
                    //bmax
                    double bmaxx=bminx+dx;
                    double bmaxy=bminy+dy;
                    double bmaxz=bminz+dz;
                    //center c
                    double cx=0.5*(bminx+bmaxx);
                    double cy=0.5*(bminy+bmaxy);
                    double cz=0.5*(bminz+bmaxz);
                    //positive half diagonal h
                    double hx=0.5*dx;
                    double hy=0.5*dy;
                    double hz=0.5*dz;
                    //test 1 AABB against minimumBounding Box of Triangle
                    //1.1: e0=(1,0,0)
                    if(!(vxmin>bmaxx+1e-7 || vxmax<bminx-1e-7)){
                        //1.2: e1=(0,1,0)
                        if(!(vymin>bmaxy+1e-7 || vymax<bminy-1e-7)){
                            if(!(vzmin>bmaxz+1e-7 || vzmax<bminz-1e-7)){

                                //test 2: test if AABB cube overlap with plane defined by tria
                                //use fast plane/AABB overlap test as described in Real-Time Rendering 4th Edition 22.10.1
                                //extent e: largest box projection on tria plane normal n
                                //tria plane normal n=(a,b,c)
                                double e=hx*absa+hy*absb+hz*absc;
                                //absolute distance of center c to plane: s
                                double s=a*cx+b*cy+c*cz+d;
                                s=abs(s);
                                //if s>e: no intersection of cube and tria plane,
                                //        that is, no intersection of cube and triangle
                                //        therefore can end testing process and NOT add i to this cell
                                //else, cube intersect with plane, need to proceed with testing process
                                if(s<=e+1e-7){  //continue testing
                                    //translated tria vertex so that cube is centered at (0,0,0)
                                    double tv0x=v1x-cx; double tv0y=v1y-cy; double tv0z=v1z-cz;
                                    double tv1x=v2x-cx; double tv1y=v2y-cy; double tv1z=v2z-cz;
                                    double tv2x=v3x-cx; double tv2y=v3y-cy; double tv2z=v3z-cz;

                                    //calculate f0,f1,f2
                                    double f0x=tv1x-tv0x; double f0y=tv1y-tv0y; double f0z=tv1z-tv0z;
                                    double f1x=tv2x-tv1x; double f1y=tv2y-tv1y; double f1z=tv2z-tv1z;
                                    double f2x=tv0x-tv2x; double f2y=tv0y-tv2y; double f2z=tv0z-tv2z;

                                    double e0x=1.0; double e0y=0.0; double e0z=0.0;
                                    double e1x=0.0; double e1y=1.0; double e1z=0.0;
                                    double e2x=0.0; double e2y=0.0; double e2z=1.0;

                                    //Test 3.1 out of 9: a00=e0xf0
                                    double a00x=e0y*f0z-e0z*f0y;
                                    double a00y=e0z*f0x-e0x*f0z;
                                    double a00z=e0x*f0y-e0y*f0x;
                                    if(tria_cube_test3(a00x,a00y,a00z,hx,hy,hz,tv0x,tv0y,tv0z,tv1x,tv1y,tv1z,tv2x,tv2y,tv2z)){
                                        //Test 3.2 out of 9: a01=e0xf1
                                        double a01x=e0y*f1z-e0z*f1y;
                                        double a01y=e0z*f1x-e0x*f1z;
                                        double a01z=e0x*f1y-e0y*f1x;
                                        if(tria_cube_test3(a01x,a01y,a01z,hx,hy,hz,tv0x,tv0y,tv0z,tv1x,tv1y,tv1z,tv2x,tv2y,tv2z)){
                                            //Test 3.3 out of 9: a02=e0xf2
                                            double a02x=e0y*f2z-e0z*f2y;
                                            double a02y=e0z*f2x-e0x*f2z;
                                            double a02z=e0x*f2y-e0y*f2x;
                                            if(tria_cube_test3(a02x,a02y,a02z,hx,hy,hz,tv0x,tv0y,tv0z,tv1x,tv1y,tv1z,tv2x,tv2y,tv2z)){
                                                //Test 3.4 out of 9: a10=e1xf0
                                                double a10x=e1y*f0z-e1z*f0y;
                                                double a10y=e1z*f0x-e1x*f0z;
                                                double a10z=e1x*f0y-e1y*f0x;
                                                if(tria_cube_test3(a10x,a10y,a10z,hx,hy,hz,tv0x,tv0y,tv0z,tv1x,tv1y,tv1z,tv2x,tv2y,tv2z)){
                                                    //Test 3.5 out of 9: a11=e1xf1
                                                    double a11x=e1y*f1z-e1z*f1y;
                                                    double a11y=e1z*f1x-e1x*f1z;
                                                    double a11z=e1x*f1y-e1y*f1x;
                                                    if(tria_cube_test3(a11x,a11y,a11z,hx,hy,hz,tv0x,tv0y,tv0z,tv1x,tv1y,tv1z,tv2x,tv2y,tv2z)){
                                                        //Test 3.6 out of 9: a12=e1xf2
                                                        double a12x=e1y*f2z-e1z*f2y;
                                                        double a12y=e1z*f2x-e1x*f2z;
                                                        double a12z=e1x*f2y-e1y*f2x;
                                                        if(tria_cube_test3(a12x,a12y,a12z,hx,hy,hz,tv0x,tv0y,tv0z,tv1x,tv1y,tv1z,tv2x,tv2y,tv2z)){
                                                            //Test 3.7 out of 9: a20=e2xf0
                                                            double a20x=e2y*f0z-e2z*f0y;
                                                            double a20y=e2z*f0x-e2x*f0z;
                                                            double a20z=e2x*f0y-e2y*f0x;
                                                            if(tria_cube_test3(a20x,a20y,a20z,hx,hy,hz,tv0x,tv0y,tv0z,tv1x,tv1y,tv1z,tv2x,tv2y,tv2z)){
                                                                //Test 3.8 out of 9: a21=e2xf1
                                                                double a21x=e2y*f1z-e2z*f1y;
                                                                double a21y=e2z*f1x-e2x*f1z;
                                                                double a21z=e2x*f1y-e2y*f1x;
                                                                if(tria_cube_test3(a21x,a21y,a21z,hx,hy,hz,tv0x,tv0y,tv0z,tv1x,tv1y,tv1z,tv2x,tv2y,tv2z)){
                                                                    //Test 3.9 out of 9: a22=e2xf2
                                                                    double a22x=e2y*f2z-e2z*f2y;
                                                                    double a22y=e2z*f2x-e2x*f2z;
                                                                    double a22z=e2x*f2y-e2y*f2x;
                                                                    if(tria_cube_test3(a22x,a22y,a22z,hx,hy,hz,tv0x,tv0y,tv0z,tv1x,tv1y,tv1z,tv2x,tv2y,tv2z)){
                                                                        cell_tria[nx*ny*kk+nx*jj+ii].push_back(i);
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

//Read in a .stl file and store the triangles and the triangle normals to 
//"triangles", and "normals"
void custom_shape_3d::ReadSTL(std::string file_name) {
    FILE* fp = std::fopen(file_name.c_str(), "r");
    if (fp == NULL) {
        printf("No STL file found\n");
        std::cout << "ERROR: cannot read " << file_name << std::endl;
        return;
    }
    
    triangles.clear();
    normals.clear();
    char input[80];
    for (;;) {
        fscanf(fp, "%s", input);
        if (input == std::string("endsolid")) {
            // reach end of file
            break;
        }
        for (;input != std::string("facet");) {
            fscanf(fp, "%s", input);
        }
        std::vector<std::vector<double> > triangle;
        std::vector<double> normal(3); 
        if (std::is_same<double, float>::value) {
            float nx, ny, nz;
            fscanf(fp, "%s %f %f %f\n", input, &nx, &ny, &nz);
            normal[0] = nx; normal[1] = ny; normal[2] = nz;
        }
        else 
            fscanf(fp, "%s %lf %lf %lf\n", input, &normal[0], &normal[1], &normal[2]);
        fscanf(fp, "%s %s", input, input);
        triangle.clear();
        for (int i = 0;i < 3;++i) {
            std::vector<double> p(3);
            if (std::is_same<double, float>::value) {
                float px, py, pz;
                fscanf(fp, "%s %f %f %f\n", input, &px, &py, &pz);
                p[0] = px; p[1] = py; p[2] = pz;
            }
            else
                fscanf(fp, "%s %lf %lf %lf\n", input, &p[0], &p[1], &p[2]);
            triangle.push_back(p);
        }
        fscanf(fp, "%s %s", input, input);
        triangles.push_back(triangle);
        //calculate correct normal
        double ABx=triangle[1][0]-triangle[0][0];
        double ABy=triangle[1][1]-triangle[0][1];
        double ABz=triangle[1][2]-triangle[0][2];
        double ACx=triangle[2][0]-triangle[0][0];
        double ACy=triangle[2][1]-triangle[0][1];
        double ACz=triangle[2][2]-triangle[0][2];
        normal[0]=ABy*ACz-ABz*ACy;
        normal[1]=-(ABx*ACz-ABz*ACx);
        normal[2]=ABx*ACy-ABy*ACx;
        double nlen=sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
        normal[0]=normal[0]/nlen;
        normal[1]=normal[1]/nlen;
        normal[2]=normal[2]/nlen;
        normals.push_back(normal);
    }
    fclose(fp);
    normalize_model(triangles);
}


//normalize and center the model so that the model is in [ax+0.1*(bx-ax), ax+0.9*(bx-ax)]x[ay+0.1*(by-ay),ay+0.9*(by-ay)]x[az+0.1*(bz-az), az+0.9*(bz-az)]
void custom_shape_3d::normalize_model(std::vector<std::vector<std::vector<double> > >& _triangles){
    double xmin=100000;
    double xmax=-100000;
    double ymin=100000;
    double ymax=-100000;
    double zmin=100000;
    double zmax=-100000;
    
    //loop through all triangles
    for (int ii = 0;ii < (signed int)_triangles.size();ii++) {
        std::vector<double> v0 = _triangles[ii][0];//three vertices of a triangle
        std::vector<double> v1 = _triangles[ii][1];
        std::vector<double> v2 = _triangles[ii][2];
        
        if (v0[0] < xmin) { xmin = v0[0]; } 
        if (v1[0] < xmin) { xmin = v1[0]; } 
        if (v2[0] < xmin) { xmin = v2[0]; }
         
        if (v0[0] > xmax) { xmax = v0[0]; } 
        if (v1[0] > xmax) { xmax = v1[0]; } 
        if (v2[0] > xmax) { xmax = v2[0]; }
        
        if (v0[1] < ymin) { ymin = v0[1]; } 
        if (v1[1] < ymin) { ymin = v1[1]; } 
        if (v2[1] < ymin) { ymin = v2[1]; }
         
        if (v0[1] > ymax) { ymax = v0[1]; } 
        if (v1[1] > ymax) { ymax = v1[1]; } 
        if (v2[1] > ymax) { ymax = v2[1]; }
      
        if (v0[2] < zmin) { zmin = v0[2]; } 
        if (v1[2] < zmin) { zmin = v1[2]; } 
        if (v2[2] < zmin) { zmin = v2[2]; }
        
        if (v0[2] > zmax) { zmax = v0[2]; } 
        if (v1[2] > zmax) { zmax = v1[2]; } 
        if (v2[2] > zmax) { zmax = v2[2]; }
    }
    
    double max_range=xmax-xmin;
    if(ymax-ymin>max_range){max_range=ymax-ymin;}
    if(zmax-zmin>max_range){max_range=zmax-zmin;}
    double xmid=0.5*(xmin+xmax);
    double ymid=0.5*(ymin+ymax);
    double zmid=0.5*(zmin+zmax);
    
    double min_domain_range=std::min(bx-ax,std::min(by-ay,bz-az));
    
    //loop through all triangles,
    //normalize them so that the model fix into domain box, biggest dimension of moedel is 0.8 of the size of the minimum length side of the box
    //centered at the box midpoint.
    for (int ii = 0;ii < (signed int)_triangles.size();ii++) {
        _triangles[ii][0][0]=min_domain_range*0.8*(_triangles[ii][0][0]-xmid)/max_range+0.5*(bx-ax);
        _triangles[ii][1][0]=min_domain_range*0.8*(_triangles[ii][1][0]-xmid)/max_range+0.5*(bx-ax);
        _triangles[ii][2][0]=min_domain_range*0.8*(_triangles[ii][2][0]-xmid)/max_range+0.5*(bx-ax);
        
        _triangles[ii][0][1]=min_domain_range*0.8*(_triangles[ii][0][1]-ymid)/max_range+0.5*(by-ay);
        _triangles[ii][1][1]=min_domain_range*0.8*(_triangles[ii][1][1]-ymid)/max_range+0.5*(by-ay);
        _triangles[ii][2][1]=min_domain_range*0.8*(_triangles[ii][2][1]-ymid)/max_range+0.5*(by-ay);
        
        _triangles[ii][0][2]=min_domain_range*0.8*(_triangles[ii][0][2]-zmid)/max_range+0.5*(bz-az);
        _triangles[ii][1][2]=min_domain_range*0.8*(_triangles[ii][1][2]-zmid)/max_range+0.5*(bz-az);
        _triangles[ii][2][2]=min_domain_range*0.8*(_triangles[ii][2][2]-zmid)/max_range+0.5*(bz-az);
    }
    
}


void custom_shape_3d::get_ver_tria(){
    ver_tria.clear();
    ver_tria.resize(vertices.size());
    ver_normals.resize(vertices.size());
    //i is the index of the triangle looking at: ith triangle
    for(int i=0; i<elements.size();i++){
        //jth vertex's index
        for(int j=0; j<3; j++){
            //the vertex looking at has index vti
            int vti=elements[i][j];
            //the vti'th vertex touch this triangle i
            ver_tria[vti].push_back(i);
        }
    }
    for(int i=0; i<vertices.size();i++){
        ver_normals[i].resize(3);
        ver_normals[i][0]=0.0;  ver_normals[i][1]=0.0; ver_normals[i][2]=0.0;
        int num_tria=(signed int)ver_tria[i].size();
        for(int j=0; j<num_tria; j++){
            int triai=ver_tria[i][j];
            ver_normals[i][0]+=normals[triai][0];
            ver_normals[i][1]+=normals[triai][1];
            ver_normals[i][2]+=normals[triai][2];
        }
        ver_normals[i][0]=1.0/num_tria*ver_normals[i][0];
        ver_normals[i][1]=1.0/num_tria*ver_normals[i][1];
        ver_normals[i][2]=1.0/num_tria*ver_normals[i][2];
    }
}

void custom_shape_3d::get_edge_tria(){
    edge_tria.clear();
    edge_tria.resize(num_edges);
    edge_normals.resize(num_edges);
    //i is the index of the triangle looking at: ith triangle
    for(int i=0; i<edges.size();i++){
        //jth edge's index
        for(int j=0;j<3;j++){
            //the edge looking at has index edi
            int edi=edges[i][j];
            //the edi'th edge touch this triangle i
            edge_tria[edi].push_back(i);
        }
    }
    for(int i=0; i<num_edges;i++){
        edge_normals[i].resize(3);
        edge_normals[i][0]=0.0;  edge_normals[i][1]=0.0; edge_normals[i][2]=0.0;
        int num_tria=(signed int)edge_tria[i].size();
        for(int j=0; j<num_tria; j++){
            int triai=edge_tria[i][j];
            edge_normals[i][0]+=normals[triai][0];
            edge_normals[i][1]+=normals[triai][1];
            edge_normals[i][2]+=normals[triai][2];
        }
        edge_normals[i][0]=1.0/num_tria*edge_normals[i][0];
        edge_normals[i][1]=1.0/num_tria*edge_normals[i][1];
        edge_normals[i][2]=1.0/num_tria*edge_normals[i][2];
    }
    
}

void custom_shape_3d::get_elements(){
    // construct indices of vertices for each triangle
    elements.clear();
    for (int i = 0;i < triangles.size();++i) {
        std::vector<int> element(3);
        element[0]=find_index(triangles[i][0][0],triangles[i][0][1],triangles[i][0][2]);
        element[1]=find_index(triangles[i][1][0],triangles[i][1][1],triangles[i][1][2]);
        element[2]=find_index(triangles[i][2][0],triangles[i][2][1],triangles[i][2][2]);
        elements.push_back(element);
    }
}

void custom_shape_3d::get_vertex() {
    // collect all vertices from triangle soup
    std::vector<std::vector<double> > vertex_list;
    vertex_list.clear();
    for (int i = 0;i < triangles.size();++i)
        for (int j = 0;j < 3;++j)
            vertex_list.push_back(triangles[i][j]);
    std::sort(vertex_list.begin(), vertex_list.end(), 
        [](const std::vector<double>& A, const std::vector<double>& B) -> bool {
            if (A[0] < B[0] - 1e-7) return true;
            if (A[0] > B[0] + 1e-7) return false;
            if (A[1] < B[1] - 1e-7) return true;
            if (A[1] > B[1] + 1e-7) return false;
            if (A[2] < B[2] - 1e-7) return true;
            return false;
        });

    // delete duplicated vertices
    vertices.clear();
    vertices.push_back(vertex_list[0]);
    for (int i = 1;i < (signed int)vertex_list.size();++i) {
        if (vertexCmp(vertex_list[i - 1][0],vertex_list[i - 1][1],vertex_list[i - 1][2], vertex_list[i][0],vertex_list[i][1],vertex_list[i][2])) {
            vertices.push_back(vertex_list[i]);
        }
    }
    
    double xmin=10000; double xmax=-10000;
    double ymin=10000; double ymax=-10000;
    double zmin=10000; double zmax=-10000;
    for(int i=0; i<vertices.size();i++){
        if(vertices[i][0]>xmax) xmax=vertices[i][0];
        if(vertices[i][0]<xmin) xmin=vertices[i][0];
        if(vertices[i][1]>ymax) ymax=vertices[i][1];
        if(vertices[i][1]<ymin) ymin=vertices[i][1];
        if(vertices[i][2]>zmax) zmax=vertices[i][2];
        if(vertices[i][2]<zmin) zmin=vertices[i][2];
    }
    printf("model x, y, z range: %g %g, %g %g, %g %g triangles: %d \n", xmin, xmax, ymin,ymax,zmin,zmax, (signed int)triangles.size());
    //center model to (0.5,0.5,0.5)
}

bool custom_shape_3d::vertexCmp(double x0,double y0,double z0,double x1,double y1,double z1) {
        if (x0 < x1 - 1e-7) return true;
        if (x0 > x1 + 1e-7) return false;
        if (y0 < y1 - 1e-7) return true;
        if (y0 > y1 + 1e-7) return false;
        if (z0 < z1 - 1e-7) return true;
        return false;
    }

void custom_shape_3d::CreateEdges() {
    std::map<std::pair<int, int>, int> _edge_index_mapping;
    num_edges = 0;
    _edge_index_mapping.clear();
    edges = elements;                 // _edges has same size as _elements
    for (int i = 0;i < elements.size();++i){
        for (int j = 0;j < 3;++j) {
            int idx0 = elements[i][j];
            int idx1 = elements[i][(j + 1) % 3];
            if (idx0 > idx1)
                std::swap(idx0, idx1);
            std::pair<int, int> edge = std::make_pair(idx0, idx1);
            std::map<std::pair<int, int>, int>::iterator it = _edge_index_mapping.find(edge);
            if (it == _edge_index_mapping.end()) {
                _edge_index_mapping[edge] = num_edges;
                edges[i][j] = num_edges ++;
            } else {
                edges[i][j] = it->second;
            }
        }
    }
}

int custom_shape_3d::find_index(double x, double y, double z) {
    int l = 0, r = vertices.size() - 1;
    for (;l < r;) {
        int mid = (l + r) >> 1;
        if (vertexCmp(vertices[mid][0],vertices[mid][1],vertices[mid][2], x, y, z))
            l = mid + 1;
        else
            r = mid;
    }
    return l;
}


//return the signed distance of a point (x,y,z) to the shape boundary
double custom_shape_3d::f_b_pts(double x, double y, double z) const{
    double d=1000000; //signed distance to be calculated
    int d_status;
    int d_tria; //track which triangle is the one with shortest distance: store tria index
    double p_lx; //point on mesh closest to pt, x component
    double p_ly; //point on mesh closest to pt, y component
    double p_lz; //point on mesh closest to pt, z component
    
    
    //the grid cell(pi,pj,pk) that the point is in
    int pi=(x-ax)/dx; if(x>=bx-1e-7&&x<=bx+1e-7){pi=nx-1;} 
    int pj=(y-ay)/dy; if(y>=by-1e-7&&y<=by+1e-7){pj=ny-1;}
    int pk=(z-az)/dz; if(z>=bz-1e-7&&z<=bz+1e-7){pk=nz-1;}
    int pip=pi; int pim=pi; int pjp=pj; int pjm=pj; int pkp=pk; int pkm=pk;
    
    std::vector<int> tria_list;
    std::vector<int> tria_list_unique;
    bool region_tria=false; 
    int cell_tria_ind=nx*ny*pk+nx*pj+pi;
    //if the cell the point is in this triangles
    if(cell_tria[cell_tria_ind].size()!=0){
        region_tria=true;
        //store the triangles index
        for(int triai=0;triai<cell_tria[nx*ny*pk+nx*pj+pi].size(); triai++){
            tria_list.push_back(cell_tria[nx*ny*pk+nx*pj+pi][triai]);
        }
    }
    //the increment layer to test to see if there are triangles
    int incm=0;
    while(region_tria==false){
        incm++;
        pip=pi+incm; pim=pi-incm; if(pip>=nx){pip=nx-1;} if(pim<0){pim=0;}
        pjp=pj+incm; pjm=pj-incm; if(pjp>=ny){pjp=ny-1;} if(pjm<0){pjm=0;}
        pkp=pk+incm; pkm=pk-incm; if(pkp>=nz){pkp=nz-1;} if(pkm<0){pkm=0;}
        for(int ii=pim;ii<=pip;ii++){
            for(int jj=pjm;jj<=pjp;jj++){
                for(int kk=pkm;kk<=pkp;kk++){
                    //only need to consider the outer new layer
                    if(ii==pim || ii==pip || jj==pjm || jj==pjp || kk==pkm || kk==pkp){
                        cell_tria_ind=nx*ny*kk+nx*jj+ii;
                        if(cell_tria[cell_tria_ind].size()!=0){
                            region_tria=true;
                            //store the triangles index
                            for(int triai=0;triai<cell_tria[cell_tria_ind].size();triai++){
                                tria_list.push_back(cell_tria[cell_tria_ind][triai]);
                            }
                        }
                    }
                }
            }
        }
    }
    //Now region_tria==true
    //get triangles of the bounding box of the previous region cube
    double fac0=incm+1.0;
    double aa=fac0*dx; double bb=fac0*dy; double cc=fac0*dz;
    double a=std::max(aa, std::max(bb,cc));
    double sphere_r=sqrt(3)*a;
    double cellx0=ax+pi*dx; double cellx1=cellx0+dx;
    double celly0=ay+pj*dy; double celly1=celly0+dy;
    double cellz0=az+pk*dz; double cellz1=cellz0+dz;
    
    int ih=(cellx1+sphere_r-ax)/dx; if(ih>=nx){ih=nx-1;}
    int il=(cellx0-sphere_r-ax)/dx; if(il<0){il=0;}
    int jh=(celly1+sphere_r-ay)/dy; if(jh>=ny){jh=ny-1;}
    int jl=(celly0-sphere_r-ay)/dy; if(jl<0){jl=0;}
    int kh=(cellz1+sphere_r-az)/dz; if(kh>=nz){kh=nz-1;}
    int kl=(cellz0-sphere_r-az)/dz; if(kl<0){kl=0;}
    
    for(int ii=il;ii<=ih;ii++){
        for(int jj=jl;jj<=jh;jj++){
            for(int kk=kl;kk<=kh;kk++){
                //only need to consider the outer new layer
                if(ii<pim || ii>pip || jj<pjm ||jj>pjp || kk<pkm || kk>pkp){
                    cell_tria_ind=nx*ny*kk+nx*jj+ii;
                    if(cell_tria[cell_tria_ind].size()!=0){
                        //store the triangles index
                        for(int triai=0;triai<cell_tria[cell_tria_ind].size();triai++){
                            tria_list.push_back(cell_tria[cell_tria_ind][triai]);
                        }
                    }
                }
            }
        }
    }

    //sort the tria_list and delete the duplicated ones, and get a unique list of triangles
    //default sort into an ascending array of tria id's
    std::sort(tria_list.begin(), tria_list.end());
    //delete duplicated ones 
    tria_list_unique.push_back(tria_list[0]);
    for(int i=0;i<tria_list.size();i++){
        if(tria_list_unique.back()!=tria_list[i]){
            tria_list_unique.push_back(tria_list[i]);
        }
    }
    tria_list.clear();
    
    //loop through the unique list of triangles 
    for(int i=0; i<tria_list_unique.size(); i++){
        int triai=tria_list_unique[i];
    //for(int i=0;i<triangles.size();i++){
        //int triai=i;
        double xA=triangles[triai][0][0];
        double yA=triangles[triai][0][1];
        double zA=triangles[triai][0][2];
        double xB=triangles[triai][1][0];
        double yB=triangles[triai][1][1];
        double zB=triangles[triai][1][2];
        double xC=triangles[triai][2][0];
        double yC=triangles[triai][2][1];
        double zC=triangles[triai][2][2];
        double a=normals[triai][0];
        double b=normals[triai][1];
        double c=normals[triai][2]; 
        double closestx, closesty, closestz;
        int status;
        double d_temp=f_tria_seg(x,y,z, xA,yA,zA, xB,yB,zB, xC,yC,zC, a,b,c,closestx, closesty, closestz,status);
        if(d_temp<d-1e-7){
            d_tria=triai; //which triangle
            d_status=status; //end points or in between or on edges
            d=d_temp; //shortest distance
            p_lx=closestx; //point on mesh closest to pt, x component
            p_ly=closesty; //point on mesh closest to pt, y component
            p_lz=closestz; //point on mesh closest to pt, z component
        }
    }
    tria_list_unique.clear();
    
    //decide the sign of the shortest distance: inside geometry, -; outside geometry, +;
    //the triangle looking at has index d_tria, with status d_status 0,1,2,3,4,5,6
    //if 6, directly use outward normal info of the triangle
    //if 0,1,2 (vertex), vertex normal is the average of the normals of triangle faces sharing the vertex
    //if 3,4,5 (edge), edge normal is the average of the two triagle faces sharing the edge
    
    if(d_status==6){ //closest point is inside the triangle
        //triangle looking at has index d_tria
        //in this case, use the triangle face normal directly
        double a=normals[d_tria][0];
        double b=normals[d_tria][1];
        double c=normals[d_tria][2]; 
        //p_l to p
        double plpx=x-p_lx;
        double plpy=y-p_ly;
        double plpz=z-p_lz;
        //plp dot n: if positive, the same side, outside geometry, +; 
        //if negative, opposite side, inside geometry, -
        double dot_temp=plpx*a+plpy*b+plpz*c;
        if (dot_temp<0-1e-7){
            d=-d;
        }
    }
    //if 0,1,2 (vertex), vertex normal is the average of the normals of triangle faces sharing the vertex
    else if(d_status==0 || d_status==1 || d_status==2){
        //get vertex index
        int veri=elements[d_tria][d_status];
        //weighted normal of the vertices
        double a=ver_normals[veri][0];
        double b=ver_normals[veri][1];
        double c=ver_normals[veri][2];
        //p_l to p
        double plpx=x-p_lx;
        double plpy=y-p_ly;
        double plpz=z-p_lz;
        //plp dot n: if positive, the same side, outside geometry, +; 
        //if negative, opposite side, inside geometry, -
        double dot_temp=plpx*a+plpy*b+plpz*c;
        if (dot_temp<0-1e-7){
            d=-d;
        }
    }
    //if 3,4,5 (edge), edge normal is the average of the two triagle faces sharing the edge
    else{
        //get edge index
        int edgei;
        if(d_status==3){edgei=edges[d_tria][0];}
        else if(d_status==4){edgei=edges[d_tria][1];}
        else{edgei=edges[d_tria][2];} //5
        //weighted normal of the edge
        double a=edge_normals[edgei][0];
        double b=edge_normals[edgei][1];
        double c=edge_normals[edgei][2];
        //p_l to p
        double plpx=x-p_lx;
        double plpy=y-p_ly;
        double plpz=z-p_lz;
        //plp dot n: if positive, the same side, outside geometry, +; 
        //if negative, opposite side, inside geometry, -
        double dot_temp=plpx*a+plpy*b+plpz*c;
        if (dot_temp<0-1e-7){
            d=-d;
        }
    }
    
    return d;
}





