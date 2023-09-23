#include "adf_2d.hh"

using namespace voro;

/**
 * @brief Constructor for adf_2d class.
 * @param adf_stat Pointer to the adf_stat_2d object.
 * @param shp_ Pointer to the shape_2d object.
 * @param x0_ Lower x-coordinate of the quad cell.
 * @param x1_ Upper x-coordinate of the quad cell.
 * @param y0_ Lower y-coordinate of the quad cell.
 * @param y1_ Upper y-coordinate of the quad cell.
 * @param err_tol_ Error tolerance of the subdivision rule.
 * @param max_level_ Maximum level of the quad tree allowed.
 * @param level_ Level of the quad tree.
 */
adf_2d::adf_2d(adf_stat_2d *adf_stat, shape_2d *shp_, 
    double x0_, double x1_, double y0_, double y1_, 
    double err_tol_, int max_level_, int level_)
:shp(shp_) 
{
    adf_stat->cell_ct++;
    
    child = false;
    x0=x0_; x1=x1_; y0=y0_; y1=y1_;
    ul=NULL;
    ur=NULL;
    ll=NULL;
    lr=NULL;
    level=level_;
    err_tol=err_tol_;
    max_level=max_level_;
    v0=shp->sdf(x0,y0); v1=shp->sdf(x1,y0); v2=shp->sdf(x1,y1); v3=shp->sdf(x0,y1);
    get_v45678(); //get v4, v5, v6, v7, v8
    get_biv45678(); //get bilinear values of v4, v5, v6, v7, v8
    if(subdivide()){
        adf_stat->cell_ct--;
        child=true;
        ul=new adf_2d(adf_stat, shp_, x0,0.5*(x0+x1),0.5*(y0+y1),y1,v5,v6,v4,v3,err_tol,max_level,level+1);
        ur=new adf_2d(adf_stat, shp_, 0.5*(x0+x1),x1,0.5*(y0+y1),y1,v6,v7,v2,v4,err_tol,max_level,level+1);
        ll=new adf_2d(adf_stat, shp_, x0,0.5*(x0+x1),y0,0.5*(y0+y1),v0,v8,v6,v5,err_tol,max_level,level+1);
        lr=new adf_2d(adf_stat, shp_, 0.5*(x0+x1),x1,y0,0.5*(y0+y1),v8,v1,v7,v6,err_tol,max_level,level+1);
        if(level+1>adf_stat->depth){adf_stat->depth=level+1;}
    }

    get_err(adf_stat);
}

/**
 * @brief Constructor for adf_2d class.
 * @param adf_stat Pointer to the adf_stat_2d object.
 * @param shp_ Pointer to the shape_2d object.
 * @param x0_ Lower x-coordinate of the quad cell.
 * @param x1_ Upper x-coordinate of the quad cell.
 * @param y0_ Lower y-coordinate of the quad cell.
 * @param y1_ Upper y-coordinate of the quad cell.
 * @param v0_ Value at the bottom-left corner of the quad cell.
 * @param v1_ Value at the bottom-right corner of the quad cell.
 * @param v2_ Value at the top-right corner of the quad cell.
 * @param v3_ Value at the top-left corner of the quad cell.
 * @param err_tol_ Error tolerance of the subdivision rule.
 * @param max_level_ Maximum level of the quad tree allowed.
 * @param level_ Level of the quad tree.
 */
adf_2d::adf_2d(adf_stat_2d *adf_stat, shape_2d *shp_, 
    double x0_, double x1_, double y0_, double y1_, 
    double v0_, double v1_, double v2_, double v3_, 
    double err_tol_, int max_level_, int level_)
:shp(shp_) 
{
    adf_stat->cell_ct++;
    
    child = false;
    x0=x0_; x1=x1_; y0=y0_; y1=y1_;
    ul=NULL;
    ur=NULL;
    ll=NULL;
    lr=NULL;
    level=level_;
    err_tol=err_tol_;
    max_level=max_level_;
    v0=v0_; v1=v1_; v2=v2_; v3=v3_;
    get_v45678(); //get v4, v5, v6, v7, v8
    get_biv45678(); //get bilinear values of v4, v5, v6, v7, v8
    if(subdivide()){
        adf_stat->cell_ct--;
        child=true;
        ul=new adf_2d(adf_stat, shp_, x0,0.5*(x0+x1),0.5*(y0+y1),y1,v5,v6,v4,v3,err_tol,max_level,level+1);
        ur=new adf_2d(adf_stat, shp_, 0.5*(x0+x1),x1,0.5*(y0+y1),y1,v6,v7,v2,v4,err_tol,max_level,level+1);
        ll=new adf_2d(adf_stat, shp_, x0,0.5*(x0+x1),y0,0.5*(y0+y1),v0,v8,v6,v5,err_tol,max_level,level+1);
        lr=new adf_2d(adf_stat, shp_, 0.5*(x0+x1),x1,y0,0.5*(y0+y1),v8,v1,v7,v6,err_tol,max_level,level+1);
        if(level+1>adf_stat->depth){adf_stat->depth=level+1;}
    }

}

/**
 * @brief Calculate and update the maximum error in the quad tree.
 * @param adf_stat Pointer to the adf_stat_2d object.
 */
void adf_2d::get_err(adf_stat_2d *adf_stat){
    if(child==false){
        double v4err=abs(v4-biv4);
        double v5err=abs(v5-biv5);
        double v6err=abs(v6-biv6);
        double v7err=abs(v7-biv7);
        double v8err=abs(v8-biv8);
        double maxverr=std::max(std::max(std::max(std::max(v4err, v5err), v6err), v7err),v8err);
        if(maxverr>adf_stat->err){adf_stat->err=maxverr;}
    }
    else{
        ul->get_err(adf_stat);
        ur->get_err(adf_stat);
        ll->get_err(adf_stat);
        lr->get_err(adf_stat);
    }
}

/**
 * @brief Calculate and update the values v4, v5, v6, v7, v8 based on the shape's signed distance.
 */
void adf_2d::get_v45678(){
    double midx=0.5*(x0+x1);
    double midy=0.5*(y0+y1);
    v4=shp->sdf(midx,y1);
    v5=shp->sdf(x0,midy);
    v6=shp->sdf(midx,midy);
    v7=shp->sdf(x1,midy);
    v8=shp->sdf(midx,y0);
}

/**
 * @brief Print the quad tree to a text file.
 * 
 * @param fp The file path of the output text file.
 */
void adf_2d::print_tree(const char *fp){
    FILE *outFile = fopen(fp, "a");
    print_tree_func(outFile);
    fclose(outFile);
}

/**
 * @brief Recursive function to print the quad tree to a text file.
 * 
 * @param[out] outFile Pointer to the output file.
 */
void adf_2d::print_tree_func(FILE *&outFile){
    if(child==false){
        fprintf(outFile, "%g %g\n%g %g\n\n", x0, y0, x1, y0);
        fprintf(outFile, "%g %g\n%g %g\n\n", x1, y0, x1, y1);
        fprintf(outFile, "%g %g\n%g %g\n\n", x1, y1, x0, y1);
        fprintf(outFile, "%g %g\n%g %g\n\n", x0, y1, x0, y0);
    }
    else{
        ul->print_tree_func(outFile);
        ur->print_tree_func(outFile);
        ll->print_tree_func(outFile);
        lr->print_tree_func(outFile);
    }
}

/**
 * @brief Check if a point is within the quad cell.
 * 
 * @param x The x-coordinate of the point.
 * @param y The y-coordinate of the point.
 * @return True if the point is within the quad cell, false otherwise.
 */
bool adf_2d::point_in_cell(double x, double y){
    if(x>=x0 && x<=x1 && y>=y0 && y<=y1){
        return true;
    }
    else{
        return false;
    }
}

/**
 * @brief Get the value of the point within the quad cell using bilinear interpolation.
 * 
 * @param x The x-coordinate of the point.
 * @param y The y-coordinate of the point.
 * @return The interpolated value of the point within the quad cell.
 * @throws std::invalid_argument if the point is not within the cell's domain range.
 */
double adf_2d::getVal(double x, double y){
    if(point_in_cell(x,y)){
        if(child==false){
            return bilinear_i(x, y);
        }
        else{ //child == true
            if(ul->point_in_cell(x,y)){
                return ul->getVal(x,y);
            }
            else if(ur->point_in_cell(x,y)){
                return ur->getVal(x,y);
            }
            else if(ll->point_in_cell(x,y)){
                return ll->getVal(x,y);
            }
            else{ //lr->point_in_cell(x,y)
                return lr->getVal(x,y);
            }
        }
    }
    else{
        throw std::invalid_argument( "point (x, y) is not in cell domain range" );
    }
}

/**
 * @brief Calculate the bilinear interpolated values (biv4, biv5, biv6, biv7, biv8) within the quad cell.
 */
void adf_2d::get_biv45678(){
    double midx=0.5*(x0+x1);
    double midy=0.5*(y0+y1);
    biv4=bilinear_i(midx,y1);
    biv5=bilinear_i(x0,midy);
    biv6=bilinear_i(midx,midy);
    biv7=bilinear_i(x1,midy);
    biv8=bilinear_i(midx,y0);
}

/**
 * @brief Determine whether to subdivide the quad cell based on the maximum error of its values at the midpoints of the edges and at the cell center.
 * @return True if the quad cell should be subdivided, false otherwise.
 */
bool adf_2d::subdivide(){
    bool sub=false;
    if(level<max_level){
        double v4err=abs(v4-biv4);
        double v5err=abs(v5-biv5);
        double v6err=abs(v6-biv6);
        double v7err=abs(v7-biv7);
        double v8err=abs(v8-biv8);
        double maxverr=std::max(std::max(std::max(std::max(v4err, v5err), v6err), v7err),v8err);
        if(maxverr>err_tol){
            sub=true;
        }
    }
    return sub;
}

/**
 * @brief Calculate the transformed epsilon value for a given x-coordinate. Used for bilinear interpolation.
 * @param x x-coordinate.
 * @return The epsilon value.
 */
double adf_2d::epsilon(double x){
    double epsilon_=-1.0+2.0*(x-x0)/(x1-x0);
    return epsilon_;
}

/**
 * @brief Calculate the transformed eta value for a given y-coordinate. Used for bilinear interpolation.
 * @param y y-coordinate.
 * @return The eta value.
 */
double adf_2d::eta(double y){
    double eta_=-1.0+2.0*(y-y0)/(y1-y0);
    return eta_;
}

/**
 * @brief Perform bilinear interpolation within a cell and return the interpolated value at the specified point (x, y).
 * @param x The x-coordinate of the point.
 * @param y The y-coordinate of the point.
 * @return The interpolated value at the specified point.
 */
double adf_2d::bilinear_i(double x, double y){
    double epsilon_i=epsilon(x);
    double eta_i=eta(y);
    double m_1_epsilon_i=1.0-epsilon_i;
    double p_1_epsilon_i=1.0+epsilon_i;
    double m_1_eta_i=1.0-eta_i;
    double p_1_eta_i=1.0+eta_i;
    double vphi0=v0*(m_1_epsilon_i)*(m_1_eta_i);
    double vphi1=v1*(p_1_epsilon_i)*(m_1_eta_i);
    double vphi2=v2*(p_1_epsilon_i)*(p_1_eta_i);
    double vphi3=v3*(m_1_epsilon_i)*(p_1_eta_i);
    double v_sum=0.25*(vphi0+vphi1+vphi2+vphi3);
    return v_sum;
}

