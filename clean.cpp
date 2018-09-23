#include "clean.h"
#include <igl/barycenter.h>
#include <igl/boundary_facets.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/copyleft/tetgen/cdt.h>
#include <igl/winding_number.h>
#include <igl/unique_simplices.h>
#include <igl/remove_unreferenced.h>
#include <igl/writeOBJ.h>
#include <igl/writeMESH.h>
#include <igl/copyleft/cgal/remesh_self_intersections.h>
#include <igl/REDRUM.h>
#include <iostream>

bool clean(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & CV,
  Eigen::MatrixXi & CF,
  Eigen::VectorXi & IM)
{
  using namespace igl;
  using namespace igl::copyleft::tetgen;
  using namespace igl::copyleft::cgal;
  using namespace Eigen;
  using namespace std;
  //writeOBJ("VF.obj",V,F);
  const auto & validate_IM = [](
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXd & CV,
    const Eigen::VectorXi & IM)
  {
    assert(IM.size() >= CV.rows());
    for(int i = 0;i<CV.rows();i++)
    {
      if(IM(i)<V.rows() && IM(i)>= 0)
      {
        double diff = (V.row(IM(i))-CV.row(i)).norm();
        if(diff>1e-6)
        {
          cout<<i<<": "<<IM(i)<<" "<<diff<<endl;
        }
      }
    }
  };
  {
    MatrixXi _1;
    VectorXi _2;
    cout<<"clean: remesh_self_intersections"<<endl;
    remesh_self_intersections(V,F,{false,false,false},CV,CF,_1,_2,IM);
    for_each(CF.data(),CF.data()+CF.size(),[&IM](int & a){a=IM(a);});
    //validate_IM(V,CV,IM);
    cout<<"clean: remove_unreferenced"<<endl;
    {
      MatrixXi oldCF = CF;
      unique_simplices(oldCF,CF);
    }
    MatrixXd oldCV = CV;
    MatrixXi oldCF = CF;
    VectorXi nIM;
    remove_unreferenced(oldCV,oldCF,CV,CF,nIM);
    // reindex nIM through IM
    for_each(IM.data(),IM.data()+IM.size(),[&nIM](int & a){a=a>=0?nIM(a):a;});
    //validate_IM(V,CV,IM);
  }
  MatrixXd TV;
  MatrixXi TT;
  {
    MatrixXi _1;
    // c  convex hull
    // Y  no boundary steiners
    // p  polygon input
    // T1e-16  sometimes helps tetgen
    cout<<"clean: tetrahedralize"<<endl;
    writeOBJ("CVCF.obj",CV,CF);
    CDTParam params;
    params.flags = "CYT1e-16";
    params.use_bounding_box = true;
    if(cdt(CV,CF,params,TV,TT,_1) != 0)
    {
      cout<<REDRUM("CDT failed.")<<endl;
      return false;
    }
    //writeMESH("TVTT.mesh",TV,TT,MatrixXi());
  }
  {
    MatrixXd BC;
    barycenter(TV,TT,BC);
    VectorXd W;
    cout<<"clean: winding_number"<<endl;
    winding_number(V,F,BC,W);
    W = W.array().abs();
    const double thresh = 0.5;
    const int count = (W.array()>thresh).cast<int>().sum();
    MatrixXi CT(count,TT.cols());
    int c = 0;
    for(int t = 0;t<TT.rows();t++)
    {
      if(W(t)>thresh)
      {
        CT.row(c++) = TT.row(t);
      }
    }
    assert(c==count);
    boundary_facets(CT,CF);
    //writeMESH("CVCTCF.mesh",TV,CT,CF);
    cout<<"clean: remove_unreferenced"<<endl;
    // Force all original vertices to be referenced
    MatrixXi FF = F;
    for_each(FF.data(),FF.data()+FF.size(),[&IM](int & a){a=IM(a);});
    int ncf = CF.rows();
    MatrixXi ref(ncf+FF.rows(),3);
    ref<<CF,FF;
    VectorXi nIM;
    remove_unreferenced(TV,ref,CV,CF,nIM);
    // Only keep boundary faces
    CF.conservativeResize(ncf,3);
    cout<<"clean: IM.minCoeff(): "<<IM.minCoeff()<<endl;
    // reindex nIM through IM
    for_each(IM.data(),IM.data()+IM.size(),[&nIM](int & a){a=a>=0?nIM(a):a;});
    cout<<"clean: IM.minCoeff(): "<<IM.minCoeff()<<endl;
    //validate_IM(V,CV,IM);
  }
  return true;
}

