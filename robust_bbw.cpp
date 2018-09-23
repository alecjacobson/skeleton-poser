#include "robust_bbw.h"
#include "clean.h"

#include <igl/writeDMAT.h>
#include <igl/writeOBJ.h>
#include <igl/writeMESH.h>
#include <igl/copyleft/tetgen/mesh_with_skeleton.h>
#include <igl/boundary_conditions.h>
#include <igl/normalize_row_sums.h>
#include <igl/volume.h>
#include <igl/bbw.h>
#include <igl/REDRUM.h>
#ifdef WITH_MATLAB
#include <igl/matlab/MatlabWorkspace.h>
#endif
#include <iostream>

bool robust_bbw(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & C,
  const Eigen::MatrixXi & BE,
  Eigen::MatrixXd & W)
{
  using namespace igl;
  using namespace igl::copyleft::tetgen;
  using namespace Eigen;
  using namespace std;
  // clean mesh
  MatrixXd CV;
  MatrixXi CF;
  VectorXi IM;
  if(!clean(V,F,CV,CF,IM))
  {
    return false;
  }

#ifdef WITH_MATLAB
  igl::matlab::MatlabWorkspace mw;
  mw.save(V,"V");
  mw.save_index(F,"F");
  mw.save(CV,"CV");
  mw.save_index(CF,"CF");
  mw.save_index(IM,"IM");
#endif

  MatrixXd TV;
  MatrixXi TT;
  // compute tet-mesh
  {
    MatrixXi _1;
#ifdef VERBOSE
    cerr<<"mesh_with_skeleton"<<endl;
    //writeOBJ("mesh_with_skeleton.obj",CV,CF);
#endif
    if(!mesh_with_skeleton(CV,CF,C,{},BE,{},10,"pq1.5Y",TV,TT,_1))
    {
      cout<<REDRUM("tetgen failed.")<<endl;
      return false;
    }
    //writeMESH("mesh_with_skeleton.mesh",TV,TT,MatrixXi());
  }
  cout<<"TV.block ~= V? "<<TV.block(0,0,V.rows(),V.cols()).isApprox(V)<<endl;
  cout<<"TV.block ~= CV? "<<TV.block(0,0,CV.rows(),CV.cols()).isApprox(CV)<<endl;
  // Finally, tetgen may have still included some insanely small tets.
  // Just ignore these during weight computation (and hope they don't isolate
  // any vertices).
  {
    const MatrixXi oldTT = TT;
    VectorXd vol;
    volume(TV,TT,vol);
    const double thresh = 1e-17;
    const int count = (vol.array()>thresh).cast<int>().sum();
    TT.resize(count,oldTT.cols());
    int c = 0;
    for(int t = 0;t<oldTT.rows();t++)
    {
      if(vol(t)>thresh)
      {
        TT.row(c++) = oldTT.row(t);
      }
    }
  }

  // compute weights
  VectorXi b;
  MatrixXd bc;
  if(!boundary_conditions(TV,TT,C,{},BE,{},b,bc))
  {
    cout<<REDRUM("boundary_conditions failed.")<<endl;
    return false;
  }
  // compute BBW
  // Default bbw data and flags
  BBWData bbw_data;
  bbw_data.verbosity = 1;
  bbw_data.active_set_params.max_iter = 10;
#ifdef VERBOSE
  cerr<<"bbw"<<endl;
#endif
  // Weights matrix
  Eigen::MatrixXd TW;
  if(!igl::bbw(TV,TT,b,bc,bbw_data,TW))
  {
    return false;
  }
  // Normalize weights to sum to one
  normalize_row_sums(TW,TW);

#ifdef WITH_MATLAB
  mw.save(TV,"TV");
  mw.save_index(TT,"TT");
  mw.save(TW,"TW");
#endif

  W.resize(V.rows(),TW.cols());
  for(int i = 0;i<V.rows();i++)
  {
    W.row(i) = TW.row(IM(i));
  }

#ifdef WITH_MATLAB
  mw.save(W,"W");
  mw.write("robust_bbw.mat");
#endif

  return true;
}
