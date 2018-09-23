#ifndef CLEAN_H
#define CLEAN_H
#include <Eigen/Core>
// Given a mesh (V,F), fix self-intersections and open boundaries
//
// Inputs:
//   V  #V by 3 list of mesh vertex positions
//   F  #F by 3 list of mesh face indices into V
// Outputs:
//   CV  #CV by 3 list of mesh vertex positions
//   CF  #CF by 3 list of mesh face indices into CF
//   IM  >=#CV list of indices into V, so that (CV,IM(F)) produces the same
//     mesh as (V,F)
// Returns true iff success.
//
bool clean(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & CV,
  Eigen::MatrixXi & CF,
  Eigen::VectorXi & IM);
#endif
