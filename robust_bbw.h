#ifndef ROBUST_BBW_H
#define ROBUST_BBW_H
#include <Eigen/Core>
// Given a possibly artifact-ridden mesh (V,F) and a skeleton (C,BE), compute
// bounded biharmonic weights for each bone. The mesh will be cleaned,
// weights computed on the clean mesh and then mapped back to the original
// messy mesh.
//
// Inputs:
//   V  #V by 3 list of mesh vertex positions
//   F  #F by 3 list of mesh face indices into V
//   C  #C by 3 list of join positions
//   BE  #BE by 2 list of bone indices into C
// Output:
//   W  #V by #BE list of bone weights
// Returns true iff success
//
bool robust_bbw(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & C,
  const Eigen::MatrixXi & BE,
  Eigen::MatrixXd & W);
#endif
