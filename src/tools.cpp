#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  assert(estimations.size()==ground_truth.size());
  assert(estimations.size()>0);
  VectorXd rmse = VectorXd::Zero(4);
  for (int i=0; i<estimations.size(); i++) {
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();
    rmse += residual;
  }
  rmse /= estimations.size();
//  std::cout<<"rmse = "<<rmse<<std::endl;
  return rmse;
}
