#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  
  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if(estimations.size() != ground_truth.size()
     || estimations.size() == 0){
    cout << "Invalid estimation or ground_truth data" << endl;
    return rmse;
  }
  
  //accumulate squared residuals
  for(unsigned int i=0; i < estimations.size(); ++i){
    
    VectorXd residual = estimations[i] - ground_truth[i];
    
    //coefficient-wise multiplication
    residual = residual.array()*residual.array();
    rmse += residual;
  }
  
  //calculate the mean
  rmse = rmse/estimations.size();
  
  //calculate the squared root
  rmse = rmse.array().sqrt();
  
  //return the result
  return rmse;
}

VectorXd Tools::PolarToCartesian(VectorXd x, int size) {
  float const rho = x(0);
  float const phi = x(1);
  float const rho_dot = x(2);
  float const px = rho * cos(phi);
  float const py = rho * sin(phi);
  float const vx = rho_dot * cos(phi);
  float const vy = rho_dot * sin(phi);
  float const v  = sqrt(vx * vx + vy * vy);

  VectorXd cartesian = VectorXd::Zero(size);
  cartesian(0) = px;
  cartesian(1) = py;
  cartesian(2) = v;

  return cartesian;
}

VectorXd Tools::NormalizeRadians(VectorXd x, int radIndex) {

  float rad = x(radIndex);
  
  while (rad < -M_PI) rad += 2 * M_PI;
  while (rad > M_PI) rad -= 2 * M_PI;

  VectorXd normalizedPolar = VectorXd::Zero(x.size());
  normalizedPolar << x;
  normalizedPolar(radIndex) = rad;
  return normalizedPolar;
}
