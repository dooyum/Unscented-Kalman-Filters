#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);
  
  /**
   * Converts polar coordinates to cartesian coordinates.
   * @param x polar coordinates
   */
  VectorXd PolarToCartesian(VectorXd x, int size);
  
  /**
   * Normalizes the phi value of polar coordinates to be between pi and -pi.
   * @param x polar coordinates
   */
  VectorXd NormalizeRadians(VectorXd x, int radIndex);
};

#endif /* TOOLS_H_ */
