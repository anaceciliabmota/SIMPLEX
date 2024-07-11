#include "Data.h"
#include "mpsReader.h"
#include <iostream>
#include <Eigen/Dense>
#include <limits>

#define infinity std::numeric_limits<double>::infinity()

using namespace std;
using namespace Eigen;
   Data::Data(MatrixXd& mps_A, VectorXd& mps_b, VectorXd& mps_c, VectorXd& mps_ub, VectorXd& mps_lb)
   {
      A = mps_A;
      rhs = mps_b;
      fo = mps_c;
      u = mps_ub;
      l = mps_lb;
   }


