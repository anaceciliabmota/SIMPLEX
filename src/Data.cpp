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

   /*fo = new double[qtVariables];
   //depois substituido por uma função que pega os valores do arquivo
   double funcao_o[] = {-3, -1, -1, 2, -1, 1, 1, -4};
   for(int i = 0; i < qtVariables; i++){
      fo[i] = funcao_o[i];
   }
   A << 1, 0, 3, 1, -5, -2, 4, -6,
        0, 1, -2, -1, 4, 1, -3, 5;
   rhs << 7, -3;   
   u << 8, 6, 4, 15, 2, 10, 10, 3;
   l << 0, 0, 0, 0, 0, 0, 0, 0;*/
}

/*
   double funcao_o[] = {-19, -13, -12, -17, 0, 0, 0};
   for(int i = 0; i < qtVariables; i++){
      fo[i] = funcao_o[i];
   }
   A << 3, 2, 1, 2, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 4, 3, 3, 4, 0, 0, 1;
   rhs << 225, 117, 420;   
   u << infinity, infinity, infinity, infinity, infinity, infinity, infinity;
   l << 0, 0, 0, 0, 0, 0, 0;

   */

