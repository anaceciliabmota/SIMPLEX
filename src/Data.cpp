#include "data.h"
#include <iostream>
#include <Eigen/Dense>
#include <limits>

#define infinity std::numeric_limits<double>::infinity()

using namespace std;
using namespace Eigen;
Data::Data(int qtConstraints, int qtVariables)
    : A(qtConstraints, qtVariables), 
      rhs(qtConstraints), 
      u(qtVariables), 
      l(qtVariables)
{
   fo = new double[qtVariables];
   //depois substituido por uma função que pega os valores do arquivo
   double funcao_o[] = {-19, -13, -12, -17, 0, 0, 0};
   for(int i = 0; i < qtVariables; i++){
      fo[i] = funcao_o[i];
   }
   A << 3, 2, 1, 2, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 4, 3, 3, 4, 0, 0, 1;
   rhs << 225, 117, 420;   
   u << infinity, infinity, infinity, infinity, infinity, infinity, infinity;
   l << 0, 0, 0, 0, 0, 0, 0;
}