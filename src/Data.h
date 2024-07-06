
#ifndef DATA_H
#define DATA_H

#include <iostream>
#include <Eigen/Dense>
#include "mpsReader.h"

using namespace std;
using namespace Eigen;

class Data{

public:
        Data(MatrixXd& mps_A, VectorXd& mps_b, VectorXd& mps_c, VectorXd& mps_ul, VectorXd& mps_lb);
        //~Data();

        MatrixXd getMatrixA(){return A;}
        VectorXd getRHS(){return rhs;}  
        VectorXd getVectorU(){return u;}
        VectorXd getVectorL(){return l;}
        VectorXd getFO(){return fo;}
        void setMatrixA(MatrixXd& newA){A = newA;}
        void setVectorU(VectorXd& newU){u = newU;}
        void setVectorL(VectorXd& newL){l = newL;}
        void setFO(VectorXd& newFO){fo = newFO;}
        
private:
        MatrixXd A;
        VectorXd fo;
        VectorXd rhs;
        VectorXd u;
        VectorXd l;
  
};
#endif