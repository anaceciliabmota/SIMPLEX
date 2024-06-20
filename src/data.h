#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

class Data{

public:
        Data(int qtConstraints, int qtVariables);
        ~Data();

        MatrixXd getMatrixA(){return A;}
        VectorXd getRHS(){return rhs;}  
        double getVectorU(int index){return u(index);}
        double getVectorL(int index){return l(index);}
        double getFO(int index){return fo[index];}

private:
        
        MatrixXd A;
        double * fo;
        VectorXd rhs;
        VectorXd u;
        VectorXd l;
  
};
