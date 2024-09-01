#ifndef BASIS_H
#define BASIS_H

#include "../Eigen/Dense"
#include "../Eigen/LU"
#include "../Eigen/Sparse"
#include <vector>
#include <utility>
#include "Data.h"

//#include "Eigen/src/Core/Matrix.h"
//#include <suitesparse/umfpack.h> 

using namespace Eigen;
using namespace std;

class Basis {
public:
    Basis(MatrixXd& B);
    
    ~Basis();   
    void getD(VectorXd& d, VectorXd& a, int cont, MatrixXd variaveis_basicas, Data& data);
    void getP(VectorXd& p, VectorXd& c, int cont);
    void addElement(pair<int, VectorXd>& p);
    void loadB(Data& data, MatrixXd &variaveis_basicas);
   
    vector<pair<int, VectorXd>> v;
private:
    void * Symbolic;
    void * Numeric;
    double * null;
    SparseMatrix<double> B;

};

#endif // BASIS_H