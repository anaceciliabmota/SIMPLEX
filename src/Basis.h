#ifndef BASIS_H
#define BASIS_H

#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/Sparse>
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
    void getDi(VectorXd& d, VectorXd&a, int cont);
    void getPi(VectorXd& pi, VectorXd& c, int cont);
    void getD(VectorXd& d, VectorXd& a);
    void getP(VectorXd& p, VectorXd& c);
    void addElement(pair<int, VectorXd>& p);
    void loadB(Data& data, MatrixXd &variaveis_basicas);
   
    vector<pair<int, VectorXd>> v;
private:
    void * Symbolic;
    void * Numeric;
    double * null;
    //MatrixXd B;
    SparseMatrix<double> B;

};

#endif // BASIS_H
