#include "Basis.h"
#include <iostream>

#include <string>
#include <limits>
#include <numeric>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <suitesparse/umfpack.h>

using namespace std;
using namespace Eigen;

Basis::Basis(MatrixXd &B)
{
    null = nullptr;
    this->B = B.sparseView();

    (void)umfpack_di_symbolic(this->B.rows(), this->B.cols(), this->B.outerIndexPtr(), this->B.innerIndexPtr(), this->B.valuePtr(), &Symbolic, null, null);
    (void)umfpack_di_numeric(this->B.outerIndexPtr(), this->B.innerIndexPtr(), this->B.valuePtr(), Symbolic, &Numeric, null, null);
}

Basis::~Basis()
{
    // Cleanup UMFPACK
    umfpack_di_free_numeric(&Numeric);
    umfpack_di_free_symbolic(&Symbolic);
}

void Basis::getP(VectorXd &p, VectorXd &c)
{
    VectorXd v_1(c.rows());
    VectorXd v_i = c;

    for (int i = this->v.size() - 1; i >= 0; i--)
    {
        int indice = this->v[i].first;
        double soma = 0;
        for (int j = 0; j < v_i.size(); j++)
        {
            if (j != indice)
            {
                v_1(j) = v_i(j);
                soma += this->v[i].second(j) * v_1(j);
            }
        }
        v_1(indice) = (v_i(indice) - soma) / this->v[i].second(indice);
        v_i = v_1;
    }
    // p = this->B.transpose().colPivHouseholderQr().solve(v_i);

    (void)umfpack_di_solve(UMFPACK_A, this->B.outerIndexPtr(), this->B.innerIndexPtr(), this->B.valuePtr(), p.data(), v_i.data(), Numeric, null, null);
}

void Basis::getD(VectorXd &d, VectorXd &a)
{
    // VectorXd v_1 = this->B.colPivHouseholderQr().solve(a);
    VectorXd v_1(a.rows());

    (void)umfpack_di_solve(UMFPACK_A, this->B.outerIndexPtr(), this->B.innerIndexPtr(), this->B.valuePtr(), v_1.data(), a.data(), Numeric, null, null);

    for (int i = 0; i < this->v.size(); i++)
    {
        VectorXd v_i(v_1.size());
        int indice = this->v[i].first;
        v_i(indice) = v_1(indice) / this->v[i].second(indice);
        for (int j = 0; j < v_1.size(); j++)
        {
            if (j != indice)
            {
                v_i(j) = v_1(j) - this->v[i].second(j) * v_i(indice);
            }
        }
        v_1 = v_i;
    }
    d = v_1;
}

void Basis::addElement(pair<int, VectorXd> &p)
{
    v.push_back(p);
}

void Basis::loadB(Data &data, MatrixXd &variaveis_basicas)
{
    MatrixXd Bb(data.getMatrixA()->rows(), data.getMatrixA()->rows());
    for (int i = 0; i < variaveis_basicas.rows(); i++)
    {
        Bb.col(i) = data.getMatrixA()->col(variaveis_basicas(i, 0));
    }
    this->B = Bb.sparseView();
    this->v.clear();

    umfpack_di_free_symbolic(&Symbolic);
    umfpack_di_free_numeric(&Numeric);

    (void)umfpack_di_symbolic(this->B.rows(), this->B.cols(), this->B.outerIndexPtr(), this->B.innerIndexPtr(), this->B.valuePtr(), &Symbolic, null, null);
    (void)umfpack_di_numeric(this->B.outerIndexPtr(), this->B.innerIndexPtr(), this->B.valuePtr(), Symbolic, &Numeric, null, null);
}
