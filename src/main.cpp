#include <iostream>
#include <algorithm>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

void loadB(MatrixXd& B, VectorXd& u, int l)
{
    for (int i = 0; i < B.rows(); i++)
    {
        if (i != l)
        {
            double coeficient = -1 * u(i) / u(l);
            for (int j = 0; j < B.cols(); j++)
                B(i, j) = B(i, j) + coeficient * B(l, j);
            u(i) = u(i) + coeficient * u(l);
        }
        else
        {
            for (int j = 0; j < B.cols(); j++)
                B(l, j) = B(l, j) / u(l);
            u(l) = u(l) / u(l);
        }
    }
}
void findC(VectorXd& c, double fo[], MatrixXd& vb, int n)
{
    // n = quant. de variaveis da base
    for (int i = 0; i < n; i++)
    {
        int index = static_cast<int>(vb(i, 0));
        c(i) = fo[index];
        
    }
}
void calculateP(VectorXd& c, MatrixXd& B, VectorXd& p)
{
    p = c.transpose() * B; //por padrao, c é coluna, então a transposicao é para transformá-lo em uma linha
    cout << p << endl;
}

// A[j] é uma coluna e nao uma linha
void calculateU(MatrixXd &B, MatrixXd& A, VectorXd& u, int j)
{
    u = B * A.col(j);
}

void calculateReducedC(double fo[], VectorXd& reduced_cost, MatrixXd& A, VectorXd& p, MatrixXd& vb){
    reduced_cost.setZero();
    for(int i = 0; i < A.cols(); i++){
        if(!(vb.array() == i).any()){
            //significa que é nao basica
            reduced_cost(i) = fo[i] - p.transpose() * A.col(i);
        }
    }
    cout << reduced_cost << endl;
}

bool isOptimal(VectorXd& reduced_cost){
    for(int i = 0; i < reduced_cost.rows(); i++){
        if(reduced_cost(i) > 0){
            return false;
        }
    }
    return true;
    //para escolher j: reduced_cost.maxCoeff();
}


int main()
{
    ////passando dados manualmente enquanto nao tenho a leitura de arquivos/////
    MatrixXd A(3, 7);
    A << 3, 2, 1, 2, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 4, 3, 3, 4, 0, 0, 1;
    MatrixXd B(3, 3);
    B << 3, 1, 0, 1, 1, 0, 4, 3, 1;
    
    double fo[] = {19, 13, 12, 17, 0, 0, 0}; // função objetivo
    VectorXd rhs(3);
    rhs << 225, 117, 420;

    MatrixXd variaveis_basicas(3, 2);
    variaveis_basicas.col(0) << 0, 2, 6;  //linha 1 da matriz serão os indices da matriz -> variavel -1
    variaveis_basicas.col(1) = B.inverse() * rhs;
    
    
    // calculando a inversa
    B = B.inverse();
    cout << "Matrix B inversa \n"
         << B << endl;

    VectorXd c(B.cols()); // c = coeficientes da fo das variaveis basicas
    VectorXd p(B.cols()); // p = duais
    VectorXd u(B.cols()); // u = B^-1*Aj
    VectorXd reduced_cost(A.cols());
    
    findC(c, fo, variaveis_basicas, A.rows());
    calculateP(c, B, p);
    calculateReducedC(fo, reduced_cost, A, p, variaveis_basicas);
    
    // matriz e u para teste de loadB (exemplo livro bertsimas)
    MatrixXd teste(3, 3);
    teste << 1, 2, 3, -2, 3, 1, 4, -3, -2;
    cout << "Here is the input matrix teste \n"
         << teste << endl;
    u << -4, 2, 2;
    loadB(teste, u, 2);
    cout << "Here is the input matrix teste \n"
         << teste << endl;
    
    return 0;
}
