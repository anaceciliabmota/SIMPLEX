#include <iostream>
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
void findC(VectorXd& c, double fo[], int vb[], int n)
{
    // n = quant. de variaveis da base
    for (int i = 0; i < n; i++)
    {
        c(i) = fo[vb[i]];
        
    }
}
void calculateP(VectorXd& c, MatrixXd& B, VectorXd& p)
{
    p = c.transpose() * B; //por padrao, c é coluna, então a transposicao é para transformá-lo em uma linha
    
}
// A[j] é uma coluna e nao uma linha

void calculateU(MatrixXd &B, MatrixXd& A, VectorXd& u, int j)
{
    
    u = B * A.col(j);
}

int main()
{
    ////passando dados manualmente enquanto nao tenho a leitura de arquivos/////
    MatrixXd A(3, 7);
    A << 3, 2, 1, 2, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 4, 3, 3, 4, 0, 0, 1;
    MatrixXd B(3, 3);
    B << 3, 1, 0, 1, 1, 0, 4, 3, 1;
    int variaveis_basicas[] = {0, 2, 6};     // indices da matriz -> variavel -1
    double fo[] = {19, 13, 12, 17, 0, 0, 0}; // função objetivo

    // calculando a inversa
    B = B.inverse();
    cout << "Matrix B inversa \n"
         << B << endl;

    VectorXd c(3); // c = coeficientes da fo das variaveis basicas
    VectorXd p(3); // p = duais
    VectorXd u(3); // u = B^-1*Aj
    
    findC(c, fo, variaveis_basicas, B.rows()); // achar c a partir de quais variaveis tao na base
    cout << c << endl;
    calculateP(c, B, p);
     cout << p << endl;
    
    
    // matriz e u para teste de loadB
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
