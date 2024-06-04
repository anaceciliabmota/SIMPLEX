#include <iostream>
#include <algorithm>
#include <limits>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

//colocar exceção de u < 0

struct Solution {
    bool isOptimal;
    double z;
    MatrixXd variaveis_basicas;
} typedef Solution;

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
    //cout << p << endl;
}

// A[j] é uma coluna e nao uma linha
void calculateU(MatrixXd &B, MatrixXd& A, VectorXd& u, int j)
{
    u = B * A.col(j);
    //!any u > 0, nao tem solução otima
}

void calculateReducedC(double fo[], VectorXd& reduced_cost, MatrixXd& A, VectorXd& p, MatrixXd& vb){
    reduced_cost.setZero();
    for(int i = 0; i < A.cols(); i++){
        if(!(vb.col(0).array() == i).any()){
            //significa que é nao basica
            reduced_cost(i) = fo[i] - p.transpose() * A.col(i);
        }
    }
    //cout << reduced_cost << endl;
}

bool isOptimal(VectorXd& reduced_cost){
    for(int i = 0; i < reduced_cost.rows(); i++){
        if(reduced_cost(i) < 0){
            return false;
        }
    }
    return true;
    //para escolher j: reduced_cost.maxCoeff();
}

int findTeta(MatrixXd& vb, VectorXd u, double * teta){
    double min = numeric_limits<double>::infinity();
    int l;
    for(int i = 0; i < vb.rows();i++){
        if(u(i) > 0 && min > vb(i, 1)/u(i)){
            min = vb(i,1)/u(i);
            l = i; //ou l = vb(i, 0)?
        }
    }
    *teta = min;
    return l;
}
int chooseJ(VectorXd& reduced_cost){
    Index index;
    reduced_cost.minCoeff(&index);
    return static_cast<int>(index);
}

void changingVariables(MatrixXd& vb, VectorXd& u, double teta, int l, int j){
    for(int i = 0; i < vb.rows(); i++){
        if(i == l){
            vb(i, 0) = j;
            vb(i, 1) = teta;
        }else{
            vb(i, 1) = vb(i, 1) - teta*u(i);
        }
    }
}

Solution SIMPLEX(MatrixXd& A, MatrixXd& B, double fo[], MatrixXd& variaveis_basicas, int n){
    B = B.inverse();

    VectorXd c(n); // c = coeficientes da fo das variaveis basicas
    VectorXd p(n); // p = duais
    VectorXd u(n); // u = B^-1*Aj
    VectorXd reduced_cost(A.cols());
    
    Solution s;

    while(true){
        findC(c, fo, variaveis_basicas, A.rows());
        calculateP(c, B, p);
        
        calculateReducedC(fo, reduced_cost, A, p, variaveis_basicas);
        
        if(isOptimal(reduced_cost)){
            s.variaveis_basicas = variaveis_basicas;
            s.z = c.transpose() * variaveis_basicas.col(1);
            break;
        }
        else{
            int j = chooseJ(reduced_cost);
            calculateU(B, A, u , j);
            if(u.maxCoeff() <= 0){
                s.z = -1*numeric_limits<double>::infinity();
                break;
            }else{
                double teta;
                int l = findTeta(variaveis_basicas, u, &teta);
                changingVariables(variaveis_basicas, u, teta, l, j);
                loadB(B, u, l);
            }
        }
    }
    return s;
}
int main()
{
    ////passando dados manualmente enquanto nao tenho a leitura de arquivos/////
    int n;
    MatrixXd A(3, 7);
    n = A.rows();
    A << 3, 2, 1, 2, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 4, 3, 3, 4, 0, 0, 1;
    MatrixXd B(n, n);
    B << 3, 1, 0, 1, 1, 0, 4, 3, 1;
    
    double fo[] = {-19, -13, -12, -17, 0, 0, 0}; // função objetivo
    VectorXd rhs(n);
    rhs << 225, 117, 420;

    MatrixXd variaveis_basicas(n, 2);
    variaveis_basicas.col(0) << 0, 2, 6;  //linha 1 da matriz serão os indices da matriz -> variavel -1
    variaveis_basicas.col(1) = B.inverse() * rhs;


    Solution s = SIMPLEX(A, B, fo, variaveis_basicas, n);
    if(s.z != -1*numeric_limits<double>::infinity()){
        cout << "Solucao:" << endl; 
        for(int i = 0; i < n; i++){
            cout << "x" << variaveis_basicas(i, 0) + 1 << ": " << s.variaveis_basicas(i, 1) << endl;
        }
        cout << "z: " << s.z << endl;
    }else{
        cout << "Solução ótima é igual a menos infinito";
    }
    
    
    return 0;
}

/*MatrixXd A(, );
    n = A.rows();
    A << ;
    MatrixXd B(n, n);
    B << ;
    
    double fo[] = {-}; // função objetivo
    VectorXd rhs(n);
    rhs << 2;

    MatrixXd variaveis_basicas(n, 2);
    variaveis_basicas.col(0) << ;  //linha 1 da matriz serão os indices da matriz -> variavel -1
    variaveis_basicas.col(1) = B.inverse() * rhs;*/