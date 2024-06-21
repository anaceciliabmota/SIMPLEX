#include "data.h"
#include <iostream>
#include <algorithm>
#include <limits>
#include <numbers>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#define EPSILON 10e-5

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
void findC(VectorXd& c, MatrixXd& vb, Data* data, int n)
{
    // n = quant. de variaveis da base
    for (int i = 0; i < n; i++)
    {
        int index = static_cast<int>(vb(i, 0));
        c(i) = data->getFO(index);
        
    }
}
void calculateP(VectorXd& c, MatrixXd& B, VectorXd& p)
{
    p = c.transpose() * B; //por padrao, c é coluna, então a transposicao é para transformá-lo em uma linha
    //cout << p << endl;
}

// A[j] é uma coluna e nao uma linha
void calculateU(MatrixXd &B, Data* data, VectorXd& u, int j)
{
    u = B * data->getMatrixA().col(j);
    //!any u > 0, nao tem solução otima
}

void calculateReducedC( VectorXd& reduced_cost, VectorXd& p, MatrixXd& vb, Data* data){
    reduced_cost.setZero();
    for(int i = 0; i < data->getMatrixA().cols(); i++){
        if(!(vb.col(0).array() == i).any()){
            //significa que é nao basica
            reduced_cost(i) = data->getFO(i) - p.transpose() * data->getMatrixA().col(i);
        }
    }
    //cout << reduced_cost << endl;
}

bool isOptimal(VectorXd& reduced_cost){
    for(int i = 0; i < reduced_cost.rows(); i++){
        if(reduced_cost(i) < -EPSILON){
            return false;
        }
    }
    return true;
}
bool isSatisfied(Data * data, int j, double teta, int i, MatrixXd& vb, bool islower, VectorXd& u){
    //verifica se xB(t) é satisfeita para todos os upperbound e lowerbound
    for(int k = 0; k < vb.rows(); k++){
        int index = static_cast<int>(vb(k, 0));
        if(islower){
            if(data->getVectorL(index) > vb(k, 1) - teta*u(i)/*xB*/ || data->getVectorU(index) < vb(k, 1) - teta*u(i))//compara o lb de cada vb com o xj(t) correspondente
                return false;
        }else{
            if(data->getVectorL(index) > vb(k, 1) + teta*u(i)/*xB*/ || data->getVectorU(index) < vb(k, 1) + teta*u(i))//compara o lb de cada vb com o xj(t) correspondente
                return false;
        }
    }
    //verifica se xj(t) é satisfeita pra o lowerbound e upperbound de j
    if(islower){
        if(data->getVectorL(j) > teta/*xj*/ || data->getVectorU(j) < teta)//compara o lb e o ub de cada vb com o xj(t) correspondente
            return false;
    }else{
        if(data->getVectorL(j) > -1*teta/*xj*/ || data->getVectorU(j) < -1*teta)//compara o lb de cada vb com o xj(t) correspondente
            return false;
    }
    return true;
}
int findTeta(Data * data, int j, bool islower, MatrixXd& vb, VectorXd u, double * teta, bool * unbounded){
    *teta = -1*numeric_limits<double>::infinity();
    int l;
    for(int i = 0; i < vb.rows();i++){
        if(u(i) > 0){
            if(isSatisfied(data, j, vb(i, 1)/u(i)/*possivel teta*/, i, vb, islower, u)){
                *teta = max(*teta, vb(i,1)/u(i)); 
                l = i;
            }else{
                *unbounded = false;
            }
        }
    }
    return l;
}
int chooseJ(VectorXd& reduced_cost){
    for(int i = 0; i < reduced_cost.rows(); i++){
        if(reduced_cost(i) < -EPSILON){
            return i;
        }   
    }
    return 0; //pode deixar assim?  
}

bool isLower(VectorXd& p, Data* data, int j){
    int ya = p.transpose() * data->getMatrixA().col(j);
    if(data->getFO(j) > ya){
        return false;
    }
    return true;
}

void changingVariables(MatrixXd& vb, VectorXd& u, double teta, int l, int j, bool islower){
    for(int i = 0; i < vb.rows(); i++){
        if(islower){
            if(i == l){
                vb(i, 0) = j;
                vb(i, 1) = teta;
            }else{
                vb(i, 1) = vb(i, 1) - teta*u(i);
            }
        }else{
            if(i == l){
                vb(i, 0) = j;
                vb(i, 1) = -teta;
            }else{
                vb(i, 1) = vb(i, 1) + teta*u(i);
            }
        }
    }
}

Solution simplex(Data* data, MatrixXd& B, MatrixXd& variaveis_basicas){
    B = B.inverse();
    int n = data->getMatrixA().rows();
    VectorXd c(n); // c = coeficientes da fo das variaveis basicas
    VectorXd p(n); // p = duais
    VectorXd u(n); // u = B^-1*Aj
    VectorXd reduced_cost(data->getMatrixA().cols());
    
    Solution s;

    while(true){
        findC(c, variaveis_basicas, data, n);
        calculateP(c, B, p);
        calculateReducedC(reduced_cost, p, variaveis_basicas, data);
        
        if(isOptimal(reduced_cost)){
            s.variaveis_basicas = variaveis_basicas;
            s.z = c.transpose() * variaveis_basicas.col(1);
            break;
        }
        else{

            int j = chooseJ(reduced_cost);
            calculateU(B, data, u , j);
            bool islower = isLower(p, data, j);
            bool unbounded = true;
            double teta;
            int l = findTeta(data, j, islower,variaveis_basicas, u, &teta, &unbounded);
            if(unbounded){
                s.z = -1*numeric_limits<double>::infinity();
                break;
            }
            changingVariables(variaveis_basicas, u, teta, l, j, islower);
            loadB(B, u, l);
        }
    }
    return s;
}
int main()
{
    ////passando dados manualmente enquanto nao tenho a leitura de arquivos/////
    
    Data* data = new Data(3, 7);
    int n = data->getMatrixA().rows();
    
    MatrixXd B(n, n);
    B << 3, 1, 0, 1, 1, 0, 4, 3, 1;
    MatrixXd variaveis_basicas(n, 2);
    variaveis_basicas.col(0) << 0, 2, 6;  //linha 1 da matriz serão os indices da matriz -> variavel -1
    variaveis_basicas.col(1) = B.inverse() * data->getRHS();

    
    Solution s = simplex(data, B, variaveis_basicas);

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