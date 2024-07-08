#include "Data.h"
#include "mpsReader.h"
#include <iostream>
#include <algorithm>
#include <limits>
#include <numbers>
#include <vector>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#define EPSILON 10e-5


struct Solution {
  double z;
  MatrixXd variaveis_basicas;
  MatrixXd variaveis_nao_basicas;
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

void findC(VectorXd& c, MatrixXd& vb, Data* data)
{
  // n = quant. de variaveis da base
  for (int i = 0; i < vb.rows(); i++)
  {
    int index = static_cast<int>(vb(i, 0));
    c(i) = data->getFO()(index);

  }
}

void calculateP(VectorXd& c, MatrixXd& B, VectorXd& p)
{
  p = c.transpose() * B; //por padrao, c é coluna, então a transposicao é para transformá-lo em uma linha
}

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
      reduced_cost(i) = data->getFO()(i) - p.transpose() * data->getMatrixA().col(i);
    }
  }
}

int chooseJ(Data * data, MatrixXd& variaveis_nao_basicas, VectorXd& p, bool * isoptimal){

  *isoptimal = true;

  int min_index = INT_MAX;

  for(int i = 0; i < variaveis_nao_basicas.rows(); i++){

    int index = static_cast<int>(variaveis_nao_basicas(i, 0));
    double ya = p.transpose() * data->getMatrixA().col(index);  

    if( ya - data->getFO()(index) > EPSILON && variaveis_nao_basicas(i, 1) < data->getVectorU()(index)){    
      *isoptimal = false;
      if (index < min_index)
        min_index = index;

    }else if(ya - data->getFO()(index) < -EPSILON  && variaveis_nao_basicas(i, 1) > data->getVectorL()(index)){
      *isoptimal = false;
      if (index < min_index)
        min_index = index;
    }
  }
  return min_index;
}

int findTeta(Data * data, int j, MatrixXd& vb, VectorXd u, double * teta, bool * unbounded, bool is_negative, bool * restrictive){
    *teta = numeric_limits<double>::infinity();
    int l;

    if(is_negative){
        for(int i = 0; i < vb.rows();i++){
        double aux;
        if(u(i) == 0){
            aux = numeric_limits<double>::infinity();
            if(aux < *teta){
            *teta = aux;
            l = i;
            }
        }else if(u(i) > 0){
            aux = (vb(i, 1) - data->getVectorL()(static_cast<int>(vb(i, 0))))/u(i);
            if(aux < *teta){
            *teta = aux;
            l = i;
            }
        }else{
            aux = (vb(i, 1) - data->getVectorU()(static_cast<int>(vb(i, 0))))/u(i);
            if(aux < *teta){
            *teta = aux;
            l = i;
            }
        }
        }
    }else {
        for(int i = 0; i < vb.rows(); i++){
        double aux;
        if(u(i) == 0){
            aux = numeric_limits<double>::infinity();
            if(aux < *teta){
            *teta = aux;
            l = i;
            }
        }else if(u(i) > 0){
            aux = (data->getVectorU()(static_cast<int>(vb(i, 0))) - vb(i, 1))/u(i);
            if(aux < *teta){
            *teta = aux;
            l = i;
            }
        }else{
            aux = (data->getVectorL()(static_cast<int>(vb(i, 0))) - vb(i, 1))/u(i);
            if(aux < *teta){
            *teta = aux;
            l = i;
            }
        }
        }
    }

    *restrictive = (*teta < data->getVectorU()(j) - data->getVectorL()(j));
    //retorna true se a basica for a mais restritiva
    *teta = min(*teta, data->getVectorU()(j) - data->getVectorL()(j));
    if(*teta == numeric_limits<double>::infinity()){
        *unbounded = true;
    }else{
        *unbounded = false;
    }
    return l;
}

bool isCjNegative(double reduced_cost){
  if(reduced_cost > EPSILON){
    return false;
  }
  return true;
}

void changingVariables(MatrixXd& variaveis_basicas, MatrixXd& variaveis_nao_basicas, VectorXd& u, double teta, int l, int j, bool is_negative, Data * data, bool non_basic_direction, bool restrictive){

    int new_non_basic = static_cast<int>(variaveis_basicas(l, 0));

    auto linha_j = find(variaveis_nao_basicas.col(0).data(), variaveis_nao_basicas.col(0).data() + variaveis_nao_basicas.rows(), j);
    Index index_j = linha_j - variaveis_nao_basicas.col(0).data();

    for (int i = 0; i < variaveis_basicas.rows(); i++)
    {
        variaveis_basicas(i, 1) += - teta * u(i) * (is_negative ? 1 : -1);
    }

    variaveis_nao_basicas(index_j, 1) += teta * (is_negative ? 1 : -1);
    if(restrictive){
        swap(variaveis_basicas(l, 0), variaveis_nao_basicas(index_j, 0));
        swap(variaveis_basicas(l, 1), variaveis_nao_basicas(index_j, 1));
    }   
}

Solution simplex(Data* data, MatrixXd& B, MatrixXd& variaveis_basicas, MatrixXd& variaveis_nao_basicas){

  int n = data->getMatrixA().rows();
  VectorXd c(n); // c = coeficientes da fo das variaveis basicas
  VectorXd p(n); // p = duais
  VectorXd u(n); // u = B^-1*Aj
  VectorXd reduced_cost(data->getMatrixA().cols());
  Solution s;

  while(true){
    findC(c, variaveis_basicas, data);

    calculateP(c, B, p);

    calculateReducedC(reduced_cost, p, variaveis_basicas, data);

    bool isoptimal;

    int j = chooseJ(data, variaveis_nao_basicas, p, &isoptimal);
    if(isoptimal){
        s.variaveis_basicas = variaveis_basicas;
        s.variaveis_nao_basicas = variaveis_nao_basicas;
        s.z = c.transpose() * variaveis_basicas.col(1);
        VectorXd non_basic_c(variaveis_nao_basicas.rows());
        findC(non_basic_c, variaveis_nao_basicas, data);
        s.z += non_basic_c.transpose() * variaveis_nao_basicas.col(1);
        break;
    }
    else{
        calculateU(B, data, u , j);
        bool is_negative = isCjNegative(reduced_cost(j));

        bool unbounded = true, restrictive;
        double teta;
        int l = findTeta(data, j,variaveis_basicas, u, &teta, &unbounded, is_negative, &restrictive);
        if(unbounded){
            s.z = -1*numeric_limits<double>::infinity();
            break;
        }
        int new_non_basic = static_cast<int>(variaveis_basicas(l, 0));
        changingVariables(variaveis_basicas, variaveis_nao_basicas, u, teta, l, j, is_negative, data, isCjNegative(reduced_cost(new_non_basic)), restrictive);
        if(restrictive){
            loadB(B, u, l);
        }
    }
  }
  return s;
}

double defineXj(double uj, double lj){
    double xj;
    if(uj != numeric_limits<double>::infinity()){
        xj = uj;
    }else if(lj != numeric_limits<double>::infinity()){
        xj = lj;
    }else{
        xj = 0;
    }
    return xj;
}

MatrixXd PhaseOne(Data * data, MatrixXd& variaveis_basicas, MatrixXd& variaveis_nao_basicas, int n, int m, bool * is_feasible){

  //variaveis para inserir na Data do problema auxiliar
  VectorXd rhs = data->getRHS();
  VectorXd u(m + n);
  VectorXd l(m + n);
  VectorXd c(m + n); //fo
  
  //atribuindo as variaveis locais os valores das variaveis xj (nao artificiais) 
  for(int i = 0; i < n; i++){
    variaveis_nao_basicas(i, 0) = i;
    variaveis_nao_basicas(i, 1) = defineXj(data->getVectorU()(i), data->getVectorL()(i));
  }

  u.head(data->getVectorU().size()) = data->getVectorU();
  l.head(data->getVectorL().size()) = data->getVectorL();

  //definindo valor das variaveis artificiais
  for(int i = 0; i < m; i++){
    variaveis_basicas(i, 0) = i + n;
    variaveis_basicas(i, 1) = data->getRHS()(i) - data->getMatrixA().row(i)*variaveis_nao_basicas.col(1);
  }
  //definindo bounds das variaveis artificiais
  for(int i = 0; i < m; i++){
    u(n+i) = ((variaveis_basicas(i, 1) >= 0) ? numeric_limits<double>::infinity() : 0);
    l(n+i) = ((variaveis_basicas(i, 1) >= 0) ? 0 : -1*numeric_limits<double>::infinity());
  }
 
  //definindo nova funcao objetivo
  c.head(n).setZero(); //sera zero para toda variavel nao artificial

  for(int i = 0; i < m; i++){
    c(n+i) = (l(n+i) == 0 ? 1 : -1);
  }

  MatrixXd B = MatrixXd::Identity(m, m); //base do simplex

  MatrixXd A(m, n+m);
  A << data->getMatrixA(), B; //concatenando a matriz A com a matriz das variaveis artificiais

  Data data_auxiliary(A, rhs, c, u, l);

  Solution s = simplex(&data_auxiliary, B, variaveis_basicas, variaveis_nao_basicas);

  cout << s.z << endl;
  if(s.z > EPSILON){
    *is_feasible = false;
  }

  for(int i = 0; i < m; i++){
    u(n+i) = 0;
    l(n+i) = 0;
  }
  //definindo nova fo
  VectorXd fo(n+m);
  fo.head(n) = data->getFO();
  fo.tail(m).setZero();

  data->setMatrixA(A);
  data->setVectorU(u);
  data->setVectorL(l);
  data->setFO(fo);

  return B;
}

int main(int argc, char** argv)
{
    if(argc < 2){
      cout << "Too few arguments" << endl;
      return 1;
    }
    mpsReader mps(argv[1]);

    cout << mps.A << endl;
    cout << mps.ub.transpose() << endl;
    cout << mps.lb.transpose() << endl;
    cout << mps.b << endl << endl;
    for(int i = 0; i < mps.restricoes.size(); i++){
      cout << mps.restricoes[i] << " ";
    }
    cout << endl;

    
    Data data(mps.A, mps.b, mps.c, mps.ub, mps.lb);

    int n = data.getMatrixA().cols();
    int m = data.getMatrixA().rows();

    //variaveis para usar na resolução do simplex
    MatrixXd variaveis_basicas(m, 2);
    MatrixXd variaveis_nao_basicas(n, 2);

    bool is_feasible = true;
    
    MatrixXd B = PhaseOne(&data, variaveis_basicas, variaveis_nao_basicas, n, m, &is_feasible);
    
    if(is_feasible){
      Solution s = simplex(&data, B, variaveis_basicas, variaveis_nao_basicas);

      if(s.z != -1*numeric_limits<double>::infinity()){
        cout << "Solucao:" << endl; 
        for(int i = 0; i < data.getMatrixA().cols(); i++){
            if((s.variaveis_basicas.col(0).array() == i).any()){
                auto linha_i = find(s.variaveis_basicas.col(0).data(), s.variaveis_basicas.col(0).data() + s.variaveis_basicas.rows(), i);
                Index index_i = linha_i - s.variaveis_basicas.col(0).data();
                cout << "x" << s.variaveis_basicas(index_i, 0) + 1 << ": " << s.variaveis_basicas(index_i, 1) << endl;
            }else{
            auto linha_i = find(s.variaveis_nao_basicas.col(0).data(), s.variaveis_nao_basicas.col(0).data() + s.variaveis_nao_basicas.rows(), i);
                Index index_i = linha_i - s.variaveis_nao_basicas.col(0).data();
                cout << "x" << s.variaveis_nao_basicas(index_i, 0) + 1 << ": " << s.variaveis_nao_basicas(index_i, 1) << endl; 
            }
        
        }
        cout << "z: " << s.z << endl;
      }else{
          cout << "Solução ótima é igual a menos infinito" << endl;
      }
    }else{
      cout << "Solucao indefinida" << endl;
    }
    
  return 0;
}