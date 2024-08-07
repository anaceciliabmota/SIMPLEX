#include "Data.h"
#include "mpsReader.h"
#include "Basis.h"
#include <iostream>
#include <algorithm>
#include <limits>
#include <numbers>
#include <vector>
#include <time.h> 
#include <cmath>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#define EPSILON 1e-5
#define EPSILON2 1e-8
#define EPSILON3 1e-6 

struct Solution
{
  double z;
  MatrixXd variaveis_basicas;
  MatrixXd variaveis_nao_basicas;
} typedef Solution;

void calculateP(VectorXd &c, MatrixXd &B, VectorXd &p);

void calculateReducedC(VectorXd &reduced_cost, VectorXd &p, MatrixXd &vn, Data *data);

void calculateU(MatrixXd &B, Data *data, VectorXd &u, int j);

void changingVariables(MatrixXd &variaveis_basicas, MatrixXd &variaveis_nao_basicas, VectorXd &u, double teta, int l, int j, bool is_negative, Data *data, bool non_basic_direction, bool restrictive);

int chooseJ(Data *data, MatrixXd &variaveis_nao_basicas, VectorXd &reduced_cost, bool *isoptimal);

double defineXj(double uj, double lj);

void findC(VectorXd &c, MatrixXd &vb, Data *data);

int findTeta(Data *data, int j, MatrixXd &vb, VectorXd u, double *teta, bool *unbounded, bool is_negative, bool *restrictive, MatrixXd& B);

bool isCjNegative(double reduced_cost);

void loadB(MatrixXd &B, VectorXd &u, int l);

MatrixXd loadB2(MatrixXd &variaveis_basicas, Data *data);

void PhaseOne(Data *data, MatrixXd &variaveis_basicas, MatrixXd &variaveis_nao_basicas, int n, int m, bool *is_feasible, Basis& b);

void printData(Data& data, MatrixXd& variaveis_basicas, MatrixXd& variaveis_nao_basicas);

void printIteration(Data& data, MatrixXd& variaveis_basicas, MatrixXd& variaveis_nao_basicas, int l, int j, Solution& s, int cont, double teta);

void printMatrixSequential(MatrixXd& A) ;

void printSolution(Data& data, Solution& s);

void printVectorSequential(const Eigen::VectorXd& v);

Solution simplex(Data *data, Basis& b, MatrixXd &variaveis_basicas, MatrixXd &variaveis_nao_basicas, bool first_phase);

int main(int argc, char **argv)
{
  clock_t start, end;
  start = clock();
  if (argc < 2)
  {
    cout << "Too few arguments" << endl;
    return 1;
  }
  //data provisorio para teste:
  /*Data data(argv[1]);
  
  VectorXd fo = (*data.getFO());
  fo = fo*-1;
  data.setFO(fo);*/

  mpsReader mps(argv[1]);
 
  Data data(mps.A, mps.b, mps.c, mps.ub, mps.lb);

  int n = data.getMatrixA()->cols();
  int m = data.getMatrixA()->rows();

  // variaveis para usar na resolução do simplex
  MatrixXd variaveis_basicas(m, 2);
  MatrixXd variaveis_nao_basicas(n, 2);
  

  //colocando a solucao basica passada no arquivo para as variaveis
  /*for(int i = 0; i < data.indice_vb.rows(); i++){
    variaveis_basicas(i, 0) = data.indice_vb(i);
    int index = static_cast<int>(data.indice_vb(i));
    variaveis_basicas(i, 1) = data.variaveis(index);
  }
  for(int i = 0; i < data.indice_vn.rows(); i++){
    variaveis_nao_basicas(i, 0) = data.indice_vn(i);
    int index = static_cast<int>(data.indice_vn(i));
    variaveis_nao_basicas(i, 1) = data.variaveis(index);
  }*/

  //printData(data, variaveis_basicas, variaveis_nao_basicas);
  
  //cout << variaveis_basicas << endl << endl;
  //cout << variaveis_nao_basicas << endl;
  

  //MatrixXd B = loadB2(variaveis_basicas, &data);
  //cout << B << endl;
  //B = B.inverse();
  //cout << B << endl;
  MatrixXd B = MatrixXd :: Identity(m, m);
  Basis b(B);
  bool is_feasible = true;
  
  PhaseOne(&data, variaveis_basicas, variaveis_nao_basicas, n, m, &is_feasible, b);

  //printData(data, variaveis_basicas, variaveis_nao_basicas);
  //exit(0);    
  cout << "###############" << endl;
  if (is_feasible)
  {
    //cout << "feasible" << endl;
    Solution s = simplex(&data, b,  variaveis_basicas, variaveis_nao_basicas, false);

    if (s.z != -1 * numeric_limits<double>::infinity())
    {
      //printSolution(data, s);
    }
    else
    {
      cout << "Solução ótima é igual a menos infinito" << endl;
    }
  }
  else
  {
    cout << "Solucao indefinida" << endl;
  }
  end = clock();
  double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
  cout << "tempo gasto: " << fixed << time_taken;
  cout << " secs" << endl;

  //b.~Basis();
  return 0;
}

void loadB(MatrixXd &B, VectorXd &u, int l)
{
  for (int i = 0; i < B.rows(); i++)
  {
    if (i != l)
    {
      double coeficient = -1 * u(i) / u(l);
      B.row(i) = B.row(i) + coeficient * B.row(l);
      u(i) = u(i) + coeficient * u(l);
    }
    else
    {
      B.row(l) = B.row(l) / u(l);
      u(l) = u(l) / u(l);
    }
  }
}

void loadBi(MatrixXd &B, VectorXd &u, int l){
  for (int i = 0; i < B.rows(); i++)
  {
    if (i != l)
    {
      double coeficient = -1 * u(i) / u(l);
      B.row(i) = B.row(i) + coeficient * B.row(l);
      //u(i) = u(i) + coeficient * u(l);
    }
    else
    {
      B.row(l) = B.row(l) / u(l);
      //u(l) = u(l) / u(l);
    }
  }
}

MatrixXd loadB2(MatrixXd &variaveis_basicas, Data *data)
{
  MatrixXd B(data->getMatrixA()->rows(), data->getMatrixA()->rows());
  for (int i = 0; i < variaveis_basicas.rows(); i++)
  {
    B.col(i) = data->getMatrixA()->col(variaveis_basicas(i, 0));
  }
  return B;
}

void findC(VectorXd &c, MatrixXd &vb, Data *data)
{
  // n = quant. de variaveis da base
  for (int i = 0; i < vb.rows(); i++)
  {
    int index = static_cast<int>(vb(i, 0));
    c(i) = (*data->getFO())(index);
  }
}

void calculateP(VectorXd &c, MatrixXd &B, VectorXd &p)
{
  p = B.transpose().colPivHouseholderQr().solve(c);; // por padrao, c é coluna, então a transposicao é para transformá-lo em uma linha
}

void calculateU(MatrixXd &B, Data *data, VectorXd &u, int j)
{
  u = B.colPivHouseholderQr().solve(data->getMatrixA()->col(j));
  //! any u > 0, nao tem solução otima
}

void calculateReducedC(VectorXd &reduced_cost, VectorXd &p, MatrixXd &vn, Data *data)
{
  reduced_cost.setZero();
  for (int i = 0; i < vn.rows(); i++)
  {
    int index = static_cast<int>(vn(i, 0));
    reduced_cost(index) = (*data->getFO())(index) - p.transpose() * data->getMatrixA()->col(index);
  }
}

bool isOptimal(MatrixXd& vn, VectorXd& reduced_cost, Data& data){
  bool optimal = true;
  for(int i = 0; i < vn.rows(); i++){
    int index = static_cast<int>(vn(i, 0));
    if(vn(i, 1) == (*data.getVectorL())(index)){
      if(reduced_cost(index) < -EPSILON){
        optimal = false;
        break;
      }
    }else if(vn(i, 1) == (*data.getVectorU())(index)){
      if(reduced_cost(index) > EPSILON){
        optimal = false;
        break;
      }
    }else if((*data.getVectorU())(index) == numeric_limits<double>::infinity() && (*data.getVectorL())(index) == -numeric_limits<double>::infinity()){
      if(reduced_cost(index) > EPSILON || reduced_cost(index) < -EPSILON){
        optimal = false;
        break;
      }
    }
  }
  return optimal;
}

int chooseJ(Data *data, MatrixXd &variaveis_nao_basicas, VectorXd &reduced_cost, bool *isoptimal)
{
  *isoptimal = true;

  int min_index = INT_MAX;

  for (int i = 0; i < variaveis_nao_basicas.rows(); i++)
  {
    int index = static_cast<int>(variaveis_nao_basicas(i, 0));
    //cout << reduced_cost(index) << " " << variaveis_nao_basicas(i, 1) << " " <<   (*data->getVectorU())(index) << endl;
    //cout << reduced_cost(index) << " " << variaveis_nao_basicas(i, 1) << " " <<   (*data->getVectorL())(index) << endl;

    if (reduced_cost(index) < -EPSILON && (variaveis_nao_basicas(i, 1) + EPSILON) < (*data->getVectorU())(index))
    {
      *isoptimal = false;
      if (index < min_index)
        min_index = index;
    }
    else if (reduced_cost(index) > EPSILON && (variaveis_nao_basicas(i, 1) - EPSILON) > (*data->getVectorL())(index))
    {
      *isoptimal = false;
      if (index < min_index)
        min_index = index;

    }
  
  }

  return min_index;
}

int findTeta(Data *data, int j, MatrixXd &vb, VectorXd u, double *teta, bool *unbounded, bool is_negative, bool *restrictive)
{
  *teta = numeric_limits<double>::infinity();
  int l;

  if (is_negative)
  {
    for (int i = 0; i < vb.rows(); i++)
    {

      double aux;
      if (u(i) < -EPSILON)
      {
        aux = (vb(i, 1) - (*data->getVectorU())(static_cast<int>(vb(i, 0)))) / u(i);
        if(vb(i, 0) == 204 || vb(i, 0) == 0){
          //cout <<"aux " <<  vb(i, 0) << ": " << aux << endl;
        }
        if(aux < *teta){
           //double a = B(i, i)/u(i);
          //if(a < -EPSILON2 || a > EPSILON2){
              *teta = aux;
              l = i;
            
          //}
        }else if(abs(aux - *teta) <= EPSILON){
          if(vb(l, 0) > vb(i, 0)){
            l = i;
          }
        }
      }
      else if (u(i) > EPSILON)
      {
        aux = (vb(i, 1) - (*data->getVectorL())(static_cast<int>(vb(i, 0)))) / u(i);
        if(vb(i, 0) == 204 || vb(i, 0) == 0){
          //cout <<"aux " <<  vb(i, 0) << ": " << aux << endl;
        }
        if(aux < *teta){
           //double a = B(i, i)/u(i);
          //if(a < -EPSILON2 || a > EPSILON2){
            *teta = aux;
            l = i;
          //}
        }else if(abs(aux - *teta) <= EPSILON){
          if(vb(l, 0) > vb(i, 0)){
            l = i;
          }
        }
      }
      else
      {
        aux = numeric_limits<double>::infinity();
        if(vb(i, 0) == 204 || vb(i, 0) == 0){
          //cout <<"aux " <<  vb(i, 0) << ": " << aux << endl;
        }
       if(aux < *teta){
           //double a = B(i, i)/u(i);
          //if(a < -EPSILON2 || a > EPSILON2){
            *teta = aux;
            l = i;
          //}
        }else if(abs(aux - *teta) <= EPSILON){
          if(vb(l, 0) > vb(i, 0)){
            l = i;
          }
        }
      }
    }
  }
  else
  {
    for (int i = 0; i < vb.rows(); i++)
    {
      double aux;
      if (u(i) < -EPSILON)
      {
        aux = ((*data->getVectorL())(static_cast<int>(vb(i, 0))) - vb(i, 1)) / u(i);
        if(vb(i, 0) == 204 || vb(i, 0) == 0){
          //cout <<"aux " <<  vb(i, 0) << ": " << aux << endl;
        }
          if(aux < *teta){
          //double a = B(i, i)/u(i);
          //if(a < -EPSILON2 || a > EPSILON2){
            *teta = aux;
            l = i;
          //}
        }else if(abs(aux - *teta) <= EPSILON){
          if(vb(l, 0) > vb(i, 0)){
            l = i;
          }
        }
      }
      else if (u(i) > EPSILON)
      {
        aux = ((*data->getVectorU())(static_cast<int>(vb(i, 0))) - vb(i, 1)) / u(i);
        if(vb(i, 0) == 204 || vb(i, 0) == 0){
          //cout <<"aux " <<  vb(i, 0) << ": " << aux << endl;
        }
        if(aux < *teta){
           //double a = B(i, i)/u(i);
          //if(a < -EPSILON2 || a > EPSILON2){
            *teta = aux;
            l = i;
          //}
        }else if(abs(aux - *teta) <= EPSILON){
          if(vb(l, 0) > vb(i, 0)){
            l = i;
          }
        }
      }
      else
      {
        aux = numeric_limits<double>::infinity();
        if(vb(i, 0) == 204 || vb(i, 0) == 0){
          //cout <<"aux " <<  vb(i, 0) << ": " << aux << endl;
        }
        if(aux < *teta){
           //double a = B(i, i)/u(i);
          //if(a < -EPSILON2 || a > EPSILON2){
            *teta = aux;
            l = i;
          //}
        }else if(abs(aux - *teta) <= EPSILON){
          if(vb(l, 0) > vb(i, 0)){
            l = i;
          }
        }
        
      }
    }
  }

  *restrictive = (*teta < (*data->getVectorU())(j) - (*data->getVectorL())(j));
  // retorna true se a basica for a mais restritiva
  *teta = min(*teta, (*data->getVectorU())(j) - (*data->getVectorL())(j));
  if (*teta == numeric_limits<double>::infinity())
  {
    *unbounded = true;
  }
  else
  {
    *unbounded = false;
  }
  return l;
}

bool isCjNegative(double reduced_cost)
{
  if (reduced_cost > EPSILON)
  {
    return false;
  }
  return true;
}

void changingVariables(MatrixXd &variaveis_basicas, MatrixXd &variaveis_nao_basicas, VectorXd &u, double teta, int l, int j, bool is_negative, Data *data, bool non_basic_direction, bool restrictive)
{

  int new_non_basic = static_cast<int>(variaveis_basicas(l, 0));

  auto linha_j = find(variaveis_nao_basicas.col(0).data(), variaveis_nao_basicas.col(0).data() + variaveis_nao_basicas.rows(), j);
  Index index_j = linha_j - variaveis_nao_basicas.col(0).data();

  for (int i = 0; i < variaveis_basicas.rows(); i++)
  {
    variaveis_basicas(i, 1) += -teta * u(i) * (is_negative ? 1 : -1);
  }

  variaveis_nao_basicas(index_j, 1) += teta * (is_negative ? 1 : -1);

  if (restrictive)
  {
    swap(variaveis_basicas(l, 0), variaveis_nao_basicas(index_j, 0));
    swap(variaveis_basicas(l, 1), variaveis_nao_basicas(index_j, 1));
  }
}

void printIteration(Data& data, MatrixXd& variaveis_basicas, MatrixXd& variaveis_nao_basicas, int l, int j, Solution& s, int cont, double teta){
  cout << "Iteração " << cont << endl << endl;
  cout << "variavel a entrar na base nessa iteração: " << j << endl;
  int variavel_l = static_cast<int>(variaveis_basicas(l, 0));
  cout << "variavel a sair da base nessa iteração e seu indice nas bsicas: " << variavel_l << " na linha " << l << endl;
  cout << "valor de teta: " << teta << endl;
  cout << "novas variaveis basicas" << endl;
  cout << variaveis_basicas << endl;
  cout << "novas variaveis nao basicas" << endl;
  cout << variaveis_nao_basicas << endl << endl;
}

Solution simplex(Data *data, Basis& b, MatrixXd &variaveis_basicas, MatrixXd &variaveis_nao_basicas, bool first_phase)
{
  int cont = 0;
  /*B = loadB2(variaveis_basicas, data);
  B = B.inverse();*/
  //cout << "TAMANO VETOR::: " << b.v.size() << endl;  
  int n = data->getMatrixA()->rows();
  //cout <<"n: "<< n << endl;
  VectorXd c(n); // c = coeficientes da fo das variaveis basicas
  VectorXd p(n); // p = duais
  VectorXd u(n); // u = B^-1*Aj
  VectorXd reduced_cost(data->getMatrixA()->cols());
  Solution s;

  while (true)
  {
    //cout << "ITERACAO " << cont << endl;
    //cout << endl << "b.v.size: " << b.v.size() << endl;
    findC(c, variaveis_basicas, data);
    //cout << "c: " << c.transpose() << endl;
    //calculateP(c, B, p);
    b.getP(p, c);
    //cout << "p: " << p.transpose() << endl;

    calculateReducedC(reduced_cost, p, variaveis_nao_basicas, data);

    
    //cout << reduced_cost.transpose() << endl;

    bool isoptimal;

    int j = chooseJ(data, variaveis_nao_basicas, reduced_cost, &isoptimal);
    //cout << "j: " << j << endl;

    /*s.z = c.transpose() * variaveis_basicas.col(1);

    VectorXd non_basic_c(variaveis_nao_basicas.rows());
    findC(non_basic_c, variaveis_nao_basicas, data);

    s.z += non_basic_c.transpose() * variaveis_nao_basicas.col(1);*/

    
    /*if(first_phase && s.z < EPSILON){
      break;
    }*/
    if (isoptimal)
    {
      s.variaveis_basicas = variaveis_basicas;
      s.variaveis_nao_basicas = variaveis_nao_basicas;

      s.z = c.transpose() * variaveis_basicas.col(1);

      VectorXd non_basic_c(variaveis_nao_basicas.rows());
      findC(non_basic_c, variaveis_nao_basicas, data);

      s.z += non_basic_c.transpose() * variaveis_nao_basicas.col(1);
      break;
    }
    else
    {
      //calculateU(B, data, u, j);
      VectorXd a = data->getMatrixA()->col(j);
      //cout << a.transpose() << endl;
      b.getD(u, a);
      //cout << "U: " << u.transpose() << endl;
      
     
      //cout << "0";
      bool is_negative = isCjNegative(reduced_cost(j));

      bool unbounded = true, restrictive;
      double teta;
      //cout << "1" ;
      int l = findTeta(data, j, variaveis_basicas, u, &teta, &unbounded, is_negative, &restrictive);
      //cout << "2";
      pair <int, VectorXd> pair(l, u);
     
     
      //cout << b.v[cont].first <<  " " << b.v[cont].second.transpose() << endl;
      //cout << "l: " << l << endl;
      if (unbounded)
      {
        s.z = -1 * numeric_limits<double>::infinity();
        break;
      }
      int new_non_basic = static_cast<int>(variaveis_basicas(l, 0));
      //cout << first_phase << " " << cont << " " << s.z << " " << j << " " << new_non_basic << endl;
      changingVariables(variaveis_basicas, variaveis_nao_basicas, u, teta, l, j, is_negative, data, isCjNegative(reduced_cost(new_non_basic)), restrictive);
      
      if (restrictive)
      { 
        if(b.v.size() == 20){
          b.loadB(*data, variaveis_basicas);
        }else{
          b.addElement(pair);
        }
        
        //cout << pair.first << " " << pair.second.transpose() << endl;
        //B.col(l) = data->getMatrixA()->col(j);
        //loadB(B, u, l);
        //VectorXd aux = B * variaveis_basicas.col(1);
        /*
        if(/*!(aux).isApprox(*data->getRHS(), EPSILON3)*//* cont > 20){
          B = loadB2(variaveis_basicas, data);
          B = B.inverse();
          cont = 0;
        }
        B = loadB2(variaveis_basicas, data);
        B = B.inverse();  */
      }
       //printIteration(*data, variaveis_basicas, variaveis_nao_basicas, l, j, s, cont, teta);
    }
  }
  cout << s.z << endl;
  //cout << "NUMERO DE ITERACOES: " << cont << endl;
  return s;
}

double defineXj(double uj, double lj)
{
  double xj;
  if (uj != numeric_limits<double>::infinity())
  {
    xj = uj;
  }
  else if (lj != -numeric_limits<double>::infinity())
  {
    xj = lj;
  }
  else
  {
    xj = 0;
  }
  return xj;
}

void PhaseOne(Data *data, MatrixXd &variaveis_basicas, MatrixXd &variaveis_nao_basicas, int n, int m, bool *is_feasible, Basis& b)
{

  // variaveis para inserir na Data do problema auxiliar
  VectorXd rhs = *data->getRHS();
  VectorXd u(m + n);
  VectorXd l(m + n);
  VectorXd c(m + n); // fo

  // atribuindo as variaveis locais os valores das variaveis xj (nao artificiais)
  for (int i = 0; i < n; i++)
  {
    variaveis_nao_basicas(i, 0) = i;
    variaveis_nao_basicas(i, 1) = defineXj((*data->getVectorU())(i), (*data->getVectorL())(i));
  }

  u.head(data->getVectorU()->size()) = *data->getVectorU();
  l.head(data->getVectorL()->size()) = *data->getVectorL();

  // definindo valor das variaveis artificiais
  for (int i = 0; i < m; i++)
  {
    variaveis_basicas(i, 0) = i + n;
    variaveis_basicas(i, 1) = (*data->getRHS())(i) - data->getMatrixA()->row(i) * variaveis_nao_basicas.col(1);
  }
  // definindo bounds das variaveis artificiais
  for (int i = 0; i < m; i++)
  {
    u(n + i) = ((variaveis_basicas(i, 1) < -EPSILON) ? 0 : numeric_limits<double>::infinity());
    l(n + i) = ((variaveis_basicas(i, 1) < -EPSILON) ? -numeric_limits<double>::infinity() : 0);
  }

  // definindo nova funcao objetivo
  c.head(n).setZero(); // sera zero para toda variavel nao artificial

  for (int i = 0; i < m; i++)
  {
    c(n + i) = (l(n + i) < -EPSILON ? -1 : 1);
  }

  MatrixXd B = MatrixXd::Identity(m, m); // base do simplex

  MatrixXd A(m, n + m);
  A << *data->getMatrixA(), B; // concatenando a matriz A com a matriz das variaveis artificiais

  Data data_auxiliary(A, rhs, c, u, l);
  /*
  cout << data_auxiliary.getMatrixA() << endl;
  cout << data_auxiliary.getVectorL().transpose() << endl;
  cout << data_auxiliary.getVectorU().transpose() << endl;*/
  /*
  cout << "SOlucao inicial da primeira fase" << endl;
  cout << B << endl;
  cout << "vb: " << variaveis_basicas << endl;
  cout << "vn: " << variaveis_nao_basicas << endl;*/

  //printData(data_auxiliary, variaveis_basicas, variaveis_nao_basicas);

  Solution s = simplex(&data_auxiliary, b, variaveis_basicas, variaveis_nao_basicas, true);
  //cout << s.z << endl;
  if (s.z > EPSILON)
  {
    *is_feasible = false;
  }

  for (int i = 0; i < m; i++)
  {
    u(n + i) = 0;
    l(n + i) = 0;
  }
  // definindo nova fo
  VectorXd fo(n + m);
  fo.head(n) = *data->getFO();
  fo.tail(m).setZero();

  data->setMatrixA(A);
  data->setVectorU(u);
  data->setVectorL(l);
  data->setFO(fo);

}

void printMatrixSequential(MatrixXd& A) {
    ostringstream oss;
    oss << "";
    for (int i = 0; i < A.rows(); ++i) {
        oss << " ";
        for (int j = 0; j < A.cols(); ++j) {
            oss << A(i, j);
            if (i != A.rows() - 1 || j != A.cols() - 1) {
                oss << " ";
            }
        }
        oss << endl;
    }
    std::cout << oss.str() << std::endl;
}

void printVectorSequential(const Eigen::VectorXd& v) {
    std::ostringstream oss;
    oss << "";
    for (int i = 0; i < v.size(); ++i) {
        oss << v(i);
        if (i != v.size() - 1) {
            oss << " ";
        }
    }
    std::cout << oss.str() << std::endl;
}

void printData(Data& data, MatrixXd& variaveis_basicas, MatrixXd& variaveis_nao_basicas){

  cout << "m n" << endl;
  cout << data.getMatrixA()->rows() << " " << data.getMatrixA()->cols()<< endl;
  cout << "A" << endl;
  /*for(int i = 0; i < data.getMatrixA().rows(); i++){
    cout << data.getMatrixA().row(i) << endl;
  }*/
  MatrixXd *A = data.getMatrixA();
  printMatrixSequential(*A);
  cout << "b" << endl;
  VectorXd *b = data.getRHS();
  printVectorSequential(*b);
  cout << "c " << endl;
  VectorXd *fo = data.getFO();
  printVectorSequential(*fo);
   cout << "l"<< endl;
  VectorXd *l = data.getVectorL();
  printVectorSequential(*l);
  cout << " u" << endl;
  VectorXd *u = data.getVectorU();
  printVectorSequential(*u);
  cout << "x" << endl;
  vector<double> variaveis(data.getMatrixA()->cols());
  for (int i = 0; i < variaveis_basicas.rows(); i++)
  {
    variaveis[static_cast<int>(variaveis_basicas(i, 0))] = variaveis_basicas(i, 1);
  }
  for (int i = 0; i < variaveis_nao_basicas.rows(); i++)
  {
    variaveis[static_cast<int>(variaveis_nao_basicas(i, 0))] = variaveis_nao_basicas(i, 1);
  }
  cout << "  ";
  for(int i = 0; i < variaveis.size(); i++){
    cout << variaveis[i] << " ";
  }
  cout << endl << "indice das basicas " << endl;
  VectorXd vb = variaveis_basicas.col(0);
  printVectorSequential(vb);
  cout << "indice das nao basicas " << endl;
  VectorXd vn = variaveis_nao_basicas.col(0);
  printVectorSequential(vn);
}

void printSolution(Data& data, Solution& s){
  cout << "Solucao:" << endl;
      for (int i = 0; i < data.getMatrixA()->cols(); i++)
      {
        if ((s.variaveis_basicas.col(0).array() == i).any())
        {
          auto linha_i = find(s.variaveis_basicas.col(0).data(), s.variaveis_basicas.col(0).data() + s.variaveis_basicas.rows(), i);
          Index index_i = linha_i - s.variaveis_basicas.col(0).data();
          cout << "x" << s.variaveis_basicas(index_i, 0) + 1 << ": " << s.variaveis_basicas(index_i, 1) << endl;
        }
        else
        {
          auto linha_i = find(s.variaveis_nao_basicas.col(0).data(), s.variaveis_nao_basicas.col(0).data() + s.variaveis_nao_basicas.rows(), i);
          Index index_i = linha_i - s.variaveis_nao_basicas.col(0).data();
          cout << "x" << s.variaveis_nao_basicas(index_i, 0) + 1 << ": " << s.variaveis_nao_basicas(index_i, 1) << endl;
        }
      }
      cout << "z: " << s.z << endl;
}
