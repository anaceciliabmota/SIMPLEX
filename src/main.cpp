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
  #include <algorithm>
  #include "../Eigen/Dense"
  #include "../Eigen/Sparse"
  #include <iomanip>

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

  struct PrimeiraFase
  {
    bool is_first_phase;
    vector<int> P;
    vector<int> Q;
    VectorXd ub;
    VectorXd lb;
    VectorXd fo;

  } typedef PrimeiraFase;

  void calculateReducedC(VectorXd &reduced_cost, VectorXd &p, MatrixXd &vn, Data *data, PrimeiraFase& aux);

  void changingVariables(MatrixXd &variaveis_basicas, MatrixXd &variaveis_nao_basicas, VectorXd &u, double teta, int l, int j, bool is_negative,  bool restrictive);

  int chooseJ(MatrixXd &variaveis_nao_basicas, VectorXd &reduced_cost, bool *isoptimal, Data * data);

  double defineXj(double uj, double lj);

  void findC(VectorXd &c, MatrixXd &vb, PrimeiraFase& aux);

  int findTeta(Data& data, int j, MatrixXd &vb, VectorXd& u, double *teta, bool *unbounded, bool is_negative, bool *restrictive, PrimeiraFase& aux, int cont);

  bool isCjNegative(double reduced_cost);

  MatrixXd loadB2(MatrixXd &variaveis_basicas, Data *data);

  Basis PhaseOne(Data *data, MatrixXd &variaveis_basicas, MatrixXd &variaveis_nao_basicas, int n, int m, bool *is_feasible, int * iteration);

  void printIteration(Data& data, MatrixXd& variaveis_basicas, MatrixXd& variaveis_nao_basicas, int l, int j, Solution& s, int cont, double teta);

  void printSolution(int n, Solution& s);

  void simplex(Solution * s, Data *data, Basis * b, MatrixXd &variaveis_basicas, MatrixXd &variaveis_nao_basicas, PrimeiraFase * aux, int * iteration);

  int main(int argc, char **argv)
  {
    clock_t start, end;
    
    if (argc < 2)
    {
      cout << "Too few arguments" << endl;
      return 1;
    }

    cout << "Reading problem data from " << argv[1] << "...";
    //data provisorio para teste:
    /*Data data(argv[1]);
    
    VectorXd fo = (*data.getFO());
    fo = fo*-1;
    data.setFO(fo);*/

    mpsReader mps(argv[1]);
  
    Data data(mps.A, mps.b, mps.c, mps.ub, mps.lb);

    int n = data.getMatrixA()->cols(); 
    int m = data.getMatrixA()->rows();

    //variaveis para usar na resolução do simplex  
    MatrixXd variaveis_basicas(m, 2);
    MatrixXd variaveis_nao_basicas(n-m, 2);

    int iteration = 0;

    bool is_feasible = true;
    start = clock();
    Basis b = PhaseOne(&data, variaveis_basicas, variaveis_nao_basicas, n, m, &is_feasible, &iteration);

    Solution s;

    PrimeiraFase aux;
    aux.is_first_phase = false;
    aux.lb = *data.getVectorL();
    aux.ub = *data.getVectorU();
    aux.fo = *data.getFO();

    if (is_feasible)
    {
      simplex(&s, &data, &b,  variaveis_basicas, variaveis_nao_basicas, &aux, &iteration);
      end = clock();

      if (s.z != -1 * numeric_limits<double>::infinity())
      {
        //printSolution(n, s);
        cout << "z: " << s.z << endl;
        cout << "OPTIMAL MPS SOLUTION FOUND" << endl;
      }
      else
      {
        cout << "UNBOUNDED SOLUTION" << endl;
      }
    }
    else
    {
      cout << "UNDEFINED SOLUTION" << endl;
    }
    
    double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
    cout << "TIME USED: " << fixed << time_taken;
    cout << " secs" << endl;
    
    return 0;
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

  void findC(VectorXd &c, MatrixXd &vb, PrimeiraFase& aux)
  {
    // n = quant. de variaveis da base
    for (int i = 0; i < vb.rows(); i++)
    {
      int index = static_cast<int>(vb(i, 0));
      c(i) = aux.fo(index);
    }
  }


  void calculateReducedC(VectorXd &reduced_cost, VectorXd &p, MatrixXd &vn, Data *data, PrimeiraFase& aux)
  { 
    reduced_cost = aux.fo - data->getMatrixA()->transpose() * p;
  }

  int chooseJ(MatrixXd &variaveis_nao_basicas, VectorXd &reduced_cost, bool *isoptimal, Data * data)
  {
    *isoptimal = true;

    int min_index = INT_MAX;

    for (int i = 0; i < variaveis_nao_basicas.rows(); i++)
    {
      int index = static_cast<int>(variaveis_nao_basicas(i, 0));
      if (reduced_cost(index) < -EPSILON && (variaveis_nao_basicas(i, 1) + EPSILON) < (*data->getVectorU())(index))
      {
        *isoptimal = false;
        if (index < min_index){
          min_index = index;
        }
          
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

  int findTeta(Data& data, int j, MatrixXd &vb, VectorXd& u, double *teta, bool *unbounded, bool is_negative, bool *restrictive, PrimeiraFase& aux, int cont)
  {
    VectorXd guarda_ub = aux.ub;
    VectorXd guarda_lb = aux.lb;
    if(aux.is_first_phase){
      for(int i = 0; i < aux.P.size(); i++){
        aux.ub(aux.P[i]) = (*data.getVectorL())(aux.P[i]);

      }

      for(int i = 0; i < aux.Q.size(); i++){
        aux.lb(aux.Q[i]) = (*data.getVectorU())(aux.Q[i]);
      }
    }

    int sinal = is_negative ? -1 : 1;

    vector<double> tetas;
    *teta = numeric_limits<double>::infinity();
    int l = -1;

    for(int i = 0; i < vb.rows(); i++){
      double temp;
      if(abs(u(i)) <= EPSILON){
        temp = numeric_limits<double>::infinity();
      }
      else if (u(i) * -sinal < -EPSILON)
      {
        temp = (aux.ub(static_cast<int>(vb(i, 0))) - vb(i, 1)) / abs(u(i));
      }
      else if(u(i) * -sinal > EPSILON)
      {
        temp = (vb(i, 1) - aux.lb(static_cast<int>(vb(i, 0)))) / abs(u(i));
      }
      tetas.push_back(temp);

      if(temp == *teta){
        if( l != -1 && vb(l, 0) > vb(i, 0)){
          l = i;
        }
      }else if(temp < *teta){
        *teta = temp;
        l = i;
      }
    }
    

    *restrictive = (*teta < aux.ub(j) - aux.lb(j) - EPSILON);

    if(!*restrictive){
      *teta =  aux.ub(j) - aux.lb(j);
    }
    if (*teta == numeric_limits<double>::infinity())
    {
      *unbounded = true;
    }
    else
    {
      *unbounded = false;
    }
    

    aux.lb = guarda_lb;
    aux.ub = guarda_ub;
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

  void changingVariables(MatrixXd &variaveis_basicas, MatrixXd &variaveis_nao_basicas, VectorXd &u, double teta, int l, int j, bool is_negative, bool restrictive)
  {

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

  void simplex(Solution * s,Data *data, Basis * b, MatrixXd &variaveis_basicas, MatrixXd &variaveis_nao_basicas, PrimeiraFase * aux, int * iteration)
  {
    int cont = 0;
    int n = data->getMatrixA()->rows();
    VectorXd c(n); // c = coeficientes da fo das variaveis basicas
    VectorXd p(n); // p = duais
    VectorXd u(n); // u = B^-1*Aj
    VectorXd reduced_cost(n); // reduced_cost = c - A.transpose * p

    while (true)
    { 
      bool isoptimal;

      //os coeficientes das variaveis basicas
      findC(c, variaveis_basicas, *aux);

      //acha o vetor das duais
      b->getP(p, c, cont);

      //calcula o custo reduzido
      calculateReducedC(reduced_cost, p, variaveis_nao_basicas, data, *aux);

      //determina a variavel que entra na base
      int j = chooseJ(variaveis_nao_basicas, reduced_cost, &isoptimal, data);

      //código para calcular o custo a cada iteração
      s->z = c.transpose() * variaveis_basicas.col(1);

      VectorXd non_basic_c(variaveis_nao_basicas.rows());
      findC(non_basic_c, variaveis_nao_basicas, *aux);

      s->z += non_basic_c.transpose() * variaveis_nao_basicas.col(1);


      if (isoptimal)
      {
        //passa as variaveis basicas para a solução
        s->variaveis_basicas = variaveis_basicas;
        s->variaveis_nao_basicas = variaveis_nao_basicas;

        //calcula a fo considerando as basicas
        s->z = c.transpose() * variaveis_basicas.col(1);

        //calcula a fo considerando as nao basicas
        VectorXd non_basic_c(variaveis_nao_basicas.rows());
        findC(non_basic_c, variaveis_nao_basicas, *aux);
        s->z += non_basic_c.transpose() * variaveis_nao_basicas.col(1);

        cout << *iteration << ": obj = " << setprecision(10) << s->z << endl;
        
        break;

      }
      else
      {
        bool unbounded = true, restrictive;
        double teta;

        //acha o vetor direção d
        VectorXd a = data->getMatrixA()->col(j);
        b->getD(u, a, cont, variaveis_basicas, *data);
        
        bool is_negative = isCjNegative(reduced_cost(j));

        //acha a variavel que sai
        int l = findTeta(*data, j, variaveis_basicas, u, &teta, &unbounded, is_negative, &restrictive, *aux, cont);

        
    
        if (unbounded)
        {
          s->z = -1 * numeric_limits<double>::infinity();
          break;
        }
        
        //muda as variaveis
        changingVariables(variaveis_basicas, variaveis_nao_basicas, u, teta, l, j, is_negative, restrictive);
        
        //altera a base
        if (restrictive)
        { 
          int new_non_basic = static_cast<int>(variaveis_basicas(l, 0));
          pair <int, VectorXd> pair(l, u);
          if(b->v.size() >= 13){
            b->loadB(*data, variaveis_basicas);
          }else{
            b->addElement(pair);
          }
        }else{
          b->loadB(*data, variaveis_basicas);
        }
        //etapa adicional da primeira fase

          double infeasibility = 0;

          aux->P.clear();
          aux->Q.clear();
          for(int i = 0; i < variaveis_basicas.rows(); i++){

            int index = static_cast<int>(variaveis_basicas(i, 0));

            if(variaveis_basicas(i, 1) + EPSILON < (*data->getVectorL())(index)){
              infeasibility += (*data->getVectorL())(index) - variaveis_basicas(i, 1);
              aux->P.push_back(index);
              aux->lb(index) = -numeric_limits<double>::infinity();
            }else{
              aux->lb(index) = (*data->getVectorL())(index); 
            }

            if(variaveis_basicas(i, 1) - EPSILON >  (*data->getVectorU())(index)){
              infeasibility += (variaveis_basicas(i, 1) - (*data->getVectorU())(index));
              aux->Q.push_back(index);
              aux->ub(index) = numeric_limits<double>::infinity();
            }else{
              aux->ub(index) = (*data->getVectorU())(index); 
            }
          }
          //printa os custos a cada 100 iterações
          if((*iteration)%1000 == 0){
            cout << *iteration << ": obj = " << setprecision(10) << s->z << " inf = " << infeasibility << endl;
          }
          if(aux->is_first_phase){
            if(infeasibility < EPSILON){
              break; 
            }

            VectorXd new_fo(data->getMatrixA()->cols());
            new_fo.setZero();
            for(int i = 0; i < aux->P.size(); i++){
              new_fo(aux->P[i]) = -1;
            }
            for(int i = 0; i < aux->Q.size(); i++){
              new_fo(aux->Q[i]) = 1;
            }
            aux->fo = new_fo;
          }
          
      }
      (*iteration)++;
    }
    
  }

  double defineXj(double uj, double lj)
  {
    double xj;
    if (uj == numeric_limits<double>::infinity() && lj == -numeric_limits<double>::infinity())
    {
      xj = 0;
    }
    else if (lj == -numeric_limits<double>::infinity())
    {
      xj = uj;
    }
    else
    {
      xj = lj;
    }
    return xj;
  }

  Basis PhaseOne(Data *data, MatrixXd &variaveis_basicas, MatrixXd &variaveis_nao_basicas, int n, int m, bool *is_feasible, int * iteration)
  {
    //variaveis nao basicas sao as n-m primeiras
    PrimeiraFase aux;
    aux.is_first_phase = true;

    VectorXd c(n);
    c.setZero();
    aux.ub = *data->getVectorU();
    aux.lb = *data->getVectorL();

    for(int i = 0; i < n-m; i++){
      variaveis_nao_basicas(i, 0) = i;
      variaveis_nao_basicas(i, 1) = defineXj((*data->getVectorU())(i), (*data->getVectorL())(i));
    }
    
    for(int i = 0; i < m; i++)
      variaveis_basicas(i, 0) = i + (n-m);

    
    SparseMatrix<double> A = *data->getMatrixA();
    //MatrixXd A = A_sparse;
    variaveis_basicas.col(1) =  A.leftCols(n-m) * variaveis_nao_basicas.col(1);


    for(int i = 0; i < variaveis_basicas.rows(); i++){
      int index = static_cast<int>(variaveis_basicas(i, 0));
      if(variaveis_basicas(i, 1) < (*data->getVectorL())(index)){
        aux.lb(index) = -numeric_limits<double>::infinity();
        aux.P.push_back(index);
      }else if(variaveis_basicas(i, 1) > (*data->getVectorU())(index)){
        aux.ub(index) = numeric_limits<double>::infinity();
        aux.Q.push_back(index);
      }
    }
    //definindo a fo
    for(int i = 0; i < aux.P.size(); i++){
      c(aux.P[i]) = -1;
    }
    for(int i = 0; i < aux.Q.size(); i++){
      c(aux.Q[i]) = 1;
    }
    aux.fo = c;

    MatrixXd B = A.rightCols(m);
    Basis b(B);

    Solution s1;
    simplex(&s1, data, &b,  variaveis_basicas, variaveis_nao_basicas, &aux, iteration); 
    
    
    return b;
  }

  void printSolution(int n, Solution& s){
    cout << "Solucao:" << endl;
        for (int i = 0; i < n; i++)
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