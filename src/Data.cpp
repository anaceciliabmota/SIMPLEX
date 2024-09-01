 
#include "Data.h"

#define infinity std::numeric_limits<double>::infinity()

using namespace std;
using namespace Eigen;
Data::Data(MatrixXd& mps_A, VectorXd& mps_b, VectorXd& mps_c, VectorXd& mps_ub, VectorXd& mps_lb)
   {
      A = mps_A.sparseView();
      rhs = mps_b;
      fo = mps_c;
      u = mps_ub;
      l = mps_lb;
  }

Data::Data(std::string fileName) {
    int m, n, col;
    double number;
    //std::string line;
    std::ifstream readFile(fileName);

    if (readFile.is_open()) {
      getDimensions(readFile, &m, &n);
      
      MatrixXd A_dense = LeMatrix(readFile, m, n);
      A = A_dense.sparseView();
      //le vetor c
      fo = LeVetor(readFile, n);

      //le vetor l
      l = LeVetor(readFile, n);

      u = LeVetor(readFile, n);

      rhs = LeVetor(readFile, m);

      indice_vb = LeVetor(readFile, m);

      indice_vn = LeVetor(readFile, n - m);

      variaveis = LeVetor(readFile, n);

      readFile.close();
    } else {
        std::cerr << "Unable to open file " << fileName << std::endl;
    }

}

void Data::getDimensions(ifstream &readFile, int* m, int* n){
      string line;
      getline(readFile, line);
      std::istringstream thisLine(line);
      thisLine >> *m >> *n;
      // Redimensiona as matrizes e vetores

}

MatrixXd Data::LeMatrix(ifstream &readFile, int m, int n){
   // Leitura das linhas da matriz A
      string line;
      double number;
      int col;
      MatrixXd matrizA(m, n);
      for (int row = 0; row < m; ++row) {
         getline(readFile, line);
         std::istringstream thisLine(line);
         cout << thisLine.str() << endl;
         col = 0;
         while (thisLine >> number) {
            matrizA(row, col) = number;
            col++;
         }
      }
      return matrizA;
}

VectorXd Data::LeVetor(ifstream &readFile, int dim){
   VectorXd vetor(dim);
   string number;
   int col = 0;
   string line;
   getline(readFile, line);
   std::istringstream thisLine(line);
   cout << thisLine.str() << endl;
   while(thisLine >> number){
      if(number == "inf"){
         vetor(col) = infinity;

      }
      else if(number == "-inf"){
         vetor(col) = -infinity;
      }
      else{
         vetor(col) = safe_stod(number);
      }
      //cout << number << endl;
      col++;
   }
   return vetor;
}


double Data::safe_stod(const std::string& str) {
    try {
        size_t pos;
        double val = std::stod(str, &pos);
        // Verifica se a string inteira foi convertida
        if (pos != str.length()) {
            throw std::invalid_argument("A string contém caracteres inválidos.");
        }
        return val;
    } catch (const std::invalid_argument& e) {
        std::cerr << "Erro: argumento inválido. A string fornecida não é um número válido: " << str << std::endl;
        throw; // Lança a exceção novamente para ser tratada no nível superior, se necessário
    } catch (const std::out_of_range& e) {
        std::cerr << "Erro: número fora do intervalo permitido: " << str << std::endl;
        throw; // Lança a exceção novamente para ser tratada no nível superior, se necessário
    }
}






