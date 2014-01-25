#include <iostream>
#include <vector>
#include "RInside.h"
#include "rinside.BBFit.h"

using namespace std;

class AlphaBeta {
 private:
  double RRalpha; 
  double RRbeta;
  double RAalpha; 
  double RAbeta;
  double AAalpha; 
  double AAbeta;
  string mpileup_f;
  RInside* RInst;
 public:
  explicit AlphaBeta(string mpileup_f1, RInside* R2);
  ~AlphaBeta();
  int parseMpileupOP();
  int max3(int a, int b, int c);
  int processLine(string l, int& ref_count, int& alt_count);
  int estimateAlphaBeta(std::vector<int>& k1, std::vector<int>& n1, RInside& Rinst, double& alpha, double& beta);
  std::string get_mp_name();
  int get_Alpha_Beta(std::vector<double>& ab);
  int writeKN(std::vector<int>, std::vector<int>, std::string);
};
