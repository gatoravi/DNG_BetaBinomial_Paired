#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include "AlphaBeta.h"

#define RR_AF_cutoff 0.14
#define AA_AF_cutoff 0.86

//For estimating alpha, beta in the RR case, use same number of sites as the AA case.

int AlphaBeta::estimateAlphaBeta(std::vector<int>& k1, std::vector<int>& n1, RInside& Rinst, double& alpha, double& beta)
{
  BBFit(k1, n1, Rinst, alpha, beta);
  return 0;
}

int AlphaBeta::max3(int a, int b, int c)
{
  if(a>b){
    if(a>c)
      return a;
    else
      return c;
  }
  else if(b>c)
    return b;
  return c;
}

int AlphaBeta::processLine(string l, int& ref_count, int& alt_count)
{
  istringstream iss(l);
  string chr, pos, ref, pileup, field;
  int RD, a_count, c_count, g_count, t_count;
    
  iss >> chr;
  iss >> pos;
  iss >> field;
  RD = atoi(field.c_str());
  iss >> ref;
  iss >> pileup;

  iss >> field;
  a_count = atoi(field.c_str());
  iss >> field;
  c_count = atoi(field.c_str());
  iss >> field;
  g_count = atoi(field.c_str());
  iss >> field;
  t_count = atoi(field.c_str());
  
  if(ref == "a" || ref == "A") {
    ref_count = a_count;
    alt_count = max3(c_count, g_count, t_count);
  }
  else if(ref == "c" || ref == "C") {
    ref_count = c_count;
    alt_count = max3(a_count, g_count, t_count);
  }
  else if(ref == "g" || ref == "G") {
    ref_count = g_count;
    alt_count = max3(a_count, c_count, t_count);
  }
  else if(ref == "t" || ref == "T") {
    ref_count = t_count;
    alt_count = max3(a_count, c_count, g_count);
  }
  return 0;
}

int AlphaBeta::writeKN(std::vector<int> n, std::vector<int> k, string op_f)
{
  if(n.size() == k.size()) {			
    ofstream fout(op_f.c_str());
    if(!fout.good()) {
      cerr<<"\nUnable to open "<<op_f<<". Exiting.!";
      exit(EXIT_FAILURE);
    }
    for(int j = 0; j<n.size(); j++) {
      fout<<k[j]<<"\t"<<n[j]<<"\n";
    }
    fout.close();
  }
}

int AlphaBeta::parseMpileupOP()
{
  ifstream fin(mpileup_f.c_str());
  if(!fin.good()) {
    cerr<<endl<<"Unable to open "<<mpileup_f<<". Exiting !";
    exit(EXIT_FAILURE);
  }
    
    
  int ref_c, alt_c;
  float AF;
  //typedef std::pair<int, int> n_k_pair;
  //n_k_pair nk1;
  //std::vector<n_k_pair> RR_AFs, RA_AFs, AA_AFs;
  std::vector<int> RR_n, RA_n, AA_n;
  std::vector<int> RR_k, RA_k, AA_k;
  while(fin.good()){
    string l; /* chr pos RD ref pileup a c g t */
    getline(fin, l);
    processLine(l, ref_c, alt_c);

    if(ref_c + alt_c != 0) { 
      AF = (float)alt_c/(float)(ref_c + alt_c);
    }
    else 
      AF = 0;
    //cout<<"\nl"<<l<<"\trc\t"<<ref_c<<"\tac\t"<<alt_c<<"\tAF"<<AF;
    //nk1 = std::make_pair(ref_c + alt_c, alt_c);

    if(AF <= RR_AF_cutoff) {
      //RR_AFs.push_back(nk1);
      RR_n.push_back(ref_c + alt_c);
      RR_k.push_back(alt_c);
    }
    else if(AF >= AA_AF_cutoff) {
      //AA_AFs.push_back(nk1);
      //cout<<"\n\nAA\t"<<"l"<<l<<"\trc\t"<<ref_c;
      //cout<<"\tac\t"<<alt_c<<"\tAF"<<AF<<"\n";
      AA_n.push_back(ref_c + alt_c);
      AA_k.push_back(alt_c);
    }
    else { 
      //RA_AFs.push_back(nk1);
      //cout<<"\n\nRA\t"<<"l"<<l<<"\trc\t"<<ref_c;
      //cout<<"\tac\t"<<alt_c<<"\tAF"<<AF<<"\n";
      RA_n.push_back(ref_c + alt_c);
      RA_k.push_back(alt_c);
    }
  }
  fin.close();
  //#ifdef DEBUG
  cerr<<endl<<mpileup_f;
  cerr<<endl<<"RR n size "<<RR_n.size()<<" RR k size "<<RR_k.size();
  cerr<<endl<<"RA n size "<<RA_n.size()<<" RA k size "<<RA_k.size();
  cerr<<endl<<"AA n size "<<AA_n.size()<<" AA k size "<<AA_k.size();
  //#endif
  RR_n.resize(AA_n.size()); // change size of RR to that of AA
  RR_k.resize(AA_n.size());
  //#ifdef DEBUG
  cerr<<endl<<"RR(resized) n size "<<RR_n.size()<<" RR k size "<<RR_k.size();  
  cerr<<endl;
  //#endif
  estimateAlphaBeta(RA_k, RA_n, *RInst, RAalpha, RAbeta);
  estimateAlphaBeta(AA_k, AA_n, *RInst, AAalpha, AAbeta);
  estimateAlphaBeta(RR_k, RR_n, *RInst, RRalpha, RRbeta);
#ifdef DEBUG
  writeKN(RA_n, RA_k, "RA_k_n.txt");
  writeKN(AA_n, AA_k, "AA_k_n.txt");
  writeKN(RR_n, RR_k, "RR_k_n.txt");
#endif
  //#ifdef DEBUG
  cerr<<endl<<"RR alpha beta "<<RRalpha<<"\t"<<RRbeta;
  cerr<<endl<<"RA alpha beta "<<RAalpha<<"\t"<<RAbeta;
  cerr<<endl<<"AA alpha beta "<<AAalpha<<"\t"<<AAbeta;
  //#endif
  return 0;
}

string AlphaBeta::get_mp_name()
{
  return mpileup_f;
}

int AlphaBeta::get_Alpha_Beta(std::vector<double>& ab)
{
  ab.push_back(RRalpha);
  ab.push_back(RRbeta);
  ab.push_back(RAalpha);
  ab.push_back(RAbeta);
  ab.push_back(AAalpha);
  ab.push_back(AAbeta);
  return 0;
}

AlphaBeta::AlphaBeta(string mpileup_f1, RInside* R2)
{
  mpileup_f = mpileup_f1;
  RInst = R2;
  parseMpileupOP();
}

AlphaBeta::~AlphaBeta()
{

}
