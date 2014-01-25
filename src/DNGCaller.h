#ifndef DNGCALLER_H
#define DNGCALLER_H

#include <iostream>
#include <string>
#include <vector>
#include "parser.h"
#include "lookup.h"
#include "BBPaired.h"
using namespace std;

//ordering of chromosomes to be set here
enum chromosomes { chr1,  chr2,  chr3, chr4, chr5, chr6, chr7, chr8, chr9,
		   chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19,
		   chr20, chr21, chr22,
		   chrM, chrX, chrY };

struct t_variant_info{
  string chr;
  long pos;
  string ref;
  string alt;
  int ref_RD; // no of reads supporting ref
  int alt_RD; // no of reads supporting alt 
  void reset()
  {
    chr = "NA";
    pos = -1;
    ref = "NA";
    alt = "NA";
    ref_RD = -1;
    alt_RD = -1;
  }
  void print()
  {
    cout<<"\nV1\t"<<chr<<"\t"<<pos<<"\t"<<ref<<"\t"<<ref_RD;
    cout<<"\t"<<alt<<"\t"<<alt_RD<<endl;
  }

};


class DNGCaller
{
  public: 
    explicit DNGCaller(string tumor_mp, string normal_mp, 
	    std::vector<double> tumor_ab, 
	    std::vector<double> normal_ab,
	    string tumorID, 
	    string normalID, 
	    ArgsPassed& A1);
    ~DNGCaller();
  private:
    string tumor_mp_f;
    string normal_mp_f;
    int n_site_pass;
    int RD_cutoff;
    float PP_cutoff;
    lookup_pair_t lookupPair1;
    vector<vector<string > > tgtPair1;
    std::vector<double> normal_ab;
    std::vector<double> tumor_ab;
    std::map<string, int> chr_map;
    string normalID;
    string tumorID;
    int createChrMap();
    int iterateMpileups();
    int calcLik(std::string l_t, std::string l_n);
    int parseLine(string l, t_variant_info& vi1);
    int convertVI2paired(t_variant_info vi1, pair_t& pt1, std::vector<double>& ab, 
			 string sampleID);
    double lbetaFn(double alpha, double beta);
    double lik_bb(double alpha, double beta, int n, int k);
    int normalize_liks(double log_lik_RR, double log_lik_RA, double log_lik_AA, double* normalized_liks);
};
#endif
