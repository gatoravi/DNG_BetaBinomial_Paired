#ifndef BB_PAIRED_H
#define BB_PAIRED_H

#include<iostream>
#include "RInside.h"
#include "denovogear.h"
#include "makeLookup.h"
#include "pedParser.h"
using namespace std;

extern "C" char* bam_mpileup(int argc, char *argv[], int rd_cutoff, const char* mpileup_op_f);


struct ArgsPassed {
  float PP_cutoff;
  int RD_cutoff;
};

struct variant_info{
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


class BBPaired {
  public:
    BBPaired(int argc, char* argv[]);
  private:
    bool skip_mpileup;
    string ip_RD_f; //remove this
    string normal_bam;
    string tumor_bam;
    string tumor_mp_op;
    string normal_mp_op;
    string op_calls_f;
    string ped_f;
    string ref_f;
    string tumorID;
    string normalID;
    double alphaRR;
    double betaRR;
    double alphaRA;
    double betaRA;
    double alphaAA;
    double betaAA;
    int rd_cutoff;
    float pp_cutoff;
    pair_t tumor_pt;
    pair_t normal_pt;
    variant_info VI1;
    Pair* pairs;
    int pair_count;
    int process(int argc, char* argv[]);
    int getSampleIDFromPED();
    int parseArguments(int argc, char* argv[], ArgsPassed& A1);
    int calcBBLik(pair_t& pt);
    int parseFile();
    string parseLine(string line);
    bool get_region_reads();
    int callMpileup(string bam1, string mp_op);
};
#endif
