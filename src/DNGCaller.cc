#include "DNGCaller.h"
#include "common.h"
#include "lookup.h"
#include "denovogear.h"
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include <utility>
#include <vector>
#include <map>
#include <math.h>
#include <cstring>
#define MIN_PHRED_LIKE 255
#define MAX_PHRED_LIKE 0

DNGCaller::DNGCaller(string normal_mp, string tumor_mp,
		     std::vector<double> normal_ab1,
		     std::vector<double> tumor_ab1, 
		     string aNormalID, string aTumorID, 
		     ArgsPassed& A1)
{
  //cout<<"constructing ";
  //cout<<"\n1t size "<<tumor_ab1.size()<<"\n";
  tumorID = aTumorID;
  normalID = aNormalID;
  tumor_mp_f = tumor_mp;
  normal_mp_f = normal_mp;
  tumor_ab = tumor_ab1;
  normal_ab = normal_ab1;
  n_site_pass = 0;
  RD_cutoff = A1.RD_cutoff;
  PP_cutoff = A1.PP_cutoff;
  callMakePairedLookup(tgtPair1, lookupPair1);
  createChrMap();
  iterateMpileups();
}

DNGCaller::~DNGCaller()
{
}

int DNGCaller::createChrMap()
{
  chr_map["chr1"] = chr1;
  chr_map["chr2"] = chr2;
  chr_map["chr3"] = chr3;
  chr_map["chr4"] = chr4;
  chr_map["chr5"] = chr5;
  chr_map["chr6"] = chr6;
  chr_map["chr7"] = chr7;
  chr_map["chr8"] = chr8;
  chr_map["chr9"] = chr9;
  chr_map["chr10"] = chr10;
  chr_map["chr11"] = chr11;
  chr_map["chr12"] = chr12;
  chr_map["chr13"] = chr13;
  chr_map["chr14"] = chr14;
  chr_map["chr15"] = chr15;
  chr_map["chr16"] = chr16;
  chr_map["chr17"] = chr17;
  chr_map["chr18"] = chr18;
  chr_map["chr19"] = chr19;
  chr_map["chr20"] = chr20;
  chr_map["chr21"] = chr21;
  chr_map["chr22"] = chr22;
  chr_map["chrX"] = chrX;
  chr_map["chrY"] = chrY;
  chr_map["chrM"] = chrM;
  return 0;
}

//assumes ordering "chr1 chr2 chr3.. chr10.. chr22 .. chrX .. chrY"
int DNGCaller::iterateMpileups()
{
    //cout<<"\nt size "<<tumor_ab.size()<<"\n";
  ifstream tin1(tumor_mp_f.c_str());
  if(!tin1) {
    cerr<<"Unable to open "<<tumor_mp_f<<" .Exiting !";
    exit(EXIT_FAILURE);
  }
  ifstream nin1(normal_mp_f.c_str());
  if(!nin1) {
    cerr<<"Unable to open "<<normal_mp_f<<" .Exiting !";
    exit(EXIT_FAILURE);
  }
  string l_t, l_n;
  string t_chr, n_chr;
  long t_pos, n_pos;
  getline(tin1, l_t);
  getline(nin1, l_n);

  while(tin1.good() && nin1.good()) {
    istringstream iss_t(l_t);
    istringstream iss_n(l_n);
    iss_t >> t_chr;
    iss_t >> t_pos;
    iss_n >> n_chr;
    iss_n >> n_pos;
    //cout<<t_chr<<"\t"<<chr_map[t_chr]<<"\t"<<n_chr<<"\t"<<chr_map[n_chr]<<endl;
    if(t_chr == n_chr && t_pos == n_pos) {
      //cout<<l_n<<"\t"<<l_t<<endl;
      calcLik(l_t, l_n);
      getline(tin1, l_t);
      getline(nin1, l_n);
    }
    else if(t_chr == n_chr) {
	if(t_pos > n_pos) {
	  //cout<<t_pos<<" greater than "<<n_pos<<endl;
	  getline(nin1, l_n);
	  continue;
	}
	else {
	  //cout<<n_pos<<" greater than "<<t_pos<<endl;
	  getline(tin1, l_t);
	  continue;
	} 
    }

    else if(chr_map[t_chr] > chr_map[n_chr]) {
      //cout<<t_chr<<" greater than "<<n_chr<<endl;
      getline(nin1, l_n);
      continue;
    }
    else if(chr_map[n_chr] > chr_map[t_chr]){
      //cout<<n_chr<<" greater than "<<t_chr<<endl;
      getline(tin1, l_t);
      continue;
    } 
  }
  return 0;
}

int DNGCaller::parseLine(string l, t_variant_info& vi1)
{
  #ifdef DEBUG
  cout<<l<<endl;
  #endif
  istringstream iss(l);
  string chr, ref, pileup, field;
  long pos;
  int RD = 0, a_count = 0, c_count = 0, g_count = 0, t_count = 0, ref_count = 0;
  typedef std::pair<std::string, int> basecount;
  std::vector<basecount> basecounts;

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
    ref = "A";
    ref_count = a_count;
  }
  if(ref == "c" || ref == "C") {
    ref = "C";
    ref_count = c_count;
  }
  if(ref == "g" || ref == "G") {
    ref = "G";
    ref_count = g_count;
  }
  if(ref == "t" || ref == "T") {
    ref = "T";
    ref_count = t_count;
  }

  basecounts.push_back(std::make_pair("A", a_count));
  basecounts.push_back(std::make_pair("C", c_count));
  basecounts.push_back(std::make_pair("G", g_count));
  basecounts.push_back(std::make_pair("T", t_count));
  std::sort (basecounts.begin(), basecounts.end(), compareBaseCounts);

  /* std::cout<<endl<<"Line "<<l;
  std::cout<<endl<<basecounts[0].first<<"\t"<<basecounts[0].second;
  std::cout<<endl<<basecounts[1].first<<"\t"<<basecounts[1].second;
  std::cout<<endl<<basecounts[2].first<<"\t"<<basecounts[2].second;
  std::cout<<endl<<basecounts[3].first<<"\t"<<basecounts[3].second;
  */

  vi1.ref_RD = 0;
  vi1.alt_RD = 0;
  vi1.chr = chr;
  vi1.pos = pos;
  vi1.ref = ref;
  vi1.ref_RD = ref_count;
  if(basecounts[0].first == ref) {
    vi1.alt = basecounts[1].first;
    vi1.alt_RD = basecounts[1].second;
  }
  else {
    vi1.alt = basecounts[0].first;
    vi1.alt_RD = basecounts[0].second;
  }
  if(vi1.alt_RD == 0)
    vi1.alt = "N";
  return 0;
}


double DNGCaller::lbetaFn(double alpha, double beta)
{
  return lgamma(alpha) + lgamma(beta) - lgamma(alpha + beta);
}

double DNGCaller::lik_bb(double alpha, double beta, int n, int k)
{
  while(k>500 || n>500) {
    k = k/2;
    n = n/2;
  }
  double t1  = k + alpha;
  double t2  = n - k + beta;
  double lnum = lbetaFn(t1, t2);		
  double ldenom = lbetaFn(alpha, beta);
  double log_lik  = lnum - ldenom;
  return log_lik;
}

int DNGCaller::convertVI2paired(t_variant_info vi1, pair_t& pt1, std::vector<double>& ab, string sampleID)
{
  std::map<string, int> GT_lik;
  std::map<std::string, int>::iterator GT_it;
  GT_lik["AA"] =  MIN_PHRED_LIKE;
  GT_lik["AC"] =  MIN_PHRED_LIKE;
  GT_lik["AG"] =  MIN_PHRED_LIKE;
  GT_lik["AT"] =  MIN_PHRED_LIKE;
  GT_lik["CC"] =  MIN_PHRED_LIKE;
  GT_lik["CG"] =  MIN_PHRED_LIKE;
  GT_lik["CT"] =  MIN_PHRED_LIKE;
  GT_lik["GG"] =  MIN_PHRED_LIKE;
  GT_lik["GT"] =  MIN_PHRED_LIKE;
  GT_lik["TT"] =  MIN_PHRED_LIKE;

  string RR = vi1.ref + vi1.ref;
  string RA;
  if(vi1.ref <= vi1.alt) {
  	RA = vi1.ref + vi1.alt;
  } else {
  	RA = vi1.alt + vi1.ref;
  }	
  string AA = vi1.alt + vi1.alt;

  //std::cout<<"1"<<endl;
  //std::cout<<ab.size()<<endl;

  double log_lik_bb_RR   =  -10 * lik_bb(ab[0], ab[1], vi1.ref_RD + vi1.alt_RD, vi1.alt_RD);
  double log_lik_bb_RA   =  -10 * lik_bb(ab[2], ab[3], vi1.ref_RD + vi1.alt_RD, vi1.alt_RD); 
  double log_lik_bb_AA   =  -10	* lik_bb(ab[4], ab[5], vi1.ref_RD + vi1.alt_RD, vi1.alt_RD); 

  #ifdef DEBUG
  std::cout<<"RR "<<RR;
  std::cout<<"RA "<<RA;
  std::cout<<"AA "<<AA;
  cout<<"\nbefore norm liks "<<log_lik_bb_RR<<"\t"<<log_lik_bb_RA<<"\t"<<log_lik_bb_AA;  
  #endif
  double norm_bb_liks[3];
  normalize_liks(log_lik_bb_RR, log_lik_bb_RA, log_lik_bb_AA, norm_bb_liks);
  #ifdef DEBUG
  cout<<"\nafter norm liks "<<norm_bb_liks[0]<<"\t"<<norm_bb_liks[1]<<"\t"<<norm_bb_liks[2];
  #endif

  if(vi1.alt != "N") {
    GT_lik[RA] = norm_bb_liks[1];
    GT_lik[AA] = norm_bb_liks[2];
  }
  else {
    for(GT_it = GT_lik.begin(); GT_it != GT_lik.end(); GT_it++) {
      string GT = GT_it->first;
      std::size_t found = GT.find(vi1.ref);
      if (found != std::string::npos) {
	GT_lik[GT] = norm_bb_liks[1];
	//cout<<endl<<GT<<" contains "<<vi1.ref;
      }
    }
  }
  GT_lik[RR] = norm_bb_liks[0];
  
  strcpy(pt1.chr, vi1.chr.c_str());
  pt1.pos = vi1.pos;
  pt1.ref_base = vi1.ref[0];
  strcpy(pt1.alt, vi1.alt.c_str());
  pt1.depth = vi1.ref_RD + vi1.alt_RD;
  pt1.rms_mapQ = -1;
  //cout<<endl<<"sample ID "<<sampleID;
  pt1.id = sampleID;

  int index = 0;
  for(GT_it = GT_lik.begin(); GT_it != GT_lik.end(); GT_it++) {
  #ifdef DEBUG
   cout<<GT_it->first<<"\t"<<GT_it->second<<endl;
  #endif
    pt1.lk[index++] = GT_it->second;
  }
  #ifdef DEBUG
  cout<<endl;
  #endif
  return 0;
}

int DNGCaller::normalize_liks(double log_lik_RR, double log_lik_RA, double log_lik_AA, double* normalized_liks)
{
  double min_log_lik = min(log_lik_RR, min(log_lik_RA, log_lik_AA));
  double norm_log_lik_RR = (int)(log_lik_RR - min_log_lik + 0.499);
  double norm_log_lik_RA = (int)(log_lik_RA - min_log_lik + 0.499);
  double norm_log_lik_AA = (int)(log_lik_AA - min_log_lik + 0.499);
  			
  if(norm_log_lik_RR > 255)
    norm_log_lik_RR = MIN_PHRED_LIKE;

  if(norm_log_lik_RA > 255)
    norm_log_lik_RA = MIN_PHRED_LIKE;

  if(norm_log_lik_AA > 255)
    norm_log_lik_AA = MIN_PHRED_LIKE;

  if(norm_log_lik_RR < 0)
    norm_log_lik_RR = MAX_PHRED_LIKE;

  if(norm_log_lik_RA < 0)
    norm_log_lik_RA = MAX_PHRED_LIKE;

  if(norm_log_lik_AA < 0)
    norm_log_lik_AA = MAX_PHRED_LIKE;

  normalized_liks[0] = norm_log_lik_RR;
  normalized_liks[1] = norm_log_lik_RA;
  normalized_liks[2] = norm_log_lik_AA;

  return 0;
}


int DNGCaller::calcLik(string l_t, string l_n)
{
  //cout<<"\ntumor  "<<l_t;
  //cout<<"\nnormal "<<l_n;
  t_variant_info tvi, nvi;
  pair_t pt_t, pt_n;
  #ifdef DEBUG
  cout<<endl<<"tumor";
  #endif
  parseLine(l_t, tvi);
  convertVI2paired(tvi, pt_t, tumor_ab, tumorID);
  #ifdef DEBUG
  cout<<endl<<"normal";
  #endif
  parseLine(l_n, nvi);
  convertVI2paired(nvi, pt_n, normal_ab, normalID);
  /*
  cout<<endl<<"1 "<<tvi.alt_RD+tvi.ref_RD<<"\t"<<nvi.alt_RD+nvi.alt_RD<<endl;
  cout<<endl<<"2 "<<pt_t.depth<<"\t"<<pt_n.depth<<endl;
  */
  int flag1 = 0; // flag to indicate site specific info
  ofstream fo_vcf1;
  string op_vcf_f1 = "empty";
  pair_like(pt_t, pt_n, tgtPair1, lookupPair1, flag1, op_vcf_f1, fo_vcf1, PP_cutoff, RD_cutoff, n_site_pass, tvi.alt_RD, nvi.alt_RD);

  return 0;
}
