#include <getopt.h>
#include <fstream>

#include "pedParser.h"
#include "BBPaired.h"
#include "AlphaBeta.h"
#include "DNGCaller.h"
#include "sam.h"

#define MAX_CHAR_ARRAY 1000
//std::ofstream fout;



BBPaired::BBPaired(int argc, char* argv[])
{ 
  ip_RD_f = "empty";
  op_calls_f = "empty";
  ped_f = "empty";
  ref_f = "empty";
  tumorID = "tumor";
  normalID = "normal";
  alphaRR = -1;
  betaRR = -1;
  alphaRA = -1;
  betaRA = -1;
  alphaAA = -1;
  betaAA = -1;
  rd_cutoff = 10;
  pp_cutoff = 0.001;
  skip_mpileup = false;
  process(argc, argv);
}

int BBPaired::process(int argc, char* argv[])
{
  ArgsPassed A1;
  parseArguments(argc, argv, A1);
  getSampleIDFromPED();
  normal_mp_op = normalID + "_mp_op.txt";
  tumor_mp_op = tumorID + "_mp_op.txt";
  if(!skip_mpileup) {
    callMpileup(normal_bam, normal_mp_op);
    callMpileup(tumor_bam, tumor_mp_op);
  }
  int argc1 = 0;
  char arg1[] = "nothing";
  char* argv1[] = {arg1};

  RInside R1(argc1, argv1);
  AlphaBeta ABnormal(normal_mp_op, &R1);
  AlphaBeta ABtumor(tumor_mp_op, &R1);
  std::vector<double> normal_ab;
  std::vector<double> tumor_ab;
  ABnormal.get_Alpha_Beta(normal_ab);
  ABtumor.get_Alpha_Beta(tumor_ab);
  string normal_mp = ABnormal.get_mp_name();
  string tumor_mp = ABtumor.get_mp_name();
  DNGCaller D1(normal_mp, tumor_mp, normal_ab, tumor_ab, normalID, tumorID, A1);
  return 0;
}

int BBPaired::parseArguments(int argc, char* argv[], ArgsPassed& A1)
{
  static struct option long_options[] = {
                                          {"normal_bam", 1, 0, 0}, 
                                          {"tumor_bam", 1, 0, 1}, 
					  {"op_calls_file", 1, 0, 2},
					  {"ped", 1, 0, 3},
					  {"alphaRR", 1, 0, 4},
					  {"betaRR", 1, 0, 5},
					  {"alphaRA", 1, 0, 6},
					  {"betaRA", 1, 0, 7},
					  {"alphaAA", 1, 0, 8},
					  {"betaAA", 1, 0, 9},
					  {"ref", 1, 0, 10},
					  {"RD_cutoff", 1, 0, 11},
					  {"skip_mpileup", 1, 0, 12},
					  {"PP_cutoff", 1, 0, 13},

  };
  optind = 0;
  // Read in Command Line arguments
  while (1) {
    int option_index = 0;
    int c = getopt_long (argc, argv, "", long_options, &option_index);
    if (c == -1)
      break;
    switch(c) 
      {      
        case 0:
	  normal_bam = optarg;
	  break;

        case 1:
	  tumor_bam = optarg;
	  break;

        case 2:
	  op_calls_f = optarg;
	  break;

        case 3:
	  ped_f = optarg;
	  break;

        case 4:
	  alphaRR = atof(optarg);
	  break;

        case 5:
	  betaRR = atof(optarg);
	  break;

        case 6:
	  alphaRA = atof(optarg);
	  break;

        case 7:
	  betaRA = atof(optarg);
	  break;

        case 8:
	  alphaAA = atof(optarg);
	  break;

        case 9:
	  betaAA = atof(optarg);
	  break;

        case 10:
	  ref_f = optarg;
	  break;

        case 11:
	  rd_cutoff = atoi(optarg);
	  A1.RD_cutoff = rd_cutoff;
	  break;

        case 12:
	  if(strcmp(optarg, "true") == 0) {
	    skip_mpileup = true;
	  }
	  break;

        case 13:
	  pp_cutoff = atof(optarg);
	  A1.PP_cutoff = pp_cutoff;
	  break;
	
        default:
	  break;
      }
  }
  optind = 0;
  return 0;
}



/* calls mpileup on normal BAM - minQ 20*/
int BBPaired::callMpileup(string bam1, string mp_op)
{
  char arg1[] = "mpileup\0";
  char arg2[] = "-f\0";
  char* arg3 = strdup(ref_f.c_str());
  char* arg4 = strdup(bam1.c_str());
  char arg5[] = "-Q\0";
  char arg6[] = "20\0";
  char* argv1[] = {arg1, arg2, arg3, arg4, arg5, arg6};
  int argc1 = sizeof(argv1)/sizeof(argv1[0]);
  bam_mpileup(argc1, argv1, rd_cutoff, mp_op.c_str());
  return 0;
}


int BBPaired::getSampleIDFromPED()
{
  Trio* trios;
  Pair* pairs;
  int trio_count = 0, pair_count = 0;
  parse_ped(ped_f, &trios, &pairs, trio_count, pair_count);
  tumorID = pairs[0].tumorID;
  normalID = pairs[0].normalID;
}



/*
int BBPaired::parseFile()
{
  ifstream fin(ip_RD_f.c_str());
  if(!fin) {
    cerr<<"Unable to open file: "<<ip_RD_f<<".Exiting !"<<endl;
    exit (EXIT_FAILURE);
  }
  string line, sampleID;
  while(fin.good()) {
    getline(fin, line);
    sampleID = parseLine(line);
    if(sampleID == tumorID) // do it for just the first pair for now.
      calcBBLik(tumor_pt);
    else if(sampleID == normalID) // do it for just the first pair for now.
      calcBBLik(normal_pt);
  }
  return 0;
}
*/

/*string BBPaired::parseLine(string line) //line is of the format "chr\tpos\tref\tref_RD\talt\talt_RD"
{
  string sampleID;
  VI1.reset();
  istringstream iss(line);
  iss >> VI1.chr;
  iss >> VI1.pos;
  iss >> VI1.ref;
  iss >> VI1.ref_RD;
  iss >> VI1.alt;
  iss >> VI1.alt_RD;  
  iss >> sampleID;
  return sampleID;
}*/


 /*int BBPaired::calcBBLik(pair_t& pt)
{
  VI1.print();
  return 0;
  }*/

int fetch_func(const bam1_t *bt, void *data)
{
  
  char *name = (char *) bam1_qname(bt);
  int n=0;
  char *qseq = (char *) malloc(bt->core.l_qseq+1);
  char *s = (char *) bam1_seq(bt);
  for(n=0;n<(bt->core.l_qseq);n++) {
    char v = bam1_seqi(s,n);
    qseq[n] = bam_nt16_rev_table[v];
  }
  qseq[n] = 0;
  char *t = (char *) malloc(bt->core.l_qseq + strlen(name) + 10);
  strcpy(t, ">");
  strcat(t, name);
  strcat(t, "\n");
  strcat(t, qseq);
  data = t;
  cout<<(char*) data<<endl;
  free(t);
  free(qseq);
  return 0;
}

int pileup_func(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data)
{
  return 0;
}

bool BBPaired::get_region_reads()
{
  samfile_t *sam_in;  
  bam_plbuf_t *buf;
  sam_in = samopen(normal_bam.c_str(), "rb", 0);
  if (sam_in == 0) {
    cerr<<"Fail to open BAM file "<<normal_bam;
    return EXIT_FAILURE;
  }
  bam_index_t *idx = bam_index_load(normal_bam.c_str()); // load BAM index
  if (idx == 0) {
    cerr<<"BAM indexed file is not available for "<<normal_bam<<endl;
    return EXIT_FAILURE;
  }
  bam_header_t *data = sam_in->header;
  buf = bam_plbuf_init(pileup_func, data); // initialize pileup
  bam_fetch(sam_in->x.bam, idx, 0, 10000, 20000, buf, fetch_func);
  bam_header_destroy(data);
  bam_plbuf_destroy(buf);
  bam_index_destroy(idx);
  //fout.close();
  //samclose(sam_in);
  return true;
}
