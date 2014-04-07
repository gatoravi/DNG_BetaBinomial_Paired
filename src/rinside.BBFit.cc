// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8; -*-
//
// Simple example with data in C++ that is passed to R, processed and a result is extracted
//
// Copyright (C) 2009         Dirk Eddelbuettel 
// Copyright (C) 2010 - 2011  Dirk Eddelbuettel and Romain Francois
// Modified by Avinash Ramu, WUSTL
// GPL'ed 

#include "rinside.BBFit.h"


/* Fit a beta binomial to the distribution of read counts. k is a vector of number of alt reads at a locus, n is a vector of total number of reads at a locus. */
int BBFit(std::vector<int> k1, std::vector<int> n1, RInside& R, double& alpha, double& beta) 
{
    //    RInside R(argc, argv);                // create an embedded R instance 
    std::string load = "suppressMessages(source(\"~/src/BBFit.R\"))";
    R.parseEvalQ(load);              // load library, no return value
    R["k1"] = k1;                   
    R["n1"] = n1;                   
    std::string txt = "ab = as.vector(bbFit(k1, n1));"
	  "ab";
    Rcpp::NumericVector v1 = R.parseEval(txt);
    alpha = v1[0];
    beta = v1[1];
    return 0;
}


