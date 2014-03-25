#This file is a part of Visualization Toolkit project(https://github.com/Kitware/VTK/), we thank the authors for making their code open source.
#Program:   Visualization Toolkit
#Module:    Copyright.txt
#
#Copyright (c) 1993-2008 Ken Martin, Will Schroeder, Bill Lorensen
#All rights reserved.
#
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are met:
#
#* Redistributions of source code must retain the above copyright notice,
#this list of conditions and the following disclaimer.
#
#* Redistributions in binary form must reproduce the above copyright notice,
#this list of conditions and the following disclaimer in the documentation
#and/or other materials provided with the distribution.
#
#* Neither name of Ken Martin, Will Schroeder, or Bill Lorensen nor the names
#of any contributors may be used to endorse or promote products derived
#from this software without specific prior written permission.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ``AS IS''
#AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
#ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
#ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
#CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
#OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#=========================================================================*/#

#execute R command to find the base directory, R has to be installed for this.
find_program(R_COMMAND R DOC "R executable.")

if(R_COMMAND)
  execute_process(WORKING_DIRECTORY .
    COMMAND ${R_COMMAND} RHOME 
    OUTPUT_VARIABLE R_ROOT_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  find_path(R_INCLUDE_DIR R.h
    HINTS ${R_ROOT_DIR}
    PATHS /usr/local/lib /usr/local/lib64 /usr/share
    PATH_SUFFIXES include R/include
    DOC "Path to file R.h")
  
  find_path(RCPP_INCLUDE_DIR Rcpp.h
    HINTS ${R_ROOT_DIR} ${R_ROOT_DIR}/lib/Rcpp/include ${R_ROOT_DIR}/library/Rcpp/include
    PATHS /usr/local/lib /usr/local/lib64 /usr/share
    PATH_SUFFIXES include 
    DOC "Path to file R.h")
  
  find_path(RINSIDE_INCLUDE_DIR RInside.h
    HINTS ${R_ROOT_DIR} ${R_ROOT_DIR}/lib/RInside/include ${R_ROOT_DIR}/library/RInside/include
    PATHS /usr/local/lib /usr/local/lib64 /usr/share
    PATH_SUFFIXES include R/library/include R/include
    DOC "Path to file R.h")

  find_library(R_LIBRARY_BASE R
    HINTS ${R_ROOT_DIR}/lib
    DOC "R library (example libR.a, libR.dylib, etc.).")

  find_library(R_LIBRARY_BLAS NAMES Rblas blas
    HINTS ${R_ROOT_DIR}/lib
    DOC "Rblas library (example libRblas.a, libRblas.dylib, etc.).")

  find_library(R_LIBRARY_LAPACK NAMES Rlapack lapack
    HINTS ${R_ROOT_DIR}/lib
    DOC "Rlapack library (example libRlapack.a, libRlapack.dylib, etc.).")

  find_library(R_LIBRARY_RINSIDE NAMES RInside
    HINTS ${RINSIDE_INCLUDE_DIR}/../lib
    DOC "RInside library (example libRInside.a).")

  find_library(R_LIBRARY_RCPP NAMES Rcpp
    HINTS ${RCPP_INCLUDE_DIR}/../lib
    DOC "Rcpp library (example libRcpp.a).")

else(R_COMMAND)
  message(SEND_ERROR "FindR.cmake requires the following variables to be set: R_COMMAND")
endif(R_COMMAND)

set(R_LIBRARIES ${R_LIBRARY_BASE} ${R_LIBRARY_BLAS} ${R_LIBRARY_LAPACK})
