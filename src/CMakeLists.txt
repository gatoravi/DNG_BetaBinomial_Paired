INCLUDE(CheckIncludeFile)
INCLUDE(CheckFunctionExists)
INCLUDE(CheckLibraryExists)

#add_definitions(-DDEBUG)

set(SAMTOOLS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/samtools-0.1.18")

#RInside_cmake
set(NUM_TRUNC_CHARS 2)

set (R_PATH "R")
set (RSCRIPT_PATH "Rscript")

if (CMAKE_HOST_WIN32)
  execute_process(COMMAND ${RSCRIPT_PATH} -e "cat(.Platform$r_arch)"
    OUTPUT_VARIABLE R_ARCH)

  execute_process(COMMAND ${R_PATH} --arch ${R_ARCH} RHOME
    OUTPUT_VARIABLE R_HOME)

  string(REPLACE "\\" "/" R_HOME ${R_HOME})

  set (R_PATH ${R_HOME}/bin/R)
endif()

execute_process(COMMAND ${R_PATH} CMD config --cppflags
  OUTPUT_VARIABLE RCPPFLAGS)

if (CMAKE_HOST_WIN32)
  if (${RCPPFLAGS} MATCHES "[-][I]([^ ;])+")
    set (RCPPFLAGS ${CMAKE_MATCH_0})
  endif()
endif()

string(SUBSTRING ${RCPPFLAGS} ${NUM_TRUNC_CHARS} -1 RCPPFLAGS)
include_directories(${RCPPFLAGS})

execute_process(COMMAND ${R_PATH} CMD config --ldflags
  OUTPUT_VARIABLE RLDFLAGS)
string(LENGTH ${RLDFLAGS} RLDFLAGS_LEN)

if (${RLDFLAGS} MATCHES "[-][L]([^ ;])+")
  string(SUBSTRING ${CMAKE_MATCH_0} ${NUM_TRUNC_CHARS} -1 RLDFLAGS_L)
  string(STRIP ${RLDFLAGS_L} RLDFLAGS_L )
  link_directories(${RLDFLAGS_L} )
endif()

string(STRIP ${RLDFLAGS} RLDFLAGS_l )
if (${RLDFLAGS} MATCHES "[-][l]([^;])+")
  string(SUBSTRING ${CMAKE_MATCH_0} ${NUM_TRUNC_CHARS} -1 RLDFLAGS_l)
  string(STRIP ${RLDFLAGS_l} RLDFLAGS_l )
endif()
#if (${RLDFLAGS} MATCHES "([^;]+[.][a])([^;])+")
  #set(RBASELIB ${CMAKE_MATCH_1})
#endif()

execute_process(COMMAND ${RSCRIPT_PATH} -e "Rcpp:::CxxFlags()"
  OUTPUT_VARIABLE RCPPINCL)
string(SUBSTRING ${RCPPINCL} ${NUM_TRUNC_CHARS} -1 RCPPINCL)
include_directories(${RCPPINCL})

execute_process(COMMAND ${RSCRIPT_PATH} -e "RcppClassic:::LdFlags()"
  OUTPUT_VARIABLE RCPPLIBS)

execute_process(COMMAND ${RSCRIPT_PATH} -e "RInside:::CxxFlags()"
  OUTPUT_VARIABLE RINSIDEINCL)
string(SUBSTRING ${RINSIDEINCL} ${NUM_TRUNC_CHARS} -1 RINSIDEINCL)
include_directories(${RINSIDEINCL})

execute_process(COMMAND ${RSCRIPT_PATH} -e "RInside:::LdFlags()"
  OUTPUT_VARIABLE RINSIDELIBS)

if (CMAKE_HOST_WIN32)
  string(LENGTH "libRcpp.a" lenRcppName)
  string(LENGTH ${RCPPLIBS} lenRcppFQName)

  math(EXPR RLibPathLen ${lenRcppFQName}-${lenRcppName}-1)
  string(SUBSTRING ${RCPPLIBS} 0 ${RLibPathLen} RCPPLIBS_L)
  link_directories(${RCPPLIBS_L})

  math(EXPR RLibPathLen ${RLibPathLen}+1)
  string(SUBSTRING ${RCPPLIBS} ${RLibPathLen} -1 RCPPLIBS_l)

  #Remove the quotes
  string(SUBSTRING ${RINSIDELIBS} 1 -1 RINSIDELIBS)
  string(LENGTH ${RINSIDELIBS} lenRInsideFQNameLen)
  math(EXPR lenRInsideFQNameLen ${lenRInsideFQNameLen}-1)
  string(SUBSTRING ${RINSIDELIBS} 0 ${lenRInsideFQNameLen} RINSIDELIBS)

  string(LENGTH "libRInside.a" lenRInsideName)
  string(LENGTH ${RINSIDELIBS} lenRInsideFQName)

  math(EXPR RLibPathLen ${lenRInsideFQName}-${lenRInsideName}-1)
  string(SUBSTRING ${RINSIDELIBS} 0 ${RLibPathLen} RINSIDELIBS_L)

  math(EXPR RLibPathLen ${RLibPathLen}+1)
  string(SUBSTRING ${RINSIDELIBS} ${RLibPathLen} -1 RINSIDELIBS_l)

  link_directories(${RINSIDELIBS_L})
else()
  if (${RCPPLIBS} MATCHES "[-][L]([^ ;])+")
    string(SUBSTRING ${CMAKE_MATCH_0} ${NUM_TRUNC_CHARS} -1 RCPPLIBS_L)
    link_directories(${RCPPLIBS_L} )
  endif()

  if (${RCPPLIBS} MATCHES "[-][l][R]([^;])+")
    string(SUBSTRING ${CMAKE_MATCH_0} ${NUM_TRUNC_CHARS} -1 RCPPLIBS_l)
  endif()

  if (${RCPPLIBS} MATCHES "([^;])+[.][a]")
    set(RCPPLIBS_l ${RCPPLIBS})
  endif()

  if (${RINSIDELIBS} MATCHES "[-][L]([^ ;])+")
    string(SUBSTRING ${CMAKE_MATCH_0} ${NUM_TRUNC_CHARS} -1 RINSIDELIBS_L)
    link_directories(${RINSIDELIBS_L})
  endif()

  if (${RINSIDELIBS} MATCHES "[-][l][R]([^;])+")
    string(SUBSTRING ${CMAKE_MATCH_0} ${NUM_TRUNC_CHARS} -1 RINSIDELIBS_l)
  endif()

  if (${RINSIDELIBS} MATCHES "([^;])+[.][a]")
    set(RINSIDELIBS_l ${RINSIDELIBS})
  endif()


endif()

execute_process(COMMAND ${R_PATH} CMD config CXXFLAGS
  OUTPUT_VARIABLE RCXXFLAGS)

execute_process(COMMAND ${R_PATH} CMD config BLAS_LIBS
  OUTPUT_VARIABLE RBLAS)

execute_process(COMMAND ${R_PATH} CMD config LAPACK_LIBS
  OUTPUT_VARIABLE RLAPACK)

set(CMAKE_CXX_FLAGS "-W -Wall -pedantic -Wextra ${CMAKE_CXX_FLAGS}")

if (CMAKE_BUILD_TYPE STREQUAL "DEBUG" OR
    CMAKE_BUILD_TYPE STREQUAL "RelWithDebugInfo" )
  add_definitions("-DDEBUG")
elseif ( CMAKE_BUILD_TYPE STREQUAL "RELEASE" )
  add_definitions("-O3")
endif()
ADD_DEFINITIONS(-DHAVE_CONFIG_H
  "\"-DPACKAGE_STRING=\\\"${PROJ_NAME} ${PROJ_VERSION}\\\"\""
  "\"-DPACKAGE_BUGREPORT=\\\"${PROJ_EMAIL}\\\"\""
  )

IF(CMAKE_COMPILER_IS_GNUCC)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
ENDIF(CMAKE_COMPILER_IS_GNUCC)
IF(CMAKE_COMPILER_IS_GNUCXX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

IF(MSVC)
  ADD_DEFINITIONS(-D_CRT_SECURE_NO_DEPRECATE -D_CRT_NONSTDC_NO_DEPRECATE)
  ADD_DEFINITIONS(-D_SCL_SECURE_NO_DEPRECATE)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4100 /wd4101 /wd4189 /wd4512 /wd4275 /wd4251 /wd4706")
  IF(USE_STATIC_LIBS)
    STRING(REGEX REPLACE "/MD" "/MT"   CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}")
    STRING(REGEX REPLACE "/MD" "/MT"   CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELEASE}")
    STRING(REGEX REPLACE "/MD" "/MT"   CMAKE_C_FLAGS_MINSIZEREL "${CMAKE_C_FLAGS_RELEASE}")
    STRING(REGEX REPLACE "/MDd" "/MTd" CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}")

    STRING(REGEX REPLACE "/MD" "/MT"   CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
    STRING(REGEX REPLACE "/MD" "/MT"   CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELEASE}")
    STRING(REGEX REPLACE "/MD" "/MT"   CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_C_FLAGS_RELEASE}")
    STRING(REGEX REPLACE "/MDd" "/MTd" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
  ENDIF(USE_STATIC_LIBS)

  # Use the highest warning level for visual studio.
  SET(CMAKE_CXX_WARNING_LEVEL 4)
  IF(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
    STRING(REGEX REPLACE "/W[0-4]" "/W4"
      CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  ELSE(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
  ENDIF(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
ENDIF(MSVC)

FIND_PACKAGE(ZLIB)
IF(ZLIB_FOUND)
  INCLUDE_DIRECTORIES(${ZLIB_INCLUDE_DIR})
ENDIF(ZLIB_FOUND)

#find_package(R REQUIRED)
#IF(R_INCLUDE_DIR)
#  include_directories(${R_INCLUDE_DIR})
#  message(STATUS "AVI R INCLUDE DIR: ${R_INCLUDE_DIR}
#  ${RINSIDE_INCLUDE_DIR} ${RCPP_INCLUDE_DIR} #${R_ROOT_DIR}
#  ${R_LIBRARY_RCPP} ${R_LIBRARY_RINSIDE}")
#ENDIF(R_INCLUDE_DIR)

#find_package(Samtools REQUIRED)
message(STATUS "Samtools DIR: ${SAMTOOLS_DIR}")

ADD_SUBDIRECTORY(newmat)
ADD_SUBDIRECTORY(utils)
INCLUDE_DIRECTORIES(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/utils" "${CMAKE_CURRENT_SOURCE_DIR}/newmat")
INCLUDE_DIRECTORIES(BEFORE ${SAMTOOLS_DIR})
INCLUDE_DIRECTORIES(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/include")

set(executable "denovogear-bb")
ADD_EXECUTABLE(denovogear-bb ${PROJ_MACOSX_BUNDLE}
  denovogear.cc
  snpLike.cc indelLike.cc pairLike.cc
  makeLookup.cc
  phaser.cc
  pedParser.cc
  bcf2Qcall.cc bcf2Paired.cc
  BBPaired.cc
  AlphaBeta.cc
  rinside.BBFit.cc
  DNGCaller.cc
  )

IF(USE_STATIC_LIBS)
  SET_TARGET_PROPERTIES(denovogear-bb PROPERTIES LINK_SEARCH_END_STATIC On)
ENDIF(USE_STATIC_LIBS)

#target_link_libraries(denovogear-bb ${RBASELIB})
TARGET_LINK_LIBRARIES(denovogear-bb newmat )
#TARGET_LINK_LIBRARIES(denovogear-bb newmat dngutils)
TARGET_LINK_LIBRARIES(denovogear-bb "${SAMTOOLS_DIR}/bcftools/libbcf.a")
TARGET_LINK_LIBRARIES(denovogear-bb "${SAMTOOLS_DIR}/libbam.a")
TARGET_LINK_LIBRARIES(denovogear-bb ${ZLIB_LIBRARIES})
TARGET_LINK_LIBRARIES(denovogear-bb ${CMAKE_REQUIRED_LIBRARIES})
target_link_libraries(denovogear-bb ${RLDFLAGS_l})
target_link_libraries(denovogear-bb ${BLAS_LIBS})
target_link_libraries(denovogear-bb ${LAPACK_LIBS})
target_link_libraries(denovogear-bb ${RINSIDELIBS_l})
target_link_libraries(denovogear-bb ${RCPPLIBS_l})
#TARGET_LINK_LIBRARIES(denovogear-bb ${R_LIBRARY_BLAS})
#TARGET_LINK_LIBRARIES(denovogear-bb ${R_LIBRARY_LAPACK})
#TARGET_LINK_LIBRARIES(denovogear-bb ${R_LIBRARY_BASE})
#TARGET_LINK_LIBRARIES(denovogear-bb ${R_LIBRARY_RCPP})
#TARGET_LINK_LIBRARIES(denovogear-bb ${R_LIBRARY_RINSIDE})

INSTALL(TARGETS denovogear-bb RUNTIME DESTINATION bin ${PROJ_INSTALL_EXTRA_ARGS})
