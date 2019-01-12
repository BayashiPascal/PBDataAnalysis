// ============ PBDATAANALYSIS.H ================

#ifndef PBDATAANALYSIS_H
#define PBDATAANALYSIS_H

// ================= Include =================

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <execinfo.h>
#include <errno.h>
#include <string.h>
#include "pberr.h"
#include "pbmath.h"
#include "gset.h"

// ================= Define ==================

// ================= Data structure ===================

typedef enum KMeansClustersSeed {
  KMeansClustersSeed_Random, KMeansClustersSeed_Forgy,
  KMeansClustersSeed_PlusPlus
} KMeansClustersSeed;
#define KMeansClustersSeed_Default KMeansClustersSeed_PlusPlus

typedef struct KMeansClusters {
  GSetVecFloat _centers;
  KMeansClustersSeed _seed;
} KMeansClusters;

// ================ Functions declaration ====================

// Create a KMeansClusters for a K-means clustering initialized 
// using the 'seed' technique
KMeansClusters KMeansClustersCreateStatic(const KMeansClustersSeed seed);

// Free the memory used by a PBDAKMeansClusters
void KMeansClustersFreeStatic(KMeansClusters* const that);

// Create the set of 'K' clusters clustering the 'input' data according 
// to the K-means algorithm
// 'K' must be inferior or equal to the number of input data
// srandom() must have been called before using this function
void KMeansClustersSearch(KMeansClusters* const that,
  const GSetVecFloat* const input, const int K);

// Get the set of clusters' center for the KMeansClusters 'that'
#if BUILDMODE != 0 
inline 
#endif
const GSetVecFloat* KMeansClustersCenters(
  const KMeansClusters* const that);

// Get the 'iCluster'-th cluster's center for the KMeansClusters 'that'
#if BUILDMODE != 0 
inline 
#endif
const VecFloat* _KMeansClustersCenterFromId(
  const KMeansClusters* const that, const int iCluster);

// Get the seed of the KMeansClusters 'that'
#if BUILDMODE != 0 
inline 
#endif
KMeansClustersSeed KMeansClustersGetSeed(const KMeansClusters* const that);

// Set the seed of the KMeansClusters 'that' to 'seed'
#if BUILDMODE != 0 
inline 
#endif
void KMeansClustersSetSeed(KMeansClusters* const that,
  const KMeansClustersSeed seed);

// Get the center of the cluster including the 'input' data for the 
// KMeansClusters 'that' 
#if BUILDMODE != 0 
inline 
#endif
const VecFloat* _KMeansClustersCenterFromPos(
  const KMeansClusters* const that, const VecFloat* input);

// Get the index of the cluster including the 'input' data for the 
// KMeansClusters 'that' 
int KMeansClustersGetId(const KMeansClusters* const that, 
  const VecFloat* input);

// Get the seed of the KMeansClusters 'that'
#if BUILDMODE != 0 
inline 
#endif
int KMeansClustersGetK(const KMeansClusters* const that);

// Print the KMeansClusters 'that' on the stream 'stream'
void KMeansClustersPrintln(const KMeansClusters* const that,
  FILE* const stream);

// ================= Polymorphism ==================

#define KMeansClustersCenter(Cluster, Input) _Generic(Input, \
  VecFloat*: _KMeansClustersCenterFromPos, \
  const VecFloat*: _KMeansClustersCenterFromPos, \
  int: _KMeansClustersCenterFromId, \
  const int: _KMeansClustersCenterFromId, \
  default: PBErrInvalidPolymorphism)(Cluster, Input)

// ================ Inliner ====================

#if BUILDMODE != 0
#include "pbdataanalysis-inline.c"
#endif

#endif
