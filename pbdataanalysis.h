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
#include "pbjson.h"
#include "gset.h"
#include "gdataset.h"

// ----------------- Principal component analysis ---------------

// ================= Define ==================

// ================= Data structure ===================

typedef struct PrincipalComponentAnalysis {
  // Principal components, stored from the most influent to the less one
  GSetVecFloat _components;
} PrincipalComponentAnalysis;

// ================ Functions declaration ====================

// Create a static PrincipalComponentAnalysis
PrincipalComponentAnalysis PrincipalComponentAnalysisCreateStatic();

// Free the memory used by the static PrincipalComponentAnalysis
void PrincipalComponentAnalysisFreeStatic(
  PrincipalComponentAnalysis* const that);

// Calculate the principal components for the 'dataset' and 
// store the result into the PrincipalComponentAnalysis 'that'
void PCASearch(PrincipalComponentAnalysis* const that, 
  const GDataSetVecFloat* const dataset);
  
// Get the 'dataset' converted through the first 'nb' components of the 
// PrincipalComponentAnalysis 'that'
// Return a new data set, the dataset in arguments is not modified
// Return an empty dataset if the principal components have not yet been
// calculated (using the PCASearch function)
// Returned VecFloat have dimension 'nb'. Returned GSet has same nbsample
// as 'dataset'
// Dimension of VecFloat in 'dataset' must be equal to the size of
// components in 'that'
GDataSetVecFloat PCAConvert(const PrincipalComponentAnalysis* const that,
  const GDataSetVecFloat* const dataset, const int nb);

// Get the principal components of the PrincipalComponentAnalysis 'that'
#if BUILDMODE != 0 
static inline 
#endif
const GSetVecFloat* PCAComponents(
  const PrincipalComponentAnalysis* const that);

// Get the number of principal components of the 
// PrincipalComponentAnalysis 'that'
#if BUILDMODE != 0 
static inline 
#endif
int PCAGetNbComponents(
  const PrincipalComponentAnalysis* const that);
  
// Print the principal components of the PrincipalComponentAnalysis 'that'
// on 'stream'
void PCAPrintln(const PrincipalComponentAnalysis* const that,
  FILE* const stream);
  
// ----------------- K-means clustering ---------------

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

// Free the memory used by a KMeansClusters
void KMeansClustersFreeStatic(KMeansClusters* const that);

// Create the set of 'K' clusters clustering the 'input' data according 
// to the K-means algorithm
// 'K' must be inferior or equal to the number of input data
// srandom() must have been called before using this function
void KMeansClustersSearch(KMeansClusters* const that,
  const GSetVecFloat* const input, const int K);

// Get the set of clusters' center for the KMeansClusters 'that'
#if BUILDMODE != 0 
static inline 
#endif
const GSetVecFloat* KMeansClustersCenters(
  const KMeansClusters* const that);

// Get the 'iCluster'-th cluster's center for the KMeansClusters 'that'
#if BUILDMODE != 0 
static inline 
#endif
const VecFloat* _KMeansClustersCenterFromId(
  const KMeansClusters* const that, const int iCluster);

// Get the seed of the KMeansClusters 'that'
#if BUILDMODE != 0 
static inline 
#endif
KMeansClustersSeed KMeansClustersGetSeed(const KMeansClusters* const that);

// Set the seed of the KMeansClusters 'that' to 'seed'
#if BUILDMODE != 0 
static inline 
#endif
void KMeansClustersSetSeed(KMeansClusters* const that,
  const KMeansClustersSeed seed);

// Get the center of the cluster including the 'input' data for the 
// KMeansClusters 'that' 
#if BUILDMODE != 0 
static inline 
#endif
const VecFloat* _KMeansClustersCenterFromPos(
  const KMeansClusters* const that, const VecFloat* input);

// Get the index of the cluster including the 'input' data for the 
// KMeansClusters 'that' 
int KMeansClustersGetId(const KMeansClusters* const that, 
  const VecFloat* input);

// Get the seed of the KMeansClusters 'that'
#if BUILDMODE != 0 
static inline 
#endif
int KMeansClustersGetK(const KMeansClusters* const that);

// Print the KMeansClusters 'that' on the stream 'stream'
void KMeansClustersPrintln(const KMeansClusters* const that,
  FILE* const stream);

// Load the KMeansClusters 'that' from the stream 'stream'
bool KMeansClustersLoad(KMeansClusters* that, FILE* const stream);

// Save the KMeansClusters 'that' to the stream 'stream'
// If 'compact' equals true it saves in compact form, else it saves in 
// readable form
// Return true upon success else false
bool KMeansClustersSave(const KMeansClusters* const that, 
  FILE* const stream, const bool compact);

// Function which return the JSON encoding of 'that' 
JSONNode* KMeansClustersEncodeAsJSON(const KMeansClusters* const that);

// Function which decode from JSON encoding 'json' to 'that'
bool KMeansClustersDecodeAsJSON(KMeansClusters* that, 
  const JSONNode* const json);

// ================= Polymorphism ==================

#define KMeansClustersCenter(Cluster, Input) _Generic(Input, \
  VecFloat*: _KMeansClustersCenterFromPos, \
  const VecFloat*: _KMeansClustersCenterFromPos, \
  int: _KMeansClustersCenterFromId, \
  const int: _KMeansClustersCenterFromId, \
  default: PBErrInvalidPolymorphism)(Cluster, Input)

// ================ static inliner ====================

#if BUILDMODE != 0
#include "pbdataanalysis-inline.c"
#endif

#endif
