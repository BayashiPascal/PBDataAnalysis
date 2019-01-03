// ============ PBDATAANALYSIS_INLINE.C ================

// ================ Functions implementation ====================

// Get the set of clusters' center for the KMeansClusters 'that'
#if BUILDMODE != 0 
inline 
#endif
const GSetVecFloat* KMeansClustersCenters(
  const KMeansClusters* const that) {
#if BUILDMODE == 0
  if (that == NULL) {
    PBDataAnalysisErr->_type = PBErrTypeNullPointer;
    sprintf(PBDataAnalysisErr->_msg, "'that' is null");
    PBErrCatch(GSetErr);
  }
#endif
  return &(that->_centers);
}

// Get the 'iCluster'-th cluster's center for the KMeansClusters 'that'
#if BUILDMODE != 0 
inline 
#endif
const VecFloat* _KMeansClustersCenterFromId(
  const KMeansClusters* const that, const int iCluster) {
#if BUILDMODE == 0
  if (that == NULL) {
    PBDataAnalysisErr->_type = PBErrTypeNullPointer;
    sprintf(PBDataAnalysisErr->_msg, "'that' is null");
    PBErrCatch(GSetErr);
  }
  if (iCluster < 0 || iCluster >= GSetNbElem(&(that->_centers))) {
    PBDataAnalysisErr->_type = PBErrTypeInvalidArg;
    sprintf(PBDataAnalysisErr->_msg, "'iCluster' is invalid (0<=%d<%ld)",
      iCluster, GSetNbElem(&(that->_centers)));
    PBErrCatch(GSetErr);
  }
#endif
  return GSetGet(&(that->_centers), iCluster);
}

// Get the seed of the KMeansClusters 'that'
#if BUILDMODE != 0 
inline 
#endif
KMeansClustersSeed KMeansClustersGetSeed(
  const KMeansClusters* const that) {
#if BUILDMODE == 0
  if (that == NULL) {
    PBDataAnalysisErr->_type = PBErrTypeNullPointer;
    sprintf(PBDataAnalysisErr->_msg, "'that' is null");
    PBErrCatch(GSetErr);
  }
#endif
  return that->_seed;
}

// Set the seed of the KMeansClusters 'that' to 'seed'
#if BUILDMODE != 0 
inline 
#endif
void KMeansClustersSetSeed(KMeansClusters* const that,
  const KMeansClustersSeed seed) {
#if BUILDMODE == 0
  if (that == NULL) {
    PBDataAnalysisErr->_type = PBErrTypeNullPointer;
    sprintf(PBDataAnalysisErr->_msg, "'that' is null");
    PBErrCatch(GSetErr);
  }
#endif
  that->_seed = seed;
}

// Get the seed of the KMeansClusters 'that'
#if BUILDMODE != 0 
inline 
#endif
int KMeansClustersGetK(const KMeansClusters* const that) {
#if BUILDMODE == 0
  if (that == NULL) {
    PBDataAnalysisErr->_type = PBErrTypeNullPointer;
    sprintf(PBDataAnalysisErr->_msg, "'that' is null");
    PBErrCatch(GSetErr);
  }
#endif
  return GSetNbElem(&(that->_centers));
}

// Get the center of the cluster including the 'input' data for the 
// KMeansClusters 'that' 
#if BUILDMODE != 0 
inline 
#endif
const VecFloat* _KMeansClustersCenterFromPos(
  const KMeansClusters* const that, const VecFloat* input) {
#if BUILDMODE == 0
  if (that == NULL) {
    PBDataAnalysisErr->_type = PBErrTypeNullPointer;
    sprintf(PBDataAnalysisErr->_msg, "'that' is null");
    PBErrCatch(GSetErr);
  }
  if (input == NULL) {
    PBDataAnalysisErr->_type = PBErrTypeNullPointer;
    sprintf(PBDataAnalysisErr->_msg, "'input' is null");
    PBErrCatch(GSetErr);
  }
#endif  
  return KMeansClustersCenter(that, KMeansClustersGetId(that, input));
}

