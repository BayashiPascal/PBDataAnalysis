// ============ PBDATAANALYSIS.C ================

// ================= Include =================

#include "pbdataanalysis.h"
#if BUILDMODE == 0
#include "pbdataanalysis-inline.c"
#endif

// ================= Define ==================

// ================ Functions declaration ====================

void KMeansClustersInitRandom(KMeansClusters* const that,
  const GSetVecFloat* const input);

void KMeansClustersInitForgy(KMeansClusters* const that,
  const GSetVecFloat* const input);

// ================ Functions implementation ====================

// Create a KMeansClusters for a K-means clustering initialized 
// using the 'seed' technique
// srandom() must have been called before using this function
KMeansClusters KMeansClustersCreateStatic(const KMeansClustersSeed seed) {
  // Declare the KMeansClusters
  KMeansClusters clusters;
  // Init the properties
  clusters._seed = seed;
  clusters._centers = GSetVecFloatCreateStatic();
  // Return the KMeansClusters
  return clusters;
}

// Free the memory used by a PBDAKMeansClusters
void KMeansClustersFreeStatic(KMeansClusters* const that) {
#if BUILDMODE == 0
  if (that == NULL) {
    PBDataAnalysisErr->_type = PBErrTypeNullPointer;
    sprintf(PBDataAnalysisErr->_msg, "'that' is null");
    PBErrCatch(PBDataAnalysisErr);
  }
#endif
  // Free the memory used by the cluster centers
  while (GSetNbElem(&(that->_centers)) > 0) {
    VecFloat* v = GSetPop((GSetVecFloat*)(KMeansClustersCenters(that)));
    free(v);
  }
}

// Create the set of 'K' clusters clustering the 'input' data according 
// to the K-means algorithm
// 'K' must be inferior or equal to the number of input data
// srandom() must have been called before using this function
void KMeansClustersSearch(KMeansClusters* const that, 
  const GSetVecFloat* const input, const int K) {
#if BUILDMODE == 0
  if (that == NULL) {
    PBDataAnalysisErr->_type = PBErrTypeNullPointer;
    sprintf(PBDataAnalysisErr->_msg, "'that' is null");
    PBErrCatch(PBDataAnalysisErr);
  }
  if (input == NULL) {
    PBDataAnalysisErr->_type = PBErrTypeNullPointer;
    sprintf(PBDataAnalysisErr->_msg, "'input' is null");
    PBErrCatch(PBDataAnalysisErr);
  }
  if (K < 1) {
    PBDataAnalysisErr->_type = PBErrTypeInvalidArg;
    sprintf(PBDataAnalysisErr->_msg, "'K' is invalid (%d>=1)", K);
    PBErrCatch(PBDataAnalysisErr);
  }
  if (GSetNbElem(input) < K) {
    PBDataAnalysisErr->_type = PBErrTypeInvalidArg;
    sprintf(PBDataAnalysisErr->_msg, "'K' is invalid (%d<=%ld)", K,
      GSetNbElem(input));
    PBErrCatch(PBDataAnalysisErr);
  }
#endif
  // If there are alreay computed clusters,
  // free the memory used by the cluster centers
  while (GSetNbElem(&(that->_centers)) > 0) {
    VecFloat* v = GSetPop((GSetVecFloat*)(KMeansClustersCenters(that)));
    free(v);
  }
  // Allocate an array of sets used for computation
  GSetVecFloat* inputsByCluster = 
    PBErrMalloc(PBDataAnalysisErr, sizeof(GSetVecFloat) * K);
  // Get the dimension of the input
  long dim = VecGetDim(GSetGet(input, 0));
  // Allocate memory for the clusters' center and sets used for 
  // computation
  for (int iCenter = K; iCenter--;) {
    // The dimension of the means is the same as the one of the input 
    // data
    VecFloat* v = VecFloatCreate(dim);
    GSetAppend((GSetVecFloat*)(KMeansClustersCenters(that)), v);
    inputsByCluster[iCenter] = GSetVecFloatCreateStatic();
  }
  // Initialise the means according to the seed
  switch(KMeansClustersGetSeed(that)) {
    case KMeansClustersSeed_Random:
      KMeansClustersInitRandom(that, input);
      break;
    case KMeansClustersSeed_Forgy:
      KMeansClustersInitForgy(that, input);
      break;
    default:
      break;
  }
  // Create a vector used for computation
  VecFloat* w = VecFloatCreate(dim);
  // Loop until the clusters' center are aligned with their
  // real center according to input data
  float shift = 1.0;
  while (shift > PBMATH_EPSILON) {
    // Reset the shift
    shift = 0.0;
    // Loop on input data
    GSetIterForward iterInput = GSetIterForwardCreateStatic(input);
    do {
      // Get the input data
      VecFloat* v = GSetIterGet(&iterInput);
      // Get the id of the cluster containing this data
      int clusterId = KMeansClustersGetId(that, v);
      // Add the input to the corresponding set
      GSetAppend(inputsByCluster + clusterId, v);
    } while(GSetIterStep(&iterInput));
    // For each cluster
    for (int iCenter = K; iCenter--;) {
      // Reset the vector used for computation
      VecSetNull(w);
      // Memorize the number of input associated to this cluster
      int nbInput = GSetNbElem(inputsByCluster + iCenter);
      // If there are data in this cluster
      if (nbInput > 0) {
        // Loop on the input contained by this cluster
        GSetIterForward iterSet = 
          GSetIterForwardCreateStatic(inputsByCluster + iCenter);
        do {
          // Get the input data
          VecFloat* v = GSetIterGet(&iterSet);
          // Sum it to the temporary vector
          VecOp(w, 1.0, v, 1.0); 
        } while(GSetIterStep(&iterSet));
        // Get the center (average) of the input contained by this 
        // cluster 
        VecScale(w, 1.0 / (float)nbInput);
        // Update the shift
        shift += VecDist(KMeansClustersCenter(that, iCenter), w);
        // Update the cluster center with the center of the input
//VecPrint(KMeansClustersCenter(that, iCenter), stdout);printf(" ");
        VecCopy((VecFloat*)KMeansClustersCenter(that, iCenter), w);
        // Reset the sets of input data
        GSetFlush(inputsByCluster + iCenter);
      }
    }
//printf("\n");
  }
  // Free the memory used by the vector and sets used for computation
  free(inputsByCluster);
  VecFree(&w);
}

void KMeansClustersInitRandom(KMeansClusters* const that,
  const GSetVecFloat* const input) {
#if BUILDMODE == 0
  if (that == NULL) {
    PBDataAnalysisErr->_type = PBErrTypeNullPointer;
    sprintf(PBDataAnalysisErr->_msg, "'that' is null");
    PBErrCatch(PBDataAnalysisErr);
  }
  if (input == NULL) {
    PBDataAnalysisErr->_type = PBErrTypeNullPointer;
    sprintf(PBDataAnalysisErr->_msg, "'input' is null");
    PBErrCatch(PBDataAnalysisErr);
  }
#endif
  // Get the bounds of the input data
  GSetVecFloat bounds = GSetGetBounds(input); 
  VecFloat* boundMin = GSetGet(&bounds, 0);
  VecFloat* boundMax = GSetGet(&bounds, 1);
  // Get the dimension of the data
  int dim = VecGetDim(GSetGet(input, 0));
  // For each cluster's center
  for (int iCenter = KMeansClustersGetK(that); iCenter--;) {
    // Get the iCenter-th cluster center
    VecFloat* center = (VecFloat*)KMeansClustersCenter(that, iCenter);
    // Initialize randomly the components of the iCenter-th cluster's 
    // center
    for (int iDim = dim; iDim--;) {
      VecSet(center, iDim, VecGet(boundMin, iDim) + 
        rnd() * (VecGet(boundMax, iDim) - VecGet(boundMin, iDim)));
    }
  }
  // Free memory
  GSetFlush(&bounds);
  VecFree(&boundMin);
  VecFree(&boundMax);
}

void KMeansClustersInitForgy(KMeansClusters* const that,
  const GSetVecFloat* const input) {
#if BUILDMODE == 0
  if (that == NULL) {
    PBDataAnalysisErr->_type = PBErrTypeNullPointer;
    sprintf(PBDataAnalysisErr->_msg, "'that' is null");
    PBErrCatch(PBDataAnalysisErr);
  }
  if (input == NULL) {
    PBDataAnalysisErr->_type = PBErrTypeNullPointer;
    sprintf(PBDataAnalysisErr->_msg, "'input' is null");
    PBErrCatch(PBDataAnalysisErr);
  }
#endif
  // Create a set to select the seeds
  GSetVecFloat forgyInp = GSetVecFloatCreateStatic();
  // Get the number of inputs
  long nbInput = GSetNbElem(input);
  // Declare a variable to pick one input randomly
  VecFloat* inp = NULL;
  // For each cluster's center
  for (int iCenter = KMeansClustersGetK(that); iCenter--;) {
    // Pick an input while avoiding picking twice the same
    // Supposes K is far less than the number of inputs 
    do {
      inp = GSetGet(input, (long)round(rnd() * (nbInput - 1)));
    } while (GSetFirstElem(&forgyInp, inp) != NULL);
    // Set the center of the iCenter-th cluster to this input
    VecCopy((VecFloat*)KMeansClustersCenter(that, iCenter), inp);
  }
  // Empty the set used to select the seeds
  GSetFlush(&forgyInp);
}

// Get the index of the cluster including the 'input' data for the 
// KMeansClusters 'that' 
int KMeansClustersGetId(const KMeansClusters* const that, 
  const VecFloat* input) {
#if BUILDMODE == 0
  if (that == NULL) {
    PBDataAnalysisErr->_type = PBErrTypeNullPointer;
    sprintf(PBDataAnalysisErr->_msg, "'that' is null");
    PBErrCatch(PBDataAnalysisErr);
  }
  if (input == NULL) {
    PBDataAnalysisErr->_type = PBErrTypeNullPointer;
    sprintf(PBDataAnalysisErr->_msg, "'input' is null");
    PBErrCatch(PBDataAnalysisErr);
  }
#endif  
  // Declare variables to memorize the id and distance to the input
  float dist = 0.0;
  int id = -1;
  int iCenter = 0;
  // Loop on the clusters' center
  GSetIterForward iter = 
    GSetIterForwardCreateStatic(KMeansClustersCenters(that));
  do {
    // Get the center of the cluster
    VecFloat* center = GSetIterGet(&iter);
    // Calculate the distance to the input
    float d = VecDist(center, input);
    // If it's the first considered cluster or 
    // if the distance is nearer
    if (id == -1 || dist > d) {
      id = iCenter;
      dist = d;
      // TODO: we can stop if the distance is less than half the 
      // shortest distance between two clusters' center to make
      // this function faster
    }
    // Increment the center index
    ++iCenter;
  } while (GSetIterStep(&iter));
  // Return the id
  return id;
}