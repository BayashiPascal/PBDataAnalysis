#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "pbdataanalysis.h"

float* GetYoloBoxes(int K, int nbData, float* data) {

  // Init the random generator
  srandom(time(NULL));

  // Create the set of input data 
  GSetVecFloat input = GSetVecFloatCreateStatic();
  for (int iData = 0; iData < nbData; ++iData) {
    VecFloat* vFloat = VecFloatCreate(2);
    GSetAppend(&input, vFloat);
    for (int iInput = 0; iInput < 2; ++iInput) {
      VecSet(vFloat, iInput, data[iData * 2 + iInput]);
    }
  }

  // Create the KMeansClusters
  KMeansClustersSeed seed = KMeansClustersSeed_Forgy;
  KMeansClusters clusters = KMeansClustersCreateStatic(seed);

  // Search the K-Means
  KMeansClustersSearch(&clusters, &input, K);

  // Store the found K-Means
  float* ret = malloc(sizeof(float) * K * 2);
  if (data == NULL) {
    printf("Failed to allocate memory\n");
    exit(1);
  }
  for (int iCenter = 0; iCenter < K; ++iCenter) {
    const VecFloat* vFloat = KMeansClustersCenter(&clusters, iCenter);
    ret[iCenter * 2] = VecGet(vFloat, 0);
    ret[iCenter * 2 + 1] = VecGet(vFloat, 1);
  }

  // Free memory
  KMeansClustersFreeStatic(&clusters);
  while (GSetNbElem(&input) > 0) {
    VecFloat* v = GSetPop(&input);
    VecFree(&v);
  }
  
  // Return the result
  return ret;
}

// Arguments:
// main <width image> <height image> <K, equals to 'num' in the 
// yolo config file> <input file>
// where input file has the following format: first line is the 
// number of boxes in input, following lines are relative width
// and relative height separated by a space of each box, one per line

int main(int argc, char** argv) {

  // Check the number of arguments
  if (argc != 5) {
    printf("Invalid number of arguments (%d)\n", argc);
    exit(1);
  }
  
  // Get the image dimension
  int width = atoi(argv[1]);
  int height = atoi(argv[2]);
  
  // Get the number of clusters
  int K = atoi(argv[3]);
  
  // Get the data file name
  char* filename = argv[4];
  FILE* fp = fopen(filename, "r");
  
  // Read the number of data
  int nbData;
  if (fscanf(fp, "%d\n", &nbData) == EOF) {
    printf("Failed to read the data\n");
    exit(1);
  }
  
  // Get the data
  float* data = malloc(sizeof(float) * (nbData * 2));
  if (data == NULL) {
    printf("Failed to allocate memory\n");
    exit(1);
  }
  for (int iData = 0; iData < nbData; ++iData) {
    if (fscanf(fp, "%f %f\n", data + iData * 2, data + iData * 2 + 1) == EOF) {
      printf("Failed to read the data\n");
      exit(1);
    }
  }
  
  // Close the input file
  fclose(fp);

  // Search the clusters
  float* ret = GetYoloBoxes(K, nbData, data);
  
  // Display the results
  for (int i = 0; i < K; ++i) {
    int w = (int)round(ret[2 * i] * (float)width);
    int h = (int)round(ret[2 * i + 1] * (float)height);
    printf("%d,%d ", w, h);
  }
  printf("\n");
  
  // Free memory
  free(data);
  free(ret);

  return 0;
}

