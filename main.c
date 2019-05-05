#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "pbdataanalysis.h"

void UnitTestKMean() {
  srandom(3);
  GSetVecFloat input = GSetVecFloatCreateStatic();
  int K = 3;
  float mean = 0.0;
  float sigma = 3.0;
  Gauss gauss = GaussCreateStatic(mean, sigma);
  int nbData = 100;
  VecFloat2D tgtMean[3];
  FILE* csvFile = fopen("./kmeancluster.csv", "w");
  for (int iInput = 0; iInput < K; ++iInput) {
    tgtMean[iInput] = VecFloatCreateStatic2D();
    VecSet(tgtMean + iInput, 0, rnd() * 50.0);
    VecSet(tgtMean + iInput, 1, rnd() * 50.0);
    printf("Target #%d: ", iInput);
    VecPrint(tgtMean + iInput, stdout);
    printf("\n");
    fprintf(csvFile, "%f %f\n", 
      VecGet(tgtMean + iInput, 0),
      VecGet(tgtMean + iInput, 1));
  }
  fprintf(csvFile, "\n"); 
  for (int iData = 0; iData < nbData; ++iData) {
    for (int iInput = 0; iInput < K; ++iInput) {
      VecFloat* vFloat = VecFloatCreate(2);
      GSetAppend(&input, vFloat);
      VecSet(vFloat, 0, VecGet(tgtMean + iInput, 0) + 
        GaussRnd(&gauss));
      VecSet(vFloat, 1, VecGet(tgtMean + iInput, 1) + 
        GaussRnd(&gauss));
      fprintf(csvFile, "%f %f ", 
        VecGet(vFloat, 0), VecGet(vFloat, 1));
    }
    fprintf(csvFile, "\n"); 
  }
  fprintf(csvFile, "\n"); 
  printf("--- Seed: random\n");
  KMeansClustersSeed seed = KMeansClustersSeed_Random;
  KMeansClusters clusters = KMeansClustersCreateStatic(seed);
  KMeansClustersSearch(&clusters, &input, K);
  if (GSetNbElem(KMeansClustersCenters(&clusters)) != K) {
    PBDataAnalysisErr->_type = PBErrTypeUnitTestFailed;
    sprintf(PBDataAnalysisErr->_msg, "_GetKMeansClusterFloat NOK");
    PBErrCatch(PBDataAnalysisErr);
  }
  for (int iCenter = 0; iCenter < K; ++iCenter) {
    const VecFloat* vFloat = KMeansClustersCenter(&clusters, iCenter);
    printf("Cluster #%d: ", iCenter);
    VecPrint(vFloat, stdout);
    printf("\n");
    fprintf(csvFile, "%f %f\n", 
      VecGet(vFloat, 0), VecGet(vFloat, 1));
  }
  KMeansClustersFreeStatic(&clusters);
  printf("--- Seed: forgy\n");
  seed = KMeansClustersSeed_Forgy;
  clusters = KMeansClustersCreateStatic(seed);
  KMeansClustersSearch(&clusters, &input, K);
  if (GSetNbElem(KMeansClustersCenters(&clusters)) != K) {
    PBDataAnalysisErr->_type = PBErrTypeUnitTestFailed;
    sprintf(PBDataAnalysisErr->_msg, "_GetKMeansClusterFloat NOK");
    PBErrCatch(PBDataAnalysisErr);
  }
  for (int iCenter = 0; iCenter < K; ++iCenter) {
    const VecFloat* vFloat = KMeansClustersCenter(&clusters, iCenter);
    printf("Cluster #%d: ", iCenter);
    VecPrint(vFloat, stdout);
    printf("\n");
    fprintf(csvFile, "%f %f\n", 
      VecGet(vFloat, 0), VecGet(vFloat, 1));
  }
  KMeansClustersFreeStatic(&clusters);
  printf("--- Seed: plusplus\n");
  seed = KMeansClustersSeed_PlusPlus;
  clusters = KMeansClustersCreateStatic(seed);
  KMeansClustersSearch(&clusters, &input, K);
  if (GSetNbElem(KMeansClustersCenters(&clusters)) != K) {
    PBDataAnalysisErr->_type = PBErrTypeUnitTestFailed;
    sprintf(PBDataAnalysisErr->_msg, "_GetKMeansClusterFloat NOK");
    PBErrCatch(PBDataAnalysisErr);
  }
  for (int iCenter = 0; iCenter < K; ++iCenter) {
    const VecFloat* vFloat = KMeansClustersCenter(&clusters, iCenter);
    printf("Cluster #%d: ", iCenter);
    VecPrint(vFloat, stdout);
    printf("\n");
    fprintf(csvFile, "%f %f\n", 
      VecGet(vFloat, 0), VecGet(vFloat, 1));
  }
  FILE* fd = fopen("./kmeancluster.txt", "w");
  if (!KMeansClustersSave(&clusters, fd, false)) {
    PBDataAnalysisErr->_type = PBErrTypeUnitTestFailed;
    sprintf(PBDataAnalysisErr->_msg, "KMeansClustersSave NOK");
    PBErrCatch(PBDataAnalysisErr);
  }
  fclose(fd);
  KMeansClusters loadClusters = 
    KMeansClustersCreateStatic(KMeansClustersSeed_Default);
  fd = fopen("./kmeancluster.txt", "r");
  if (!KMeansClustersLoad(&loadClusters, fd)) {
    PBDataAnalysisErr->_type = PBErrTypeUnitTestFailed;
    sprintf(PBDataAnalysisErr->_msg, "KMeansClustersLoad NOK");
    PBErrCatch(PBDataAnalysisErr);
  }
  fclose(fd);
  if (clusters._seed != loadClusters._seed ||
    GSetNbElem(KMeansClustersCenters(&clusters)) != 
      GSetNbElem(KMeansClustersCenters(&loadClusters))) {
    PBDataAnalysisErr->_type = PBErrTypeUnitTestFailed;
    sprintf(PBDataAnalysisErr->_msg, "KMeansClustersLoad NOK");
    PBErrCatch(PBDataAnalysisErr);
  }
  GSetIterForward iter = 
    GSetIterForwardCreateStatic(KMeansClustersCenters(&clusters));
  GSetIterForward iterLoad = 
    GSetIterForwardCreateStatic(KMeansClustersCenters(&loadClusters));
  do {
    VecFloat* u = GSetIterGet(&iter);
    VecFloat* v = GSetIterGet(&iterLoad);
    if (VecIsEqual(u, v) == false) {
      PBDataAnalysisErr->_type = PBErrTypeUnitTestFailed;
      sprintf(PBDataAnalysisErr->_msg, "KMeansClustersLoad NOK");
      PBErrCatch(PBDataAnalysisErr);
    }
  } while (GSetIterStep(&iter) && GSetIterStep(&iterLoad));
  KMeansClustersFreeStatic(&loadClusters);
  KMeansClustersFreeStatic(&clusters);
  fclose(csvFile);
  while (GSetNbElem(&input) > 0) {
    VecFloat* v = GSetPop(&input);
    VecFree(&v);
  }
  printf("UnitTestKMean OK\n");
}

void UnitTestPCA() {
  char* datasetPath = "./unitTestPCA.json";
  GDataSetVecFloat dataset = 
    GDataSetVecFloatCreateStaticFromFile(datasetPath);
  PrincipalComponentAnalysis pca =
    PrincipalComponentAnalysisCreateStatic();
  PCASearch(&pca, &dataset);
  printf("Components:\n");
  PCAPrintln(&pca, stdout);
  GDataSetVecFloat convDataSet1D = PCAConvert(&pca, &dataset, 1);
  printf("Projection on 1st component:\n");
  GSetIterForward iter = 
    GSetIterForwardCreateStatic(GDSSamples(&convDataSet1D));
  do {
    VecFloat* sample = GSetIterGet(&iter);
    VecPrint(sample, stdout); printf("\n");
  } while (GSetIterStep(&iter));
  GDataSetVecFloat convDataSet2D = PCAConvert(&pca, &dataset, 2);
  printf("Projection on 2nd component:\n");
  iter = GSetIterForwardCreateStatic(GDSSamples(&convDataSet2D));
  do {
    VecFloat* sample = GSetIterGet(&iter);
    VecPrint(sample, stdout); printf("\n");
  } while (GSetIterStep(&iter));
  GDataSetVecFloatFreeStatic(&convDataSet1D);
  GDataSetVecFloatFreeStatic(&convDataSet2D);
  PrincipalComponentAnalysisFreeStatic(&pca);
  GDataSetVecFloatFreeStatic(&dataset);
  printf("UnitTestPCA OK\n");
}

void UnitTestAll() {
  UnitTestKMean();
  UnitTestPCA();
}

int main(void) {
  UnitTestAll();
  return 0;
}

