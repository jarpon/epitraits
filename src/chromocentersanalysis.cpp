#include <componentlabelling.h>
#include <dataset.h>
//#include "regionanalysis2.h"
#include <regionanalysis.h>
#include <cmath>

#define TRACE
#include <trace.h>

void chromocentersAnalysis(const VoxelMatrix<float>& originalVoxelMatrix, VoxelMatrix<float>& nucleusMask,
                           VoxelMatrix<float>& ccsMask, const string& filename, const int& numNucleus, int& totalNumCCs,
                           DataSet& nucleiDataset, DataSet& chromocentersDataset, DataSet& individualChromocentersDataset)
{
  RegionAnalysis<float> regionAnalysis;
  regionAnalysis.setRegionMatrix( nucleusMask );
  regionAnalysis.run();

  ComponentLabelling<float> componentLabelling;
  componentLabelling.apply( ccsMask );

  RegionAnalysis<float> regionAnalysisCCs;
  regionAnalysisCCs.setRegionMatrix( ccsMask );
  regionAnalysisCCs.run();

  Vector<float> centroid(3);
  Vertices<float> centroids;
  centroids = regionAnalysisCCs.computeRegionCentroids(REGION_FEATURE_CENTROID, originalVoxelMatrix);
  nucleiDataset.setValue ( "name", numNucleus, filename );//filename
  nucleiDataset.setValue ( "ccsNumber", numNucleus, componentLabelling.getNumLabels() );//number of ccs obtained in the nucleus
  nucleiDataset.setValue ( "ccsVolume", numNucleus, regionAnalysisCCs.computeRegionFeature(REGION_FEATURE_VOLUME,originalVoxelMatrix).sum() );//total ccs volume
  nucleiDataset.setValue ( "ccsRelativeVolume", numNucleus, ( regionAnalysisCCs.computeRegionFeature(REGION_FEATURE_VOLUME,originalVoxelMatrix)[0] / regionAnalysis.computeRegionFeature(REGION_FEATURE_VOLUME,originalVoxelMatrix)[0] ) );//relative volume of ccs regarding the complete nucleus
  nucleiDataset.setValue ( "ccsIntensity", numNucleus, regionAnalysisCCs.computeRegionFeature(REGION_FEATURE_INTENSITY,originalVoxelMatrix).sum() );//total absolute intensity of ccs
  nucleiDataset.setValue ( "ccsIntegratedDensity", numNucleus, regionAnalysisCCs.computeRegionFeature(REGION_FEATURE_INTEGRATED_DENSITY,originalVoxelMatrix).sum() );//integrated density of ccs taking into account real volume
  nucleiDataset.setValue ( "RHF", numNucleus, regionAnalysisCCs.computeRegionFeature(REGION_FEATURE_INTENSITY,originalVoxelMatrix).sum() / regionAnalysis.computeRegionFeature(REGION_FEATURE_INTENSITY,originalVoxelMatrix)[0] );//RHF: rate of heterochromatin (int.Density of ccs/ int.Density of nuclei)

/**///chromocenters individual information
  for (int numCC = 0; numCC < regionAnalysisCCs.numRegions(); numCC++ )
  {
    chromocentersDataset.setValue ( "name", numCC+totalNumCCs, filename );
    chromocentersDataset.setValue ( "idCC", numCC+totalNumCCs, numCC+1 );

    centroid = centroids[numCC];
    chromocentersDataset.setValue ( "centroidCoordX", numCC+totalNumCCs, centroid[X] );
    chromocentersDataset.setValue ( "centroidCoordY", numCC+totalNumCCs, centroid[Y] );
    chromocentersDataset.setValue ( "centroidCoordZ", numCC+totalNumCCs, centroid[Z] );

    float ccVolume = regionAnalysisCCs.computeRegionFeature(REGION_FEATURE_VOLUME,originalVoxelMatrix)[numCC];
    float eqRadius = regionAnalysisCCs.computeRegionFeature(REGION_FEATURE_EQUIVALENT_RADIUS,originalVoxelMatrix)[numCC];
    float ccRelativeVolume = regionAnalysisCCs.computeRegionFeature(REGION_FEATURE_VOLUME,originalVoxelMatrix)[numCC] / regionAnalysisCCs.computeRegionFeature(REGION_FEATURE_VOLUME,originalVoxelMatrix).sum();
    float flatness = regionAnalysisCCs.computeRegionFeature(REGION_FEATURE_FLATNESS,originalVoxelMatrix)[numCC];
    float elongation = regionAnalysisCCs.computeRegionFeature(REGION_FEATURE_ELONGATION,originalVoxelMatrix)[numCC];
    float sphericity = regionAnalysisCCs.computeRegionFeature(REGION_FEATURE_SPHERICITY,originalVoxelMatrix)[numCC];
    float surfaceArea = regionAnalysisCCs.computeRegionFeature(REGION_FEATURE_SURFACE_AREA,originalVoxelMatrix)[numCC];
    float ccsIntensity = regionAnalysisCCs.computeRegionFeature(REGION_FEATURE_INTENSITY,originalVoxelMatrix)[numCC];
    float ccsIntegratedDensity = regionAnalysisCCs.computeRegionFeature(REGION_FEATURE_INTEGRATED_DENSITY,originalVoxelMatrix)[numCC];
    float relativeCCsIntensity = ( regionAnalysisCCs.computeRegionFeature(REGION_FEATURE_INTENSITY,originalVoxelMatrix)[numCC] ) / regionAnalysisCCs.computeRegionFeature(REGION_FEATURE_INTENSITY,originalVoxelMatrix).sum();

    chromocentersDataset.setValue ( "equivalentRadius", numCC+totalNumCCs, eqRadius );
    chromocentersDataset.setValue ( "ccVolume", numCC+totalNumCCs, ccVolume );
    chromocentersDataset.setValue ( "ccRelativeVolume", numCC+totalNumCCs, ccRelativeVolume );
    chromocentersDataset.setValue ( "flatness", numCC+totalNumCCs, flatness );
    chromocentersDataset.setValue ( "elongation", numCC+totalNumCCs, elongation );
    chromocentersDataset.setValue ( "sphericity", numCC+totalNumCCs, sphericity );
    chromocentersDataset.setValue ( "surfaceArea", numCC+totalNumCCs, surfaceArea );
    chromocentersDataset.setValue ( "ccsIntensity", numCC+totalNumCCs, ccsIntensity );
    chromocentersDataset.setValue ( "ccsIntegratedDensity", numCC+totalNumCCs, ccsIntegratedDensity );
    //this last relativeCCsIntensity is the intensity of each cc divided by the total cc's intensity
    chromocentersDataset.setValue ( "relativeCCsIntensity", numCC+totalNumCCs, relativeCCsIntensity );

    individualChromocentersDataset.setValue ( "id", numCC, filename );
    individualChromocentersDataset.setValue ( "idCC", numCC, numCC+1 );
    individualChromocentersDataset.setValue ( "centroidCoordX", numCC, centroid[X] );
    individualChromocentersDataset.setValue ( "centroidCoordY", numCC, centroid[Y] );
    individualChromocentersDataset.setValue ( "centroidCoordZ", numCC, centroid[Z] );
    individualChromocentersDataset.setValue ( "equivalentRadius", numCC+totalNumCCs, eqRadius );
    individualChromocentersDataset.setValue ( "ccVolume", numCC, ccVolume );
    individualChromocentersDataset.setValue ( "ccRelativeVolume", numCC, ccRelativeVolume );
    individualChromocentersDataset.setValue ( "flatness", numCC, flatness );
    individualChromocentersDataset.setValue ( "elongation", numCC, elongation );
    individualChromocentersDataset.setValue ( "sphericity", numCC, sphericity );
    individualChromocentersDataset.setValue ( "surfaceArea", numCC, surfaceArea );
    individualChromocentersDataset.setValue ( "ccsIntensity", numCC, ccsIntensity );
    individualChromocentersDataset.setValue ( "ccsIntegratedDensity", numCC, ccsIntegratedDensity );
    individualChromocentersDataset.setValue ( "relativeCCsIntensity", numCC, relativeCCsIntensity );

  }

  totalNumCCs += regionAnalysisCCs.numRegions();

  //return chromocentersDataset;
}
