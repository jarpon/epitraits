#include <componentlabelling.h>
#include <dataset.h>
//#include "regionanalysis2.h"
#include <regionanalysis.h>
#include <cmath>
#include <marchingcubes.h>
#include <thresholding.h>
#include <trimesh.h>

#define TRACE
#include <trace.h>

void chromocentersAnalysis(VoxelMatrix<float>& ccsMask, const string& filename, const string& parentDir,
                           const int& numNucleus, int& totalNumCCs,
                           DataSet& nucleiDataset, DataSet& chromocentersDataset, DataSet& individualChromocentersDataset)
{
  VoxelMatrix<float> originalVoxelMatrix( parentDir + "/originals_vm/" + filename + ".vm" );
  VoxelMatrix<float> nucleusMask( parentDir + "/segmented_nuclei/" + filename + ".vm" );
  RegionAnalysis<float> regionAnalysis;
  regionAnalysis.setRegionMatrix( nucleusMask );
  regionAnalysis.run();

  ComponentLabelling<float> componentLabelling;
  componentLabelling.apply( ccsMask );

  RegionAnalysis<float> regionAnalysisCCs;
  regionAnalysisCCs.setRegionMatrix( ccsMask );
  regionAnalysisCCs.run();

  //TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + "_nucleus.tm" );
  TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + ".tm" );

  Vector<float> centroid(3);
  Vector<float> vertexTriMesh(3);

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

    MarchingCubes<float> marchingCubes;
    TriMesh<float> triMesh;
    VoxelMatrix<float> currentLabeledVM = ccsMask;
    Thresholding<float> thresholding;
    thresholding.setForeground( 1.0 );
    thresholding.setBackground( 0.0 );
    thresholding.levelSetMask( currentLabeledVM, numCC+1 );
    triMesh = marchingCubes.buildMesh( currentLabeledVM, 0.5, true );
    triMesh.scale( originalVoxelMatrix.getVoxelCalibration().getVoxelSize() );

    nucleusTriMesh.closestPoint( centroid, vertexTriMesh );
    float distanceToBorder = centroid.distance( vertexTriMesh );
    float ccVolume_tm = fabs(triMesh.volume());
    float eqRadius_tm = triMesh.equivalentRadius();

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

    chromocentersDataset.setValue ( "equivalentRadius_vm", numCC+totalNumCCs, eqRadius );
    chromocentersDataset.setValue ( "equivalentRadius_tm", numCC+totalNumCCs, eqRadius_tm );
    chromocentersDataset.setValue ( "ccVolume_vm", numCC+totalNumCCs, ccVolume );
    chromocentersDataset.setValue ( "ccVolume_tm", numCC+totalNumCCs, ccVolume_tm );
    chromocentersDataset.setValue ( "distanceToTheBorder", numCC+totalNumCCs, distanceToBorder );
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
    individualChromocentersDataset.setValue ( "equivalentRadius_vm", numCC, eqRadius );
    individualChromocentersDataset.setValue ( "equivalentRadius_tm", numCC, eqRadius_tm );
    individualChromocentersDataset.setValue ( "distanceToTheBorder", numCC, distanceToBorder );
    individualChromocentersDataset.setValue ( "ccVolume_vm", numCC, ccVolume );
    individualChromocentersDataset.setValue ( "ccVolume_tm", numCC, ccVolume_tm );
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
