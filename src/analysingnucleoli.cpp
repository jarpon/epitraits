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

void nucleoliAnalysis(VoxelMatrix<float>& nucleoliMask, const string& filename, const string& parentDir,
                           int& totalNumNucleoli,
                           DataSet& nucleoliDataset, DataSet& individualNucleoliDataset)
{
  string originalName = filename.substr( 0,filename.find_last_of("-")  );

  VoxelMatrix<float> originalVoxelMatrix( parentDir + "/originals_vm/" + originalName + ".vm" );
  VoxelMatrix<float> nucleusMask( parentDir + "/segmented_nuclei/" + filename + ".vm" );

  // for segmented masks with values higher than 1, we need them with 0/1 values
  if ( nucleusMask.max().max().max() > 1 )
  {
      Thresholding<float> thresholding;
      thresholding.setBackground( 0.0 );
      thresholding.setForeground( 1.0 );
      thresholding.apply( nucleusMask );
  }

  RegionAnalysis<float> regionAnalysis;
  regionAnalysis.setRegionMatrix( nucleusMask );
  regionAnalysis.run();

  ComponentLabelling<float> componentLabelling;
  componentLabelling.apply( nucleoliMask );

  RegionAnalysis<float> regionAnalysisNucleoli;
  regionAnalysisNucleoli.setRegionMatrix( nucleoliMask );
  regionAnalysisNucleoli.run();
  //TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + "_nucleus.tm" );
  //TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + ".tm" );
  TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/nuclei/" + filename + ".tm" );
  Vector<float> centroid(3);
  Vector<float> vertexTriMesh(3);

  // get the name of the class
  string classif = parentDir;
  classif = classif.substr(classif.find_last_of("/\\")+1,classif.length());
  Vertices<float> centroids;
  centroids = regionAnalysisNucleoli.computeRegionCentroids(REGION_FEATURE_CENTROID, originalVoxelMatrix);
//  nucleiDataset.setValue ( "name", numNucleus, filename );//filename
//  nucleiDataset.setValue ( "class", numNucleus, classif );//classification: mutant, tissue, etc.
//  nucleiDataset.setValue ( "ccsNumber", numNucleus, componentLabelling.getNumLabels() );//number of ccs obtained in the nucleus
//  nucleiDataset.setValue ( "ccsVolume", numNucleus, regionAnalysisNucleoli.computeRegionFeature(REGION_FEATURE_VOLUME,originalVoxelMatrix).sum() );//total ccs volume
//  nucleiDataset.setValue ( "ccsRelativeVolume", numNucleus, ( regionAnalysisNucleoli.computeRegionFeature(REGION_FEATURE_VOLUME,originalVoxelMatrix)[0] / regionAnalysis.computeRegionFeature(REGION_FEATURE_VOLUME,originalVoxelMatrix)[0] ) );//relative volume of ccs regarding the complete nucleus
//  nucleiDataset.setValue ( "ccsIntensity", numNucleus, regionAnalysisNucleoli.computeRegionFeature(REGION_FEATURE_INTENSITY,originalVoxelMatrix).sum() );//total absolute intensity of ccs
//  nucleiDataset.setValue ( "ccsIntegratedDensity", numNucleus, regionAnalysisNucleoli.computeRegionFeature(REGION_FEATURE_INTEGRATED_DENSITY,originalVoxelMatrix).sum() );//integrated density of ccs taking into account real volume
//  nucleiDataset.setValue ( "RHF", numNucleus, regionAnalysisNucleoli.computeRegionFeature(REGION_FEATURE_INTENSITY,originalVoxelMatrix).sum() / regionAnalysis.computeRegionFeature(REGION_FEATURE_INTENSITY,originalVoxelMatrix)[0] );//RHF: rate of heterochromatin (int.Density of ccs/ int.Density of nuclei)

/**///chromocenters individual information
  for (int numNucleolus = 0; numNucleolus < regionAnalysisNucleoli.numRegions(); numNucleolus++ )
  {

    centroid = centroids[numNucleolus];

    MarchingCubes<float> marchingCubes;
    TriMesh<float> triMesh;
    VoxelMatrix<float> currentLabeledVM = nucleoliMask;
    Thresholding<float> thresholding;
    thresholding.setForeground( 1.0 );
    thresholding.setBackground( 0.0 );
    thresholding.levelSetMask( currentLabeledVM, numNucleolus+1 );
    triMesh = marchingCubes.buildMesh( currentLabeledVM, 0.5, true );
    triMesh.scale( originalVoxelMatrix.getVoxelCalibration().getVoxelSize() );

    nucleusTriMesh.closestPoint( centroid, vertexTriMesh );
    float distanceToBorder = centroid.distance( vertexTriMesh );
    float nucleolusVolume_tm = fabs(triMesh.volume());
    float eqRadius_tm = triMesh.equivalentRadius();
    float nucleolusVolume = regionAnalysisNucleoli.computeRegionFeature(REGION_FEATURE_VOLUME,originalVoxelMatrix)[numNucleolus];
    float eqRadius = regionAnalysisNucleoli.computeRegionFeature(REGION_FEATURE_EQUIVALENT_RADIUS,originalVoxelMatrix)[numNucleolus];
    float nucleolusRelativeVolume = regionAnalysisNucleoli.computeRegionFeature(REGION_FEATURE_VOLUME,originalVoxelMatrix)[numNucleolus] / regionAnalysisNucleoli.computeRegionFeature(REGION_FEATURE_VOLUME,originalVoxelMatrix).sum();
    float flatness = regionAnalysisNucleoli.computeRegionFeature(REGION_FEATURE_FLATNESS,originalVoxelMatrix)[numNucleolus];
    float elongation = regionAnalysisNucleoli.computeRegionFeature(REGION_FEATURE_ELONGATION,originalVoxelMatrix)[numNucleolus];
    float sphericity = regionAnalysisNucleoli.computeRegionFeature(REGION_FEATURE_SPHERICITY,originalVoxelMatrix)[numNucleolus];
    float sphericity_tm =  36 * M_PI * pow(abs(triMesh.volume()) , 2) / pow( abs(triMesh.area() ), 3);
    float surfaceArea = regionAnalysisNucleoli.computeRegionFeature(REGION_FEATURE_SURFACE_AREA,originalVoxelMatrix)[numNucleolus];
//    float ccsIntensity = regionAnalysisNucleoli.computeRegionFeature(REGION_FEATURE_INTENSITY,originalVoxelMatrix)[numNucleolus];
//    float ccsIntegratedDensity = regionAnalysisNucleoli.computeRegionFeature(REGION_FEATURE_INTEGRATED_DENSITY,originalVoxelMatrix)[numNucleolus];
//    float relativeCCsIntensity = ( regionAnalysisNucleoli.computeRegionFeature(REGION_FEATURE_INTENSITY,originalVoxelMatrix)[numNucleolus] ) / regionAnalysisNucleoli.computeRegionFeature(REGION_FEATURE_INTENSITY,originalVoxelMatrix).sum();
    nucleoliDataset.setValue ( "name", numNucleolus+totalNumNucleoli, filename );
    nucleoliDataset.setValue ( "class", numNucleolus+totalNumNucleoli, classif );//classification: mutant, tissue, etc.
    nucleoliDataset.setValue ( "idCC", numNucleolus+totalNumNucleoli, numNucleolus+1 );
    nucleoliDataset.setValue ( "centroidCoordX", numNucleolus+totalNumNucleoli, centroid[X] );
    nucleoliDataset.setValue ( "centroidCoordY", numNucleolus+totalNumNucleoli, centroid[Y] );
    nucleoliDataset.setValue ( "centroidCoordZ", numNucleolus+totalNumNucleoli, centroid[Z] );
    nucleoliDataset.setValue ( "distanceToTheBorder", numNucleolus+totalNumNucleoli, distanceToBorder );
    nucleoliDataset.setValue ( "equivalentRadius_vm", numNucleolus+totalNumNucleoli, eqRadius );
    nucleoliDataset.setValue ( "equivalentRadius_tm", numNucleolus+totalNumNucleoli, eqRadius_tm );
    nucleoliDataset.setValue ( "nucleolusVolume_vm", numNucleolus+totalNumNucleoli, nucleolusVolume );
    nucleoliDataset.setValue ( "nucleolusVolume_tm", numNucleolus+totalNumNucleoli, nucleolusVolume_tm );
    nucleoliDataset.setValue ( "nucleolusRelativeVolume", numNucleolus+totalNumNucleoli, nucleolusRelativeVolume );
    nucleoliDataset.setValue ( "surfaceArea", numNucleolus+totalNumNucleoli, surfaceArea );
    nucleoliDataset.setValue ( "flatness", numNucleolus+totalNumNucleoli, flatness );
    nucleoliDataset.setValue ( "elongation", numNucleolus+totalNumNucleoli, elongation );
    nucleoliDataset.setValue ( "sphericity_vm", numNucleolus+totalNumNucleoli, sphericity );
    nucleoliDataset.setValue ( "sphericity_tm", numNucleolus+totalNumNucleoli, sphericity_tm );

//    nucleoliDataset.setValue ( "ccsIntensity", numNucleolus+totalNumNucleoli, ccsIntensity );
//    nucleoliDataset.setValue ( "ccsIntegratedDensity", numNucleolus+totalNumNucleoli, ccsIntegratedDensity );
//    //this last relativeCCsIntensity is the intensity of each cc divided by the total cc's intensity
//    nucleoliDataset.setValue ( "relativeCCsIntensity", numNucleolus+totalNumNucleoli, relativeCCsIntensity );
    individualNucleoliDataset.setValue ( "name", numNucleolus, filename );
    individualNucleoliDataset.setValue ( "class", numNucleolus, classif );//classification: mutant, tissue, etc.
    individualNucleoliDataset.setValue ( "idCC", numNucleolus, numNucleolus+1 );
    individualNucleoliDataset.setValue ( "centroidCoordX", numNucleolus, centroid[X] );
    individualNucleoliDataset.setValue ( "centroidCoordY", numNucleolus, centroid[Y] );
    individualNucleoliDataset.setValue ( "centroidCoordZ", numNucleolus, centroid[Z] );
    individualNucleoliDataset.setValue ( "distanceToTheBorder", numNucleolus, distanceToBorder );
    individualNucleoliDataset.setValue ( "equivalentRadius_vm", numNucleolus, eqRadius );
    individualNucleoliDataset.setValue ( "equivalentRadius_tm", numNucleolus, eqRadius_tm );
    individualNucleoliDataset.setValue ( "nucleolusVolume_vm", numNucleolus, nucleolusVolume );
    individualNucleoliDataset.setValue ( "nucleolusVolume_tm", numNucleolus, nucleolusVolume_tm );
    individualNucleoliDataset.setValue ( "nucleolusRelativeVolume", numNucleolus, nucleolusRelativeVolume );
    individualNucleoliDataset.setValue ( "surfaceArea", numNucleolus, surfaceArea );
    individualNucleoliDataset.setValue ( "flatness", numNucleolus, flatness );
    individualNucleoliDataset.setValue ( "elongation", numNucleolus, elongation );
    individualNucleoliDataset.setValue ( "sphericity_vm", numNucleolus, sphericity );
    individualNucleoliDataset.setValue ( "sphericity_tm", numNucleolus, sphericity_tm );
  }

  totalNumNucleoli += regionAnalysisNucleoli.numRegions();

  //return chromocentersDataset;
}

