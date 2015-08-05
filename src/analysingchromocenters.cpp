#include <componentlabelling.h>
#include <dataset.h>
//#include "regionanalysis2.h"
//#include "regionanalysis.h"
#include <regionanalysis3d.h>
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

  RegionAnalysis3D<float> regionAnalysis;
  regionAnalysis.setLabelMatrix( nucleusMask );
  regionAnalysis.setValueMatrix( originalVoxelMatrix );
  regionAnalysis.run();

  float nucleusVolume = regionAnalysis.computeRegionFeature( REGION_FEATURE_VOLUME )[0];
  float nucleusIntensity = regionAnalysis.computeRegionFeature( REGION_FEATURE_INTEGRATED_INTENSITY )[0] * nucleusVolume;

//  ComponentLabelling<float> componentLabelling;
//  componentLabelling.apply( ccsMask );

  RegionAnalysis3D<float> regionAnalysisCCs;
  regionAnalysisCCs.setLabelMatrix( ccsMask );
  regionAnalysis.setValueMatrix( originalVoxelMatrix );
  regionAnalysisCCs.run();

  //TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + "_nucleus.tm" );
  TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + ".tm" );

  Vector<float> centroid(3);
  Vector<float> vertexTriMesh(3);

  // get the name of the class
  string classif = parentDir;
  classif = classif.substr(classif.find_last_of("/\\")+1,classif.length());

  Vertices<float> centroids = regionAnalysisCCs.regionCentroids();
  Vector<float> ccsVolume = regionAnalysisCCs.computeRegionFeature( REGION_FEATURE_VOLUME );
  Vector<float> ccsEqRadius = regionAnalysisCCs.computeRegionFeature( REGION_FEATURE_EQUIVALENT_RADIUS );
  Vector<float> ccsSurfaceArea = regionAnalysisCCs.computeRegionFeature( REGION_FEATURE_SURFACE_AREA );
  //Vector<float> ccsRelativeVolume = ccsVolume.operator *( ccsVolume.sum() );
  Vector<float> ccsRelativeVolume = ccsVolume.operator *( nucleusVolume ); //relative volume of cc within the nucleus
  Vector<float> ccsFlatness = regionAnalysisCCs.computeRegionFeature( REGION_FEATURE_FLATNESS );
  Vector<float> ccsElongation = regionAnalysisCCs.computeRegionFeature( REGION_FEATURE_ELONGATION );
  Vector<float> ccsSphericity = regionAnalysisCCs.computeRegionFeature( REGION_FEATURE_COMPACTNESS );
  Vector<float> ccsIntegratedDensity = regionAnalysisCCs.computeRegionFeature( REGION_FEATURE_INTEGRATED_INTENSITY );

  nucleiDataset.setValue ( "name", numNucleus, filename );//filename
  nucleiDataset.setValue ( "class", numNucleus, classif );//classification: mutant, tissue, etc.
  nucleiDataset.setValue ( "ccsNumber", numNucleus, regionAnalysisCCs.numRegions() );//number of ccs obtained in the nucleus
  nucleiDataset.setValue ( "ccsVolume", numNucleus, ccsVolume.sum() );//total ccs volume
  nucleiDataset.setValue ( "volRHF", numNucleus, ccsVolume.sum() / nucleusVolume  );//relative volume of ccs regarding the complete nucleus
  nucleiDataset.setValue ( "ccsIntensity", numNucleus, ccsIntegratedDensity.sum()*ccsVolume.sum() );//total absolute intensity of ccs
  nucleiDataset.setValue ( "ccsIntegratedDensity", numNucleus, ccsIntegratedDensity.sum() );//integrated density of ccs taking into account real volume
  nucleiDataset.setValue ( "intRHF", numNucleus, ( ccsIntegratedDensity.sum()*ccsVolume.sum() ) / nucleusIntensity );//RHF: rate of heterochromatin (int.Density of ccs/ int.Density of nuclei)

/**///chromocenters individual information
  for (int numCC = 0; numCC < regionAnalysisCCs.numRegions(); numCC++ )
  {
    chromocentersDataset.setValue ( "name", numCC+totalNumCCs, filename );
    chromocentersDataset.setValue ( "class", numCC+totalNumCCs, classif );//classification: mutant, tissue, etc.
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

    chromocentersDataset.setValue ( "equivalentRadius_vm", numCC+totalNumCCs, ccsEqRadius[numCC] );
    chromocentersDataset.setValue ( "equivalentRadius_tm", numCC+totalNumCCs, eqRadius_tm );
    chromocentersDataset.setValue ( "ccVolume_vm", numCC+totalNumCCs, ccsVolume[numCC] );
    chromocentersDataset.setValue ( "ccVolume_tm", numCC+totalNumCCs, ccVolume_tm );
    chromocentersDataset.setValue ( "distanceToTheBorder", numCC+totalNumCCs, distanceToBorder );
    chromocentersDataset.setValue ( "ccRelativeVolume", numCC+totalNumCCs, ccsRelativeVolume[numCC] );
    chromocentersDataset.setValue ( "flatness", numCC+totalNumCCs, ccsFlatness[numCC] );
    chromocentersDataset.setValue ( "elongation", numCC+totalNumCCs, ccsElongation[numCC] );
    chromocentersDataset.setValue ( "sphericity", numCC+totalNumCCs, ccsSphericity[numCC] );
    chromocentersDataset.setValue ( "surfaceArea", numCC+totalNumCCs, ccsSurfaceArea[numCC] );
    chromocentersDataset.setValue ( "ccsIntensity", numCC+totalNumCCs, ccsIntegratedDensity[numCC] * ccsVolume[numCC] );
    chromocentersDataset.setValue ( "ccsIntegratedDensity", numCC+totalNumCCs, ccsIntegratedDensity[numCC] );
    //this last relativeCCsIntensity is the intensity of each cc divided by the total cc's intensity
    chromocentersDataset.setValue ( "relativeCCsIntensity", numCC+totalNumCCs, ccsIntegratedDensity[numCC] * ccsVolume[numCC] / ( ccsIntegratedDensity.sum() * ccsVolume.sum() ) );

    individualChromocentersDataset.setValue ( "id", numCC, filename );
    individualChromocentersDataset.setValue ( "idCC", numCC, numCC+1 );
    individualChromocentersDataset.setValue ( "centroidCoordX", numCC, centroid[X] );
    individualChromocentersDataset.setValue ( "centroidCoordY", numCC, centroid[Y] );
    individualChromocentersDataset.setValue ( "centroidCoordZ", numCC, centroid[Z] );
    individualChromocentersDataset.setValue ( "equivalentRadius_vm", numCC, ccsEqRadius[numCC] );
    individualChromocentersDataset.setValue ( "equivalentRadius_tm", numCC, eqRadius_tm );
    individualChromocentersDataset.setValue ( "distanceToTheBorder", numCC, distanceToBorder );
    individualChromocentersDataset.setValue ( "ccVolume_vm", numCC, ccsVolume[numCC] );
    individualChromocentersDataset.setValue ( "ccVolume_tm", numCC, ccVolume_tm );
    individualChromocentersDataset.setValue ( "ccRelativeVolume", numCC, ccsRelativeVolume[numCC] );
    individualChromocentersDataset.setValue ( "flatness", numCC, ccsFlatness[numCC] );
    individualChromocentersDataset.setValue ( "elongation", numCC, ccsElongation[numCC] );
    individualChromocentersDataset.setValue ( "sphericity", numCC, ccsSphericity[numCC] );
    individualChromocentersDataset.setValue ( "surfaceArea", numCC, ccsSurfaceArea[numCC] );
    individualChromocentersDataset.setValue ( "ccsIntensity", numCC, ccsIntegratedDensity[numCC] * ccsVolume[numCC] );
    individualChromocentersDataset.setValue ( "ccsIntegratedDensity", numCC, ccsIntegratedDensity[numCC] );
    individualChromocentersDataset.setValue ( "relativeCCsIntensity", numCC, ccsIntegratedDensity[numCC] * ccsVolume[numCC] / ( ccsIntegratedDensity.sum() * ccsVolume.sum() ) );

  }

  totalNumCCs += regionAnalysisCCs.numRegions();

  //return chromocentersDataset;
}
