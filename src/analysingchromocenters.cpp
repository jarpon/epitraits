#include <componentlabelling.h>
#include <dataset.h>
#include "regionanalysis2d.h"
#include "regionanalysis3d.h"
//#include <regionanalysis3d.h>
//#include <regionanalysis2d.h>
#include <cmath>
#include <marchingcubes.h>
#include <thresholding.h>
#include <trimesh.h>

#define TRACE
#include <trace.h>

//PixelMatrix<float> getMaximumIntensity(const VoxelMatrix<float>& ccsMask)
//{
//  const int size1 = ccsMask.getSize1();
//  const int size2 = ccsMask.getSize2();
//  const int size3 = ccsMask.getSize3();

//  PixelMatrix<float> maxIntensity( size1, size2 );
//  maxIntensity.setZeros();

//  for ( int k = 0; k < size3; ++k )
//  {
//    for (int i = 0; i < size1; ++i )
//    {
//      for ( int j = 0; j < size2; ++j )
//      {
//        if ( maxIntensity(i,j) < ccsMask(i,j,k) )
//          maxIntensity(i,j) = ccsMask(i,j,k);
//      }
//    }
//  }

//  EVAL(ccsMask.getVoxelCalibration().getVoxelHeight());
//  EVAL(ccsMask.getVoxelCalibration().getVoxelWidth());
//  PixelCalibration  pixelCalibration;
//  pixelCalibration.setPixelHeight( ccsMask.getVoxelCalibration().getVoxelHeight() );
//  pixelCalibration.setPixelWidth( ccsMask.getVoxelCalibration().getVoxelWidth() );
//  pixelCalibration.setLengthUnit( ccsMask.getVoxelCalibration().getLengthUnit() );

//  maxIntensity.setPixelCalibration( pixelCalibration );
//  EVAL(maxIntensity.getPixelCalibration().getPixelHeight());
//  EVAL(maxIntensity.getPixelCalibration().getPixelWidth());
//  return maxIntensity;
//}


void chromocentersAnalysis(VoxelMatrix<float>& ccsMask, const string& filename, const string& parentDir,
                           const int& numNucleus, int& totalNumCCs,
                           DataSet& nucleiDataset, DataSet& chromocentersDataset, DataSet& individualChromocentersDataset)
{
  string originalName = filename.substr( 0,filename.find_last_of("-")  );

  //VoxelMatrix<float> originalVoxelMatrix( parentDir + "/segmented_nuclei/" + filename + ".vm" );
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

//  ComponentLabelling<float> componentLabelling;
//  componentLabelling.apply( ccsMask );

  RegionAnalysis3D<float> regionAnalysisCCs;
  regionAnalysisCCs.setLabelMatrix( ccsMask );
  regionAnalysisCCs.setValueMatrix( originalVoxelMatrix );
  regionAnalysisCCs.run();
  EVAL (regionAnalysisCCs.numRegions() );
  EVAL(regionAnalysisCCs.getRegion(1).getLabel());


  int num = regionAnalysisCCs.condenseRegionLabels();
  regionAnalysisCCs.run();

  //TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + "_nucleus.tm" );
  TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + ".tm" );
  float nucleusVolume = fabs(nucleusTriMesh.volume());

  //float nucleusVolume = regionAnalysis.computeRegionFeature( REGION_FEATURE_VOLUME )[0];
  float nucleusIntensity = regionAnalysis.computeRegionFeature( REGION_FEATURE_INTEGRATED_INTENSITY )[0] * nucleusVolume;

  EVAL(nucleusVolume);

  Vector<float> centroid(3);
  Vector<float> vertexTriMesh(3);

  // get the name of the class
  string classif = parentDir;
  classif = classif.substr(classif.find_last_of("/\\")+1,classif.length());

  Vertices<float> centroids = regionAnalysisCCs.regionCentroids();
  //Vector<float> ccsVolume = regionAnalysisCCs.computeRegionFeature( REGION_FEATURE_VOLUME );
  Vector<float> ccsVolume = regionAnalysisCCs.computeRegionFeature( REGION_FEATURE_VOLUME )/pow(3.,1./2.);//SPF correction
  //Vector<float> ccsEqRadius = regionAnalysisCCs.computeRegionFeature( REGION_FEATURE_EQUIVALENT_RADIUS )/pow(3.,1./6.);//SPF correction -> sqrt(3) volume correction -> eq radius
  Vector<float> ccsEqRadius;
  ccsEqRadius.setSize( regionAnalysisCCs.numRegions() );
  RegionAnalysis2D<float> regionAnalysis2D;
//  Vector<float> regionValues;
//  regionValues = regionAnalysisCCs.RegionAnalysisBase<float>::allRegionValues( );

  for ( int i = 0; i < regionAnalysisCCs.numRegions(); ++i )
  {
    PixelMatrix<float> ccs2DMask;
    ccs2DMask = regionAnalysisCCs.getLabel2DProjection( regionAnalysisCCs.getRegions()[i].getLabel() );
//    EVAL(regionAnalysisCCs.getRegions()[i].getLabel());
//    ccs2DMask.saveAsImage( parentDir + "/analysis/" + originalName + ".tif", true );

    regionAnalysis2D.setLabelMatrix( ccs2DMask );
    //regionAnalysis2D.setValueMatrix( getMaximumIntensity( originalVoxelMatrix ) );
    regionAnalysis2D.run();
    ccsEqRadius[i] = regionAnalysis2D.computeRegionFeature( REGION_FEATURE_EQUIVALENT_RADIUS )[0];// max area correction -> eq radius
//    ccsEqRadius[i] = regionAnalysis2D.computeRegionFeature( REGION_FEATURE_EQUIVALENT_RADIUS )[0]/sqrt(3.);//SPF correction -> max area correction -> eq radius
    EVAL(ccsEqRadius[i]);
  }

  Vector<float> ccsSurfaceArea = regionAnalysisCCs.computeRegionFeature( REGION_FEATURE_SURFACE_AREA )/pow(3.,1./4.);//SPF correction
  //Vector<float> ccsRelativeVolume = ccsVolume.operator *( ccsVolume.sum() );
  //Vector<float> ccsRelativeVolume = ccsVolumeR); //relative volume of cc within the nucleus
  Vector<float> ccsFlatness = regionAnalysisCCs.computeRegionFeature( REGION_FEATURE_FLATNESS );
  Vector<float> ccsElongation = regionAnalysisCCs.computeRegionFeature( REGION_FEATURE_ELONGATION );
  Vector<float> ccsSphericity = regionAnalysisCCs.computeRegionFeature( REGION_FEATURE_COMPACTNESS );
  Vector<float> ccsIntensity = regionAnalysisCCs.computeRegionFeature( REGION_FEATURE_SUM_VALUE );
  Vector<float> ccsIntegratedDensity = regionAnalysisCCs.computeRegionFeature( REGION_FEATURE_INTEGRATED_INTENSITY );
  EVAL(ccsIntegratedDensity.sum()*ccsVolume.sum())
      EVAL(ccsIntegratedDensity.sum())

  nucleiDataset.setValue ( "name", numNucleus, filename );//filename
  nucleiDataset.setValue ( "class", numNucleus, classif );//classification: mutant, tissue, etc.
  nucleiDataset.setValue ( "ccsNumber", numNucleus, regionAnalysisCCs.numRegions() );//number of ccs obtained in the nucleus
  nucleiDataset.setValue ( "avgeCCVolume", numNucleus, ccsVolume.mean() );
  nucleiDataset.setValue ( "ccsTotalVolume", numNucleus, ccsVolume.sum() );//total ccs volume
  nucleiDataset.setValue ( "volRHF", numNucleus, ccsVolume.sum()/nucleusVolume );//relative volume of ccs regarding the complete nucleus
  nucleiDataset.setValue ( "ccsIntensity", numNucleus, ccsIntensity.sum() );//total absolute intensity of ccs
  nucleiDataset.setValue ( "ccsIntegratedDensity", numNucleus, ccsIntegratedDensity.sum() );//integrated density of ccs taking into account real volume
  nucleiDataset.setValue ( "intRHF", numNucleus, ( ccsIntensity.sum() ) / nucleusIntensity );//RHF: rate of heterochromatin (int.Density of ccs/ int.Density of nuclei)


  int numCompartments = regionAnalysisCCs.numRegions();

  float tempDistance;
  Vector<float> temp1, temp2;
  Vector<float> min, max;
  min.setSize( numCompartments );
  max.setSize( numCompartments );

  EVAL ( num );

/**///chromocenters individual information
  for (int numCC = 0; numCC < numCompartments; numCC++ )
  {
    chromocentersDataset.setValue ( "name", numCC+totalNumCCs, filename );
    chromocentersDataset.setValue ( "class", numCC+totalNumCCs, classif );//classification: mutant, tissue, etc.
    chromocentersDataset.setValue ( "idCC", numCC+totalNumCCs, numCC+1 );

    centroid = centroids[numCC];

    min[numCC] = FLT_MAX;
    max[numCC] = FLT_MIN;

    for (int j = 0; j < numCompartments; ++j)
    {
      if ( numCC != j )
      {
        {
          temp2 = centroids[j];
          tempDistance = centroid.distance(temp2);
        }

        if ( tempDistance < min[numCC] )
          min[numCC] = tempDistance;

        if ( tempDistance > max[numCC] )
          max[numCC] = tempDistance;
      }
    }

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
    float eqRadius_tm = triMesh.equivalentRadius()/pow(3.,1./6.);//SPF correction -> max area correction -> eq radius

    chromocentersDataset.setValue ( "equivalentRadius_vm", numCC+totalNumCCs, ccsEqRadius[numCC] );
    chromocentersDataset.setValue ( "equivalentRadius_tm", numCC+totalNumCCs, eqRadius_tm );
    chromocentersDataset.setValue ( "ccVolume_vm", numCC+totalNumCCs, ccsVolume[numCC] );
    chromocentersDataset.setValue ( "ccVolume_tm", numCC+totalNumCCs, ccVolume_tm );
    chromocentersDataset.setValue ( "distanceToTheBorder", numCC+totalNumCCs, distanceToBorder );
    chromocentersDataset.setValue ( "ccRelativeVolume", numCC+totalNumCCs, ccVolume_tm/nucleusVolume );
    chromocentersDataset.setValue ( "flatness", numCC+totalNumCCs, ccsFlatness[numCC] );
    chromocentersDataset.setValue ( "elongation", numCC+totalNumCCs, ccsElongation[numCC] );
    chromocentersDataset.setValue ( "sphericity", numCC+totalNumCCs, ccsSphericity[numCC] );
    chromocentersDataset.setValue ( "surfaceArea", numCC+totalNumCCs, ccsSurfaceArea[numCC] );
    chromocentersDataset.setValue ( "ccsIntensity", numCC+totalNumCCs, ccsIntensity[numCC] );
    chromocentersDataset.setValue ( "ccsIntegratedDensity", numCC+totalNumCCs, ccsIntegratedDensity[numCC] );
    //this last relativeCCsIntensity is the intensity of each cc divided by the total cc's intensity
    chromocentersDataset.setValue ( "relativeCCsIntensity", numCC+totalNumCCs, ccsIntensity[numCC] / ccsIntegratedDensity.sum() );
    chromocentersDataset.setValue ( "minDistanceToCC", numCC+totalNumCCs, min[numCC] );
    chromocentersDataset.setValue ( "maxDistanceToCC", numCC+totalNumCCs, max[numCC] );

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
    individualChromocentersDataset.setValue ( "ccRelativeVolume", numCC, ccsVolume[numCC]/ccsVolume.sum() );
    individualChromocentersDataset.setValue ( "flatness", numCC, ccsFlatness[numCC] );
    individualChromocentersDataset.setValue ( "elongation", numCC, ccsElongation[numCC] );
    individualChromocentersDataset.setValue ( "sphericity", numCC, ccsSphericity[numCC] );
    individualChromocentersDataset.setValue ( "surfaceArea", numCC, ccsSurfaceArea[numCC] );
    individualChromocentersDataset.setValue ( "ccsIntensity", numCC, ccsIntensity[numCC] );
    individualChromocentersDataset.setValue ( "ccsIntegratedDensity", numCC, ccsIntegratedDensity[numCC] );
    individualChromocentersDataset.setValue ( "relativeCCsIntensity", numCC, ccsIntegratedDensity[numCC] * ccsVolume[numCC] / ( ccsIntegratedDensity.sum() * ccsVolume.sum() ) );
    individualChromocentersDataset.setValue ( "minDistanceToCC", numCC, min[numCC] );
    individualChromocentersDataset.setValue ( "maxDistanceToCC", numCC, max[numCC] );

  }

  nucleiDataset.setValue ( "minDistanceToCC", numNucleus, chromocentersDataset.getValues<float>("minDistanceToCC").min() );
  nucleiDataset.setValue ( "maxDistanceToCC", numNucleus, chromocentersDataset.getValues<float>("minDistanceToCC").max() );

  totalNumCCs += regionAnalysisCCs.numRegions();

  //return chromocentersDataset;
}
