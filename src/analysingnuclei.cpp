#include <dataset.h>
#include <regionanalysis3d.h>
#include <marchingcubes.h>
#include "thresholding.h"
#include <trimesh.h>
#include <cmath>
#include <regionanalysis2d.h>

#define TRACE
#include <trace.h>

void nucleusAnalysis(const VoxelMatrix<float>& originalVoxelMatrix, VoxelMatrix<float>& nucleusMask,
                     const string& filename, const string& parentDir, const int& numNucleus, DataSet& nucleiDataset)
{
  //const string analysisDir = parentDir + "/analysis/";
  const string shapesDir = parentDir + "/shapes/nuclei/";

//  PixelMatrix<float> zProj( nucleusMask.getSize1(), nucleusMask.getSize2() );
//  zProj.setZeros();

//  int i, j, k;

//  for (i = 0; i < nucleusMask.getSize1(); ++i)
//    for (j = 0; j < nucleusMask.getSize2(); ++j)
//      for (k = 0; k < nucleusMask.getSize3(); ++k)
//        if ( nucleusMask[k][i][j] > zProj(i,j) )
//          zProj(i,j) = nucleusMask[k][i][j];

//  PixelCalibration pixelCalibration;
//  pixelCalibration.setPixelHeight(originalVoxelMatrix.getVoxelCalibration().getVoxelHeight() );
//  pixelCalibration.setPixelWidth( originalVoxelMatrix.getVoxelCalibration().getVoxelWidth() );
//  pixelCalibration.setLengthUnit( originalVoxelMatrix.getVoxelCalibration().getLengthUnit() );
//  zProj.setPixelCalibration( pixelCalibration );

//  PixelMatrix<float> originalZProj( originalVoxelMatrix.getSize1(), originalVoxelMatrix.getSize2() );
//  originalZProj.setZeros();
//  for (i = 0; i < originalVoxelMatrix.getSize1(); ++i)
//    for (j = 0; j < originalVoxelMatrix.getSize2(); ++j)
//      for (k = 0; k < originalVoxelMatrix.getSize3(); ++k)
//        if ( originalVoxelMatrix[k][i][j] > originalZProj(i,j) )
//          originalZProj(i,j) = originalVoxelMatrix[k][i][j];

//  originalZProj.setPixelCalibration( pixelCalibration );

//  RegionAnalysis2D<float> regionAnalysis2D;
//  regionAnalysis2D.setLabelMatrix( zProj );
//  regionAnalysis2D.setValueMatrix( originalZProj );
//  regionAnalysis2D.run();

  Thresholding<float> thresholding;
  thresholding.setForeground( 1.0 );
  thresholding.setBackground( 0.0 );
  thresholding.setThreshold(0.5);
  thresholding.apply( nucleusMask );

  RegionAnalysis3D<float> regionAnalysis;
  regionAnalysis.setLabelMatrix( nucleusMask );
  regionAnalysis.setValueMatrix( originalVoxelMatrix );
  regionAnalysis.run();

  //generate 3Dmesh
  MarchingCubes<float> marchingCubes;
  TriMesh<float> triMesh;
  triMesh = marchingCubes.buildMesh( nucleusMask, 0.5, true );
  triMesh.scale( originalVoxelMatrix.getVoxelCalibration().getVoxelSize() );
  triMesh.save( shapesDir + filename , true );

  // get the name of the class
  string classif = parentDir;
  classif = classif.substr(classif.find_last_of("/\\")+1,classif.length());

  EVAL(regionAnalysis.computeRegionFeature(REGION_FEATURE_MAJOR_AXIS));
  nucleiDataset.setValue ( "name", numNucleus, filename );//filename
  nucleiDataset.setValue ( "class", numNucleus, classif );//classification: mutant, tissue, etc.
  nucleiDataset.setValue ( "nucleusVolume_vm", numNucleus, regionAnalysis.computeRegionFeature(REGION_FEATURE_VOLUME)[0] );//nucleus volume got from the voxelmatrix
  nucleiDataset.setValue ( "nucleusVolume_tm", numNucleus, abs(triMesh.volume()) );//nucleus volume got from the trimesh
  nucleiDataset.setValue ( "equivalentRadius_vm", numNucleus, regionAnalysis.computeRegionFeature(REGION_FEATURE_EQUIVALENT_RADIUS)[0]);
  nucleiDataset.setValue ( "equivalentRadius_tm", numNucleus, abs(triMesh.equivalentRadius()) );
  nucleiDataset.setValue ( "voxelSizeUnit", numNucleus, originalVoxelMatrix.getVoxelCalibration().getLengthUnit().symbol() + "^3" );//real voxel size unit
  nucleiDataset.setValue ( "major-axis", numNucleus, regionAnalysis.computeRegionFeature(REGION_FEATURE_MAJOR_AXIS)[0] );//
  nucleiDataset.setValue ( "intermediate-axis", numNucleus, regionAnalysis.computeRegionFeature(REGION_FEATURE_INTERMEDIATE_AXIS)[0] );//
  nucleiDataset.setValue ( "minor-axis", numNucleus, regionAnalysis.computeRegionFeature(REGION_FEATURE_MINOR_AXIS)[0] );//
  nucleiDataset.setValue ( "elongation", numNucleus, regionAnalysis.computeRegionFeature(REGION_FEATURE_ELONGATION)[0] );//elongation parameter
  nucleiDataset.setValue ( "flatness", numNucleus, regionAnalysis.computeRegionFeature(REGION_FEATURE_FLATNESS)[0] );//flatness parameter
  nucleiDataset.setValue ( "sphericity_vm", numNucleus, regionAnalysis.computeRegionFeature(REGION_FEATURE_COMPACTNESS)[0] );//sphericity parameter
  nucleiDataset.setValue ( "sphericity_tm", numNucleus, ( 36 * M_PI * pow(abs(triMesh.volume()) , 2) ) / pow( abs(triMesh.area() ) , 3) );
  nucleiDataset.setValue ( "surfaceArea_vm", numNucleus, regionAnalysis.computeRegionFeature(REGION_FEATURE_SURFACE_AREA)[0] );//surface area of the nucleus
  nucleiDataset.setValue ( "surfaceArea_tm", numNucleus, abs(triMesh.area() ) );//nucleus volume got from the trimesh
  //nucleiDataset.setValue ( "areaZprojection", numNucleus, regionAnalysis2D.computeRegionFeature(REGION_FEATURE_AREA)[0] );
  //nucleiDataset.setValue ( "circularity", numNucleus, regionAnalysis2D.computeRegionFeature(REGION_FEATURE_COMPACTNESS)[0] );
  //old:nucleiDataset.setValue ( "intensity", numNucleus, regionAnalysis.computeRegionFeature(REGION_FEATURE_INTEGRATED_INTENSITY)[0] * regionAnalysis.computeRegionFeature(REGION_FEATURE_VOLUME)[0]  );//sum of intensities of the nucleus volume

  //uncomment
  //  nucleiDataset.setValue ( "intensity", numNucleus, regionAnalysis.computeRegionFeature(REGION_FEATURE_SUM_VALUE)[0] );//sum of intensities of the nucleus volume
//  nucleiDataset.setValue ( "integratedDensity", numNucleus, regionAnalysis.computeRegionFeature(REGION_FEATURE_INTEGRATED_INTENSITY)[0] );//compute the integrated density taking into account the real size of the nucleus

  //return nucleiDataset;
}
