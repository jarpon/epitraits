#include <dataset.h>
//#include "regionanalysis2.h"
//#include "regionanalysis.h"
#include <regionanalysis3d.h>
#include <marchingcubes.h>
#include "thresholding.h"
#include <trimesh.h>
#include <cmath>

#define TRACE
#include <trace.h>

void nucleusAnalysis(const VoxelMatrix<float>& originalVoxelMatrix, VoxelMatrix<float>& nucleusMask,
                     const string& filename, const string& parentDir, const int& numNucleus, DataSet& nucleiDataset)
{
  //const string analysisDir = parentDir + "/analysis/";
  const string shapesDir = parentDir + "/shapes/";

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

  nucleiDataset.setValue ( "name", numNucleus, filename );//filename
  nucleiDataset.setValue ( "class", numNucleus, classif );//classification: mutant, tissue, etc.
  nucleiDataset.setValue ( "nucleusVolume_vm", numNucleus, regionAnalysis.computeRegionFeature(REGION_FEATURE_VOLUME,originalVoxelMatrix)[0] );//nucleus volume got from the voxelmatrix
  nucleiDataset.setValue ( "nucleusVolume_tm", numNucleus, abs(triMesh.volume()) );//nucleus volume got from the trimesh
  nucleiDataset.setValue ( "equivalentRadius_vm", numNucleus, regionAnalysis.computeRegionFeature(REGION_FEATURE_EQUIVALENT_RADIUS,originalVoxelMatrix)[0]);
  nucleiDataset.setValue ( "equivalentRadius_tm", numNucleus, abs(triMesh.equivalentRadius()) );
  nucleiDataset.setValue ( "voxelSizeUnit", numNucleus, originalVoxelMatrix.getVoxelCalibration().getLengthUnit().symbol() + "^3" );//real voxel size unit
  nucleiDataset.setValue ( "flatness", numNucleus, regionAnalysis.computeRegionFeature(REGION_FEATURE_FLATNESS,originalVoxelMatrix)[0] );//flatness parameter
  nucleiDataset.setValue ( "elongation", numNucleus, regionAnalysis.computeRegionFeature(REGION_FEATURE_ELONGATION,originalVoxelMatrix)[0] );//elongation parameter
  nucleiDataset.setValue ( "sphericity_vm", numNucleus, regionAnalysis.computeRegionFeature(REGION_FEATURE_SPHERICITY,originalVoxelMatrix)[0] );//sphericity parameter
  nucleiDataset.setValue ( "sphericity_tm", numNucleus, ( 36 * M_PI * pow(abs(triMesh.volume()) , 2) ) / pow( abs(triMesh.area() ) , 3) );
  nucleiDataset.setValue ( "surfaceArea_vm", numNucleus, regionAnalysis.computeRegionFeature(REGION_FEATURE_SURFACE_AREA,originalVoxelMatrix)[0] );//surface area of the nucleus
  nucleiDataset.setValue ( "surfaceArea_tm", numNucleus, abs(triMesh.area() ) );//nucleus volume got from the trimesh
  nucleiDataset.setValue ( "intensity", numNucleus, regionAnalysis.computeRegionFeature(REGION_FEATURE_INTENSITY,originalVoxelMatrix)[0] );//sum of intensities of the nucleus volume
  nucleiDataset.setValue ( "integratedDensity", numNucleus, regionAnalysis.computeRegionFeature(REGION_FEATURE_INTEGRATED_DENSITY,originalVoxelMatrix)[0] );//compute the integrated density taking into account the real size of the nucleus

  //return nucleiDataset;
}
