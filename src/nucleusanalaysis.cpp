#include <dataset.h>
//#include "regionanalysis2.h"
#include <regionanalysis.h>
#include <marchingcubes.h>
#include <thresholding.h>
#include <trimesh.h>

#define TRACE
#include <trace.h>

void nucleusAnalysis(const VoxelMatrix<float>& originalVoxelMatrix, VoxelMatrix<float>& nucleusMask,
                     const string& filename, const string& parentDir, const int& numNucleus, DataSet& nucleiDataset)
{
  const string analysisDir = parentDir + "/analysis/";

  RegionAnalysis<float> regionAnalysis;
  regionAnalysis.setRegionMatrix( nucleusMask );
  regionAnalysis.run();

  MarchingCubes<float> marchingCubes;
  TriMesh<float> triMesh;
  Thresholding<float> thresholding;
  thresholding.setForeground( 1.0 );
  thresholding.setBackground( 0.0 );

  //generate 3Dmesh
  thresholding.levelSetMask( nucleusMask, 1 );
  triMesh = marchingCubes.buildMesh( nucleusMask, 0.5, true );
  triMesh.scale( originalVoxelMatrix.getVoxelCalibration().getVoxelSize() );
  triMesh.save( analysisDir + filename , true );

  nucleiDataset.setValue ( "name", numNucleus, filename );//filename
  nucleiDataset.setValue ( "nucleusVolume_vm", numNucleus, regionAnalysis.computeRegionFeature(REGION_FEATURE_VOLUME,originalVoxelMatrix)[0] );//nucleus volume got from the voxelmatrix
  nucleiDataset.setValue ( "nucleusVolume_tm", numNucleus, triMesh.volume() );//nucleus volume got from the trimesh
  nucleiDataset.setValue ( "voxelSizeUnit", numNucleus, originalVoxelMatrix.getVoxelCalibration().getLengthUnit().symbol() + "^3" );//real voxel size unit
  nucleiDataset.setValue ( "flatness", numNucleus, regionAnalysis.computeRegionFeature(REGION_FEATURE_FLATNESS,originalVoxelMatrix)[0] );//flatness parameter
  nucleiDataset.setValue ( "elongation", numNucleus, regionAnalysis.computeRegionFeature(REGION_FEATURE_ELONGATION,originalVoxelMatrix)[0] );//elongation parameter
  nucleiDataset.setValue ( "sphericity", numNucleus, regionAnalysis.computeRegionFeature(REGION_FEATURE_SPHERICITY,originalVoxelMatrix)[0] );//sphericity parameter
  nucleiDataset.setValue ( "surfaceArea", numNucleus, regionAnalysis.computeRegionFeature(REGION_FEATURE_SURFACE_AREA,originalVoxelMatrix)[0] );//surface area of the nucleus
  nucleiDataset.setValue ( "intensity", numNucleus, regionAnalysis.computeRegionFeature(REGION_FEATURE_INTENSITY,originalVoxelMatrix)[0] );//sum of intensities of the nucleus volume
  nucleiDataset.setValue ( "integratedDensity", numNucleus, regionAnalysis.computeRegionFeature(REGION_FEATURE_INTEGRATED_DENSITY,originalVoxelMatrix)[0] );//compute the integrated density taking into account the real size of the nucleus

  //return nucleiDataset;
}
