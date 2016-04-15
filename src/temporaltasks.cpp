#include <voxelmatrix.h>

#define TRACE
#include <trace.h>
using namespace std;


void addCalibration(const VoxelMatrix<float>& notCalibratedVM, const VoxelMatrix<float>& alreadyCalibratedVM, const string& outputDir, const string& outputFilename )
{
  VoxelMatrix<float> copyVM( notCalibratedVM );
  VoxelCalibration voxelCalibration;
  voxelCalibration.setVoxelHeight( alreadyCalibratedVM.getVoxelCalibration().getVoxelHeight() );
  voxelCalibration.setVoxelWidth( alreadyCalibratedVM.getVoxelCalibration().getVoxelWidth() );
  voxelCalibration.setVoxelDepth( alreadyCalibratedVM.getVoxelCalibration().getVoxelDepth() );
  voxelCalibration.setLengthUnit( alreadyCalibratedVM.getVoxelCalibration().getLengthUnit() );

  copyVM.setVoxelCalibration( voxelCalibration );
  copyVM.save( outputDir + outputFilename + ".vm", false );
}
