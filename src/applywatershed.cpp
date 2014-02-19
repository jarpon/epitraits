#include <iostream>
#include <voxelmatrix.h>
#include <watershedtransform.h>


void applyWatershed( VoxelMatrix<float>& nucleusMask, VoxelMatrix<float>& gradientMatrix, const string& fileName, const string& outputDirectory )
{
  WatershedTransform <float> watershedTransform;
  watershedTransform.MaskVoxelMatrixProcessing <float>::setMask( nucleusMask );
  watershedTransform.apply( gradientMatrix );
  gradientMatrix.save( outputDirectory+"/nucleus/"+fileName+"-watershed.vm", true );
}
