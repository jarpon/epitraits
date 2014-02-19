#include <iostream>
#include <gaussiangradient.h>
#include <sobelgradient.h>
#include <voxelmatrix.h>
#include <watershedtransform.h>
#include <thresholding.h>
#include <volumehistogramexpansion.h>
#include <minfilter.h>
#include <maxfilter.h>

#define TRACE
#include <trace.h>

VoxelMatrix <float> applyLabelling(const VoxelMatrix<float>& originalVoxelMatrix, const VoxelMatrix<float>& nucleusMask,
                                   const string& filename, const string& intermediateProcessesDir )
{
  VoxelMatrix<float> gradientMatrix = originalVoxelMatrix;

  //int sizeZ = originalVoxelMatrix.getSize3();

  GaussianGradient<float> gaussianGradient;
  gaussianGradient.MaskVoxelMatrixProcessing<float>::setMask( nucleusMask );
  gaussianGradient.setSigma( 1.4 );
  //for (int k = 0; k < sizeZ; ++k) gaussianGradient.apply( gradientMatrix [k] );
  gaussianGradient.apply( gradientMatrix );
  gradientMatrix.save( intermediateProcessesDir + filename + "-gradient.vm", true );

  VoxelMatrix <float> regionMatrix = gradientMatrix;
  WatershedTransform<float> watershedTransform;
  watershedTransform.MaskVoxelMatrixProcessing<float>::setMask( nucleusMask );
  watershedTransform.apply( regionMatrix );

  regionMatrix.save( intermediateProcessesDir + filename + "-watershed.vm", true );

  return regionMatrix;
}
