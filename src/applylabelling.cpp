#include <iostream>
#include <gaussiangradient.h>
#include <sobelgradient.h>
#include "voxelmatrix.h"
#include <watershedtransform.h>
#include "thresholding.h"
#include <volumehistogramexpansion.h>
#include <medianfilter.h>

#define TRACE
#include <trace.h>

VoxelMatrix <float> applyLabelling(const VoxelMatrix<float>& originalVoxelMatrix, const VoxelMatrix<float>& nucleusMask,
                                   const string& filename, const string& intermediateProcessesDir )
{
  VoxelMatrix<float> gradientMatrix = originalVoxelMatrix;
  VoxelMatrix<float> nucleusMaskCopy = nucleusMask;

//  MedianFilter<float> medianFilter;
//  medianFilter.setHalfSize( 2 );
//  medianFilter.setNumIterations( 2 );
//  medianFilter.apply( gradientMatrix );

  //int sizeZ = originalVoxelMatrix.getSize3();
  if ( nucleusMaskCopy.max().max().max() > 1 )
  {
      Thresholding<float> thresholding;
      thresholding.setBackground( 0.0 );
      thresholding.setForeground( 1.0 );
      thresholding.setThreshold( 1.0 );
      thresholding.apply( nucleusMaskCopy );
  }

  GaussianGradient<float> gaussianGradient;
  gaussianGradient.MaskVoxelMatrixProcessing<float>::setMask( nucleusMaskCopy );

  //for 8bits-images
  //gaussianGradient.setSigma( 1.3 );

  //for 16bits-images
  gaussianGradient.setSigma( 2 );

  //for (int k = 0; k < sizeZ; ++k) gaussianGradient.apply( gradientMatrix [k] );
  gaussianGradient.apply( gradientMatrix );
  //gradientMatrix.save( intermediateProcessesDir + filename + "-gradient.vm", true );
  //gradientMatrix.operator /=( 16 );

  VoxelMatrix <float> regionMatrix = gradientMatrix;
  WatershedTransform<float> watershedTransform;
  watershedTransform.MaskVoxelMatrixProcessing<float>::setMask( nucleusMaskCopy );
  watershedTransform.apply( regionMatrix );

 // regionMatrix.save( intermediateProcessesDir + filename + "-watershed.vm", true );

  return regionMatrix;
}
