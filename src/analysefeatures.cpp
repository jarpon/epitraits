#include <iostream>
#include <otsuthresholding.h>
#include <regionanalysis.h>
#include <thresholding.h>
#include <voxelmatrix.h>
#include <watershedtransform.h>


#define TRACE
#include <trace.h>
{/*
VoxelMatrix <float> analyseFeatures( VoxelMatrix<float>& voxelMatrix, VoxelMatrix<float>& nucleusMask, VoxelMatrix<float>& gradientMatrix, const string& fileName, const string& outputDirectory )

  WatershedTransform <float> watershedTransform;
  watershedTransform.MaskVoxelMatrixProcessing <float>::setMask( nucleusMask );
  watershedTransform.apply( gradientMatrix );
  gradientMatrix.save( outputDirectory+"/nucleus/"+fileName+"-watershed.vm", true );

  VoxelMatrix <float> featureMatrix;
  RegionAnalysis <float> regionAnalysis;
  regionAnalysis.setRegionMatrix( gradientMatrix );
  regionAnalysis.run();
  regionAnalysis.mapRegionFeature( featureMatrix, REGION_FEATURE_AVERAGE_VALUE, voxelMatrix );
  //regionAnalysis.mapRegionFeature( featureMatrix, REGION_FEATURE_SIZE, voxelMatrix );
  featureMatrix.save( outputDirectory+"/nucleus/"+fileName+"-average_b.vm", true );


  OtsuThresholding <float> otsuThresholding;
  Vector <float> featureValues = regionAnalysis.allRegionValues( featureMatrix );
  Vector <unsigned int> histogram = featureValues.histogram( 0, 1, 256 );
  //EVAL( otsuThresholding.computeThreshold( histogram ) );

  Thresholding <float> thresholding;
  thresholding.setBackground( 0.0 );
  thresholding.setForeground( 255.0 );
  thresholding.setThreshold( otsuThresholding.computeThreshold(histogram) );
  thresholding.apply( featureMatrix );
  featureMatrix.save( outputDirectory+"/chromocenters/"+fileName+"-cc_b.vm", true );


  return featureMatrix;
}*/
