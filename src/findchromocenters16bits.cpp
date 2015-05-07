#include <iostream>
#include <componentlabelling.h>
#include <otsuthresholding.h>
//#include "regionanalysis2.h"
#include <regionanalysis.h>
#include <thresholding.h>
#include <voxelmatrix.h>
#include <watershedtransform.h>
#include <voxelmatrixdilatation.h>
#include <voxelmatrixerosion.h>
#include <holesfilling.h>
#include <gaussiangradient.h>

#include <cmath>


#define TRACE
#include <trace.h>

VoxelMatrix <float> findCCs16bits(const VoxelMatrix<float>& originalVoxelMatrix, VoxelMatrix<float>& nucleusMask,
                            const string& filename, const string& intermediateProcessesDir )
{
  int sizeZ = originalVoxelMatrix.getSize3();

  VoxelMatrix<float> gradientMatrix = originalVoxelMatrix;
  GaussianGradient<float> gaussianGradient;
  gaussianGradient.MaskVoxelMatrixProcessing<float>::setMask( nucleusMask );
  gaussianGradient.setSigma( 5 );
  for (int k = 0; k < sizeZ; ++k) gaussianGradient.apply( gradientMatrix [k] );
  //gaussianGradient.apply( gradientMatrix );
  gradientMatrix.save( intermediateProcessesDir + filename + "-gradient.vm", true );
  //gradientMatrix.operator /=( 16 );

  VoxelMatrix <float> regionMatrix = gradientMatrix;
  WatershedTransform<float> watershedTransform;
  watershedTransform.MaskVoxelMatrixProcessing<float>::setMask( nucleusMask );
  watershedTransform.apply( regionMatrix );

  regionMatrix.save( intermediateProcessesDir + filename + "-watershed.vm", true );
//  VoxelMatrix<float> reducedVoxelMatrix = originalVoxelMatrix;
//  Thresholding<float> firstThresholding;
//  firstThresholding.setBackground( 0.0 );
//  firstThresholding.setThreshold( 1000.0 );
//  firstThresholding.applyAlternative( reducedVoxelMatrix );


  //const int sizeZ = originalVoxelMatrix.getSize3();

  VoxelMatrix<float> rangeMask, copyVoxelMatrix = originalVoxelMatrix;

  RegionAnalysis<float> firstRegionAnalysis;
  firstRegionAnalysis.setRegionMatrix( regionMatrix );
  firstRegionAnalysis.run();

//  firstRegionAnalysis.mapRegionFeature( copyVoxelMatrix, REGION_FEATURE_MAXIMUM_VALUE, originalVoxelMatrix );
//  copyVoxelMatrix.save( intermediateProcessesDir + filename + "-average.vm", true );
  firstRegionAnalysis.mapRegionFeature( copyVoxelMatrix, REGION_FEATURE_CONTRAST, originalVoxelMatrix );
  copyVoxelMatrix.save( intermediateProcessesDir + filename + "-1stContrast.vm", true );

  VoxelMatrix<float> newNucleusMask = copyVoxelMatrix;
  regionMatrix = originalVoxelMatrix;

  Thresholding<float> firstThresholding;
  firstThresholding.setBackground( 0.0 );
  firstThresholding.setThreshold( newNucleusMask.max().max().max()/2 );
  firstThresholding.apply( newNucleusMask );
  newNucleusMask.save( intermediateProcessesDir + filename + "-mask.vm", true );

  VoxelMatrix<float> structElement;
  structElement.setSize(11,11,11);
  structElement.setOnes();

  VoxelMatrixDilatation<float> voxelDilatation;
  voxelDilatation.setStructElt( structElement );
  voxelDilatation.apply( newNucleusMask );



  GaussianGradient<float> gaussianGradient2;
  gaussianGradient2.MaskVoxelMatrixProcessing<float>::setMask( newNucleusMask );
  gaussianGradient2.setSigma( 5 );
  for (int k = 0; k < sizeZ; ++k) gaussianGradient.apply( regionMatrix [k] );

  WatershedTransform<float> watershedTransform2;
  watershedTransform2.MaskVoxelMatrixProcessing<float>::setMask( newNucleusMask );
  watershedTransform2.apply( regionMatrix );

  RegionAnalysis<float> regionAnalysis;
  regionAnalysis.setRegionMatrix( regionMatrix );
  regionAnalysis.run();

  regionAnalysis.mapRegionFeature( rangeMask, REGION_FEATURE_AVERAGE_VALUE, copyVoxelMatrix );
  rangeMask.save( intermediateProcessesDir + filename + "-average.vm", true );
//  regionAnalysis.mapRegionFeature( rangeMask, REGION_FEATURE_CONTRAST, copyVoxelMatrix );
//  rangeMask.save( intermediateProcessesDir + filename + "-contrast.vm", true );
//  regionAnalysis.mapRegionFeature( rangeMask, REGION_FEATURE_MAXIMUM_VALUE, copyVoxelMatrix );
//  rangeMask.save( intermediateProcessesDir + filename + "-max.vm", true );

    VoxelMatrix<float> ccsMask = rangeMask;

    OtsuThresholding<float> otsuThresholding;

    //Vector <float> featureValues = regionAnalysis.computeRegionFeature( REGION_FEATURE_CONTRAST, copyVoxelMatrix );
  //  Vector <float> featureValues = regionAnalysis.computeRegionFeature( REGION_FEATURE_CONTRAST, rangeMask );
//    Vector <float> featureValues = regionAnalysis.computeRegionFeature( REGION_FEATURE_AVERAGE_VALUE, rangeMask );
//    Vector <unsigned int> histogram = featureValues.histogram( 0, 1, histogram.getSize()-1 );

    //EVAL(featureValues);
//    EVAL(histogram);
  Thresholding<float> thresholding;
  thresholding.setBackground( 0.0 );



  thresholding.setThreshold( otsuThresholding.computeThreshold(rangeMask) );
  thresholding.applyAlternative( ccsMask );



//  HolesFillingf holesFilling;
//  int sizeZ = originalVoxelMatrix.getSize3();

//  //to obtain a better filling, it's applied to each 2D slice instead of the complete 3D stack
//  for (int k = 0; k < sizeZ; ++k)  holesFilling.apply( ccsMask[k] );
//  holesFilling.apply( ccsMask );

  //for (int k = 0; k < sizeZ; ++k)  holesFilling.apply( nucleusMask[k] );
  //labeling the image
  ComponentLabelling<float> componentLabelling;
  componentLabelling.apply( ccsMask );


//  VoxelMatrix<float> structElement;
//  structElement.setSize(3,3,3);
//  structElement.setOnes();

//  VoxelMatrixDilatation<float> voxelDilatation;
//  voxelDilatation.setStructElt( structElement );
//  voxelDilatation.apply( ccsMask );

// //  VoxelMatrix<float> structElement2;
// //  structElement2.setSize(1,1,1);
// //  structElement2.setOnes();

//  VoxelMatrixErosion<float> voxelErosion;
//  voxelErosion.setStructElt( structElement );
//  voxelErosion.apply( ccsMask );

  ccsMask.setVoxelCalibration( originalVoxelMatrix.getVoxelCalibration() );

  //ccsMask.save( chromocentersDir + filename + ".vm", true );

  return ccsMask;
}

