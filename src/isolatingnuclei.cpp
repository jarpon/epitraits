#include <componentlabelling.h>
#include <holesfilling.h>
#include <otsuthresholding.h>
#include <thresholding.h>
#include <voxelmatrix.h>
#include <volumehistogramexpansion.h>
//#include <maxfilter.h>
//#include <minfilter.h>
#include <voxelmatrixerosion.h>
#include <voxelmatrixdilatation.h>
#include <regionanalysis.h>
#include <watershedtransform.h>
#include <gaussiangradient.h>


#define TRACE
#include <trace.h>
using namespace std;

VoxelMatrix<float> isolateNuclei(const VoxelMatrix<float>& originalVoxelMatrix)
{
  ENTER("isolateNuclei");

  VoxelMatrix<float> nucleiMask = originalVoxelMatrix;

//  VolumeHistogramExpansion<float> normalizerHistogram;
//  normalizerHistogram.apply( nucleusMask );

  OtsuThresholding<float> otsuThresholding;
  otsuThresholding.setForegroundIsAbove( true );
  otsuThresholding.setForeground( 1.0 );
  otsuThresholding.setBackground( 0.0 );
  otsuThresholding.apply( nucleiMask );
  EVAL(nucleiMask.getSize());
//  nucleiMask.save( "/home/jarpon/Desktop/test0.vm", true );

  HolesFillingf holesFilling;

  int sizeZ = originalVoxelMatrix.getSize3();


  ComponentLabelling<float> componentLabelling;
  componentLabelling.apply( nucleiMask );

  const int numComponents = componentLabelling.getNumLabels();


  Thresholding<float> thresholding;
  thresholding.setForegroundIsAbove( true );
  thresholding.setForeground( 1.0 );
  thresholding.setBackground( 0.0 );
//  thresholding.apply( nucleusMask );  //for images already segmented and comment the next threshold

  //a final threshold is applied in each slide to leave only the biggest shape
  //(wich is labeled with 1, the other labeled shapes are set to 0)
  //thresholding.setThreshold(otsuThresholding.computeThreshold(histogram)  );
 // thresholding.setThreshold( 0.5 );

  RegionAnalysis<float> regionAnalysis;
  regionAnalysis.setRegionMatrix( nucleiMask );
  regionAnalysis.run();
  regionAnalysis.mapRegionFeature( nucleiMask, REGION_FEATURE_VOLUME, originalVoxelMatrix );
  nucleiMask.save( "/home/jarpon/Desktop/test1.vm", true );

  // manual threshold -> just looking to the minimum volume of nuclei
  thresholding.setThreshold( 20 );
  thresholding.apply( nucleiMask );
//  nucleiMask.save( "/home/jarpon/Desktop/test2.vm", true );

//  //process to improve nuclei segmentation when the nucleoli touch the envelope
//  VoxelMatrix<float> structElement3;
//  structElement3.setSize(7,7,7);
//  structElement3.setOnes();

//  VoxelMatrixDilatation<float> voxelDilatation2;
//  voxelDilatation2.setStructElt( structElement3 );
//  voxelDilatation2.apply( nucleiMask );

//  VoxelMatrix<float> structElement4;
//  structElement4.setSize(3,3,3);
//  structElement4.setOnes();

//  VoxelMatrixErosion<float> voxelErosion2;
//  voxelErosion2.setStructElt( structElement4 );
//  voxelErosion2.apply( nucleiMask );


  //another 3D filling
  //holesFilling.apply( nucleiMask );


  if ( nucleiMask.max().max().max() > 100 )
  {
      VoxelMatrix<float> gradientMatrix = originalVoxelMatrix;

      GaussianGradient<float> gaussianGradient;
      gaussianGradient.MaskVoxelMatrixProcessing<float>::setMask( nucleiMask );
      gaussianGradient.setSigma( 1.4 );
      gaussianGradient.apply( gradientMatrix );

      VoxelMatrix <float> regionMatrix = gradientMatrix;
      WatershedTransform<float> watershedTransform;
      watershedTransform.MaskVoxelMatrixProcessing<float>::setMask( nucleiMask );
      watershedTransform.apply( regionMatrix );
      regionAnalysis.setRegionMatrix( regionMatrix );
      regionAnalysis.run();
      regionAnalysis.mapRegionFeature( nucleiMask, REGION_FEATURE_CONTRAST, originalVoxelMatrix );

      Vector <float> featureValues = regionAnalysis.computeRegionFeature( REGION_FEATURE_CONTRAST, originalVoxelMatrix );
      Vector <unsigned int> histogram = featureValues.histogram( 0, 1, histogram.getSize()-1 );

      thresholding.setThreshold( otsuThresholding.computeThreshold(histogram) );
      thresholding.apply( nucleiMask );
      //thresholding.applyAlternative( nucleiMask );
  //    nucleiMask.save( "/home/jarpon/Desktop/test3.vm", true );
  }

  nucleiMask.setVoxelCalibration( originalVoxelMatrix.getVoxelCalibration() );
//  VoxelMatrix<float> newOriginal = originalVoxelMatrix;
//  newOriginal.applyMask( originalVoxelMatrix, nucleiMask );
  //newOriginal.save( "/home/jarpon/Desktop/testVM.vm", true );


  LEAVE();
  return nucleiMask;
}
