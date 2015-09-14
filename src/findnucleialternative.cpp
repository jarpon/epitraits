#include <componentlabelling.h>
#include <holesfilling.h>
#include "otsuthresholding.h"
#include "thresholding.h"
#include "voxelmatrix2.h"
#include <volumehistogramexpansion.h>
#include <voxelmatrixerosion.h>
#include <voxelmatrixdilatation.h>
#include <image.h>
#include <medianfilter.h>
#include <gaussiangradient.h>
#include "regionanalysis.h"
#include <watershedtransform.h>
#include <gaussiangradient.h>

#define TRACE
#include <trace.h>
using namespace std;

/* using Z projection
VoxelMatrix<float> findNucleusAlternative(const VoxelMatrix<float>& originalVoxelMatrix)
{
  VoxelMatrix<float> nucleusMask = originalVoxelMatrix;
  VoxelMatrix<float> nucleusZprojectionVol = originalVoxelMatrix;
  PixelMatrix<float> nucleusZprojection;
  Image image;
  nucleusZprojection = nucleusMask.getZProjection();
  image.make32bits();
  nucleusZprojection.convertToImage(image);
  image.writeFile("nucleus","/home/javier/Desktop/",true);

//  MedianFilter<float> medianFilter2;
//  medianFilter2.setHalfSize( 2 );
//  medianFilter2.apply( nucleusZprojection );

  for ( int z = 0; z < originalVoxelMatrix.getSize3(); ++z )
    nucleusZprojectionVol[z] = nucleusZprojection;


  Thresholding<float> thresholding;
  thresholding.setBackground( 0.0 );
  thresholding.setForeground( 1.0 );
  thresholding.apply( nucleusZprojectionVol );
  nucleusZprojectionVol.save ( "/home/javier/Desktop/22.vm", true );

  VoxelMatrix<float> structElement;
  structElement.setSize(5,5,5);
  structElement.setOnes();

  VoxelMatrixDilatation<float> voxelDilatation;
  voxelDilatation.setStructElt( structElement );
  voxelDilatation.apply( nucleusZprojectionVol );

  VoxelMatrix <float> regionMatrix = nucleusMask;
  WatershedTransform<float> watershedTransform;
  watershedTransform.MaskVoxelMatrixProcessing<float>::setMask( nucleusZprojectionVol );
  watershedTransform.apply( regionMatrix );
  regionMatrix.save ( "/home/javier/Desktop/gaussian2.vm", true );

  RegionAnalysis<float> regionAnalysis;
  regionAnalysis.setRegionMatrix( regionMatrix );
  regionAnalysis.run();
  regionAnalysis.mapRegionFeature( nucleusMask, REGION_FEATURE_MEDIAN_VALUE, originalVoxelMatrix );
  nucleusMask.save( "/home/javier/Desktop/nucleus-contrast.vm", true );

  Vector <float> featureValues = regionAnalysis.computeRegionFeature( REGION_FEATURE_MEDIAN_VALUE, originalVoxelMatrix );
  Vector <unsigned int> histogram = featureValues.histogram( 0, 1, histogram.getSize()-1 );


  OtsuThresholding<float> otsuThresholding;

  otsuThresholding.setForegroundIsAbove( true );
  otsuThresholding.setForeground( 1.0 );
  otsuThresholding.setBackground( 0.0 );
  otsuThresholding.setThreshold( otsuThresholding.computeThreshold( regionMatrix ) );
  otsuThresholding.apply( nucleusMask );
  //otsuThresholding.apply(nucleusZprojection);
  regionMatrix.save( "/home/javier/Desktop/nucleus-contrast.vm", true );

// // for originals in 16bits
//  Thresholding<float> thresholding2;
//  thresholding2.setForegroundIsAbove( true );
//  thresholding2.setForeground( 1.0 );
//  thresholding2.setBackground( 0.0 );
// // thresholding2.setThreshold( 300.0 );
//  thresholding2.setThreshold( nucleusMask.max().max().max()/10 );
//  thresholding2.apply( nucleusMask );
//  nucleusMask.save ( "/home/javier/Desktop/gaussian.vm", true );

  HolesFillingf holesFilling;
  int sizeZ = originalVoxelMatrix.getSize3();

  //to obtain a better filling, it's applied to each 2D slice instead of the complete 3D stack
  for (int k = 0; k < sizeZ; ++k)  holesFilling.apply( nucleusMask[k] );
  holesFilling.apply( nucleusMask );
  //for (int k = 0; k < sizeZ; ++k)  holesFilling.apply( nucleusMask[k] );

  ComponentLabelling<float> componentLabelling;
  componentLabelling.apply( nucleusMask );
  //ComponentLabelling<float> componentLabelling2;
  //componentLabelling2.apply( nucleusZprojection );

  const int numComponents = componentLabelling.getNumLabels();
  EVAL(numComponents);
  if ( numComponents > 1 )
  {
    const Vector<unsigned int> histogram = nucleusMask.histogram( 0, 1, numComponents+1 );
    const Vector<unsigned int> histogram2 = histogram.copy( 1, histogram.getSize()-1 );

    Thresholding<float> thresholding;
    thresholding.setForegroundIsAbove( true );
  //  thresholding.setThreshold( 1.0 );  //for images already segmented
    thresholding.setForeground( 1.0 );
    thresholding.setBackground( 0.0 );
  //  thresholding.apply( nucleusMask );  //for images already segmented and comment the next threshold

    //a final threshold is applied in each slide to leave only the biggest shape
    //(wich is labeled with 1, the other labeled shapes are set to 0)
    for (int k = 0; k < sizeZ; ++k)  thresholding.levelSetMask( nucleusMask[k], histogram2.maxPos()+1 );
  }

//  MedianFilter<float> medianFilter;
//  medianFilter.setHalfSize( 2 );
//  //medianFilter.setNumIterations( 2 );
//  for (int k = 0; k < sizeZ; ++k) medianFilter.apply( nucleusMask[k] );
//  medianFilter.apply(nucleusZprojection);
//  holesFilling.apply(nucleusZprojection);
//  nucleusMask.save ( "/home/javier/Desktop/gaussian2.vm", true );

  //nucleusZprojection.convertToImage(image);
  //image.writeFile("nucleusSegmented","/home/javier/Desktop/",true);


//  VoxelMatrix<float> structElement1;
//  structElement1.setSize(1,1,1);
//  structElement1.setOnes();

//  VoxelMatrixDilatation<float> voxelDilatation3;
//  voxelDilatation3.setStructElt( structElement1 );
//  voxelDilatation3.apply( nucleusMask );

  //another 3D filling
  //for (int k = 0; k < sizeZ; ++k)  holesFilling.apply( nucleusMask[k] );
  VoxelMatrix<float> nucleusFill = nucleusMask;
  nucleusMask.fillIt(nucleusFill);
  //nucleusMask = nucleusFill;
  nucleusFill.save ( "/home/javier/Desktop/gaussian2.vm", true );

  nucleusMask.save ( "/home/javier/Desktop/gaussian3.vm", true );

//  //process to improve nuclei segmentation when the nucleoli touch the envelope
  VoxelMatrix<float> structElement3;
  structElement3.setSize(5,5,5);
  structElement3.setOnes();

  VoxelMatrixDilatation<float> voxelDilatation2;
  voxelDilatation2.setStructElt( structElement3 );
  voxelDilatation2.apply( nucleusMask );

  VoxelMatrix<float> structElement4;
  structElement4.setSize(3,3,3);
  structElement4.setOnes();

  VoxelMatrixErosion<float> voxelErosion2;
  voxelErosion2.setStructElt( structElement4 );
  voxelErosion2.apply( nucleusMask );

  holesFilling.apply( nucleusMask );
  nucleusMask.fillIt(nucleusMask);
  nucleusMask.setVoxelCalibration( originalVoxelMatrix.getVoxelCalibration() );

  return nucleusMask;
}
*/

VoxelMatrix<float> findNucleusAlternative(const VoxelMatrix<float>& originalVoxelMatrix)
{
  VoxelMatrix<float> nucleusMask = originalVoxelMatrix;

  MedianFilter<float> medianFilter2;
  medianFilter2.setHalfSize( 2 );
  medianFilter2.apply( nucleusMask );

  OtsuThresholding<float> otsuThresholding;
  Thresholding<float> thresholding;
  thresholding.setBackground( 0.0 );
  thresholding.setForeground( 1.0 );
  thresholding.setThreshold( otsuThresholding.computeThreshold( nucleusMask ) );
  thresholding.apply( nucleusMask );
  nucleusMask.save ( "/home/javier/Desktop/nucleus.vm", true );

  VoxelMatrix<float> structElement;
  structElement.setSize(9,9,9);
  structElement.setOnes();

  VoxelMatrix<float> dilated = nucleusMask;
  VoxelMatrix<float> eroded = nucleusMask;

  VoxelMatrixDilatation<float> voxelDilatation;
  voxelDilatation.setStructElt( structElement );
  voxelDilatation.apply( dilated );
  dilated.save ( "/home/javier/Desktop/dilated.vm", true );

  VoxelMatrixErosion<float> voxelErosion;
  voxelErosion.setStructElt( structElement );
  voxelErosion.apply( eroded );
  eroded.save ( "/home/javier/Desktop/eroded.vm", true );

  VoxelMatrix<float> tophat = dilated;
  tophat.operator -=(eroded);
  tophat.save ( "/home/javier/Desktop/tophat.vm", true );

//  VoxelMatrix <float> regionMatrix = nucleusMask;
//  WatershedTransform<float> watershedTransform;
//  watershedTransform.MaskVoxelMatrixProcessing<float>::setMask( nucleusZprojectionVol );
//  watershedTransform.apply( regionMatrix );
//  regionMatrix.save ( "/home/javier/Desktop/gaussian2.vm", true );

//  RegionAnalysis<float> regionAnalysis;
//  regionAnalysis.setRegionMatrix( regionMatrix );
//  regionAnalysis.run();
//  regionAnalysis.mapRegionFeature( nucleusMask, REGION_FEATURE_MEDIAN_VALUE, originalVoxelMatrix );
//  nucleusMask.save( "/home/javier/Desktop/nucleus-contrast.vm", true );

//  Vector <float> featureValues = regionAnalysis.computeRegionFeature( REGION_FEATURE_MEDIAN_VALUE, originalVoxelMatrix );
//  Vector <unsigned int> histogram = featureValues.histogram( 0, 1, histogram.getSize()-1 );


  OtsuThresholding<float> otsuThresholding2;

  otsuThresholding2.setForegroundIsAbove( true );
  otsuThresholding2.setForeground( 1.0 );
  otsuThresholding2.setBackground( 0.0 );

  otsuThresholding2.apply( nucleusMask );
  //otsuThresholding.apply(nucleusZprojection);
  nucleusMask.save( "/home/javier/Desktop/nucleus-contrast.vm", true );

// // for originals in 16bits
//  Thresholding<float> thresholding2;
//  thresholding2.setForegroundIsAbove( true );
//  thresholding2.setForeground( 1.0 );
//  thresholding2.setBackground( 0.0 );
// // thresholding2.setThreshold( 300.0 );
//  thresholding2.setThreshold( nucleusMask.max().max().max()/10 );
//  thresholding2.apply( nucleusMask );
//  nucleusMask.save ( "/home/javier/Desktop/gaussian.vm", true );

  HolesFillingf holesFilling;
  int sizeZ = originalVoxelMatrix.getSize3();

  //to obtain a better filling, it's applied to each 2D slice instead of the complete 3D stack
  for (int k = 0; k < sizeZ; ++k)  holesFilling.apply( nucleusMask[k] );
  holesFilling.apply( nucleusMask );
  //for (int k = 0; k < sizeZ; ++k)  holesFilling.apply( nucleusMask[k] );

  ComponentLabelling<float> componentLabelling;
  componentLabelling.apply( nucleusMask );
  //ComponentLabelling<float> componentLabelling2;
  //componentLabelling2.apply( nucleusZprojection );

  const int numComponents = componentLabelling.getNumLabels();
  EVAL(numComponents);
  if ( numComponents > 1 )
  {
    const Vector<unsigned int> histogram = nucleusMask.histogram( 0, 1, numComponents+1 );
    const Vector<unsigned int> histogram2 = histogram.copy( 1, histogram.getSize()-1 );

    Thresholding<float> thresholding;
    thresholding.setForegroundIsAbove( true );
  //  thresholding.setThreshold( 1.0 );  //for images already segmented
    thresholding.setForeground( 1.0 );
    thresholding.setBackground( 0.0 );
  //  thresholding.apply( nucleusMask );  //for images already segmented and comment the next threshold

    //a final threshold is applied in each slide to leave only the biggest shape
    //(wich is labeled with 1, the other labeled shapes are set to 0)
    for (int k = 0; k < sizeZ; ++k)  thresholding.levelSetMask( nucleusMask[k], histogram2.maxPos()+1 );
  }

//  MedianFilter<float> medianFilter;
//  medianFilter.setHalfSize( 2 );
//  //medianFilter.setNumIterations( 2 );
//  for (int k = 0; k < sizeZ; ++k) medianFilter.apply( nucleusMask[k] );
//  medianFilter.apply(nucleusZprojection);
//  holesFilling.apply(nucleusZprojection);
//  nucleusMask.save ( "/home/javier/Desktop/gaussian2.vm", true );

  //nucleusZprojection.convertToImage(image);
  //image.writeFile("nucleusSegmented","/home/javier/Desktop/",true);


//  VoxelMatrix<float> structElement1;
//  structElement1.setSize(1,1,1);
//  structElement1.setOnes();

//  VoxelMatrixDilatation<float> voxelDilatation3;
//  voxelDilatation3.setStructElt( structElement1 );
//  voxelDilatation3.apply( nucleusMask );

  //another 3D filling
  //for (int k = 0; k < sizeZ; ++k)  holesFilling.apply( nucleusMask[k] );
  VoxelMatrix<float> nucleusFill = nucleusMask;
  
// this was uncommented 12 june
//  nucleusMask.fillIt(nucleusFill);
//  //nucleusMask = nucleusFill;
//  nucleusFill.save ( "/home/javier/Desktop/gaussian2.vm", true );

  //nucleusMask.save ( "/home/javier/Desktop/gaussian3.vm", true );

//  //process to improve nuclei segmentation when the nucleoli touch the envelope
  VoxelMatrix<float> structElement3;
  structElement3.setSize(5,5,5);
  structElement3.setOnes();

  VoxelMatrixDilatation<float> voxelDilatation2;
  voxelDilatation2.setStructElt( structElement3 );
  voxelDilatation2.apply( nucleusMask );

  VoxelMatrix<float> structElement4;
  structElement4.setSize(3,3,3);
  structElement4.setOnes();

  VoxelMatrixErosion<float> voxelErosion2;
  voxelErosion2.setStructElt( structElement4 );
  voxelErosion2.apply( nucleusMask );

  holesFilling.apply( nucleusMask );
// this was uncommented 12 june
//  nucleusMask.fillIt(nucleusMask);
  nucleusMask.setVoxelCalibration( originalVoxelMatrix.getVoxelCalibration() );

  return nucleusMask;
}
