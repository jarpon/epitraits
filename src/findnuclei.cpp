#include <componentlabelling.h>
#include <holesfilling.h>
#include "otsuthresholding.h"
#include "thresholding.h"
#include "voxelmatrix.h"
#include <volumehistogramexpansion.h>
//#include <maxfilter.h>
//#include <minfilter.h>
#include <voxelmatrixerosion.h>
#include <voxelmatrixdilatation.h>
#include <image.h>
#include <medianfilter.h>
#include <gaussiangradient.h>

#define TRACE
#include <trace.h>
using namespace std;

VoxelMatrix<float> findNucleus(const VoxelMatrix<float>& originalVoxelMatrix)
{
  VoxelMatrix<float> nucleusMask = originalVoxelMatrix;


  MedianFilter<float> medianFilter2;
  medianFilter2.setHalfSize( 2 );
  medianFilter2.apply( nucleusMask );


  OtsuThresholding<float> otsuThresholding;
  otsuThresholding.setForegroundIsAbove( true );
  otsuThresholding.setForeground( 1.0 );
  otsuThresholding.setBackground( 0.0 );
  otsuThresholding.apply( nucleusMask );


// // for originals in 16bits
//  Thresholding<float> thresholding2;
//  thresholding2.setForegroundIsAbove( true );
//  thresholding2.setForeground( 1.0 );
//  thresholding2.setBackground( 0.0 );
// // thresholding2.setThreshold( 300.0 );
//  thresholding2.setThreshold( nucleusMask.max().max().max()/10 );
//  thresholding2.apply( nucleusMask );
//  //nucleusMask.save ( "/home/javier/Desktop/gaussian.vm", true );

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
  if ( numComponents == 0 )
  {
    nucleusMask = originalVoxelMatrix;
    otsuThresholding.setThreshold( otsuThresholding.getThreshold()*0.8 );
    otsuThresholding.apply( nucleusMask );
  }

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

 //this was uncommented 12 june 15
//  VoxelMatrix<float> nucleusFill = nucleusMask;
//  nucleusMask.fillIt(nucleusFill);
//  //nucleusMask = nucleusFill;


  //process to improve nuclei segmentation when the nucleoli touch the envelope
  VoxelMatrix<float> structElement4;
//  structElement4.setSize(3,3,3);
  structElement4.setSize(13,13,3);
  structElement4.setOnes();

  VoxelMatrixErosion<float> voxelErosion2;
  voxelErosion2.setStructElt( structElement4 );
  voxelErosion2.apply( nucleusMask );

  VoxelMatrix<float> structElement3;
  structElement3.setSize(3,3,3);
  structElement3.setOnes();

  VoxelMatrixDilatation<float> voxelDilatation2;
  voxelDilatation2.setStructElt( structElement3 );
  voxelDilatation2.apply( nucleusMask );

  //16bits
//  VoxelMatrix<float> structElement3;
//  structElement3.setSize(9,9,9);
//  structElement3.setOnes();

//  VoxelMatrixDilatation<float> voxelDilatation2;
//  voxelDilatation2.setStructElt( structElement3 );
//  voxelDilatation2.apply( nucleusMask );

//  VoxelMatrix<float> structElement4;
//  structElement4.setSize(9,9,9);
//  structElement4.setOnes();

//  VoxelMatrixErosion<float> voxelErosion2;
//  voxelErosion2.setStructElt( structElement4 );
//  voxelErosion2.apply( nucleusMask );


  holesFilling.apply( nucleusMask );
 //this was uncommented 12 june 15
//  nucleusMask.fillIt(nucleusMask);
  nucleusMask.setVoxelCalibration( originalVoxelMatrix.getVoxelCalibration() );

  return nucleusMask;
}
