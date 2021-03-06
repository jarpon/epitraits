//#include <componentlabelling.h>
//#include <holesfilling.h>
//#include "otsuthresholding.h"
//#include "thresholding.h"
//#include "voxelmatrix2.h"
//#include <volumehistogramexpansion.h>
////#include <maxfilter.h>
////#include <minfilter.h>
//#include <voxelmatrixerosion.h>
//#include <voxelmatrixdilatation.h>
//#include <image.h>
//#include <medianfilter.h>
//#include <gaussiangradient.h>
//#include <sstream>

//#define TRACE
//#include <trace.h>
//using namespace std;

//VoxelMatrix<float> findNucleus(const VoxelMatrix<float>& originalVoxelMatrix)
//{
//  VoxelMatrix<float> nucleusMask = originalVoxelMatrix;
//  string nucleiDir = "/media/jarpon/WD-MIN-1TB/data-javier/maria/more/c0/segmented_nuclei/MK67_2 leaf002_1_c0.tif.vm";
//  string filename = "MK67_2 leaf002_1_c0.tif.vm";

////  MedianFilter<float> medianFilter2;
////  medianFilter2.setHalfSize( 2 );
////  medianFilter2.apply( nucleusMask );


//  OtsuThresholding<float> otsuThresholding;
//  otsuThresholding.setForegroundIsAbove( true );
//  otsuThresholding.setForeground( 1.0 );
//  otsuThresholding.setBackground( 0.0 );
//  otsuThresholding.apply( nucleusMask );

//  const int size1 = nucleusMask.getSize1();
//  const int size2 = nucleusMask.getSize2();
//  const int size3 = nucleusMask.getSize3();

//  HolesFillingf holesFilling;

//  //to obtain a better filling, it's applied to each 2D slice instead of the complete 3D stack
//  for (int k = 0; k < size3; ++k)  holesFilling.apply( nucleusMask[k] );
//  holesFilling.apply( nucleusMask );

////this was uncommented 12 june
////  VoxelMatrix<float> nucleusFill = nucleusMask;
////  nucleusMask.fillIt(nucleusFill);

////  /*process to improve nuclei segmentation when the nucleoli touch the envelope
////   */
////  VoxelMatrix<float> structElement;
////  structElement.setSize(5,5,5);
////  structElement.setOnes();
////  //different structure element due to problems of segmentation (a kind of smoothing)
////  VoxelMatrix<float> structElement2;
////  structElement2.setSize(3,3,3);
////  structElement2.setOnes();

////  VoxelMatrixDilatation<float> voxelDilatation;
////  voxelDilatation.setStructElt( structElement );
////  voxelDilatation.apply( nucleusMask );
////  VoxelMatrixErosion<float> voxelErosion;
////  voxelErosion.setStructElt( structElement2 );
////  voxelErosion.apply( nucleusMask );

////  holesFilling.apply( nucleusMask );
////  nucleusMask.fillIt(nucleusMask);

//  /* label each region
//   */
//  ComponentLabelling<float> componentLabelling;
//  componentLabelling.apply( nucleusMask );

//  /* check if some label touches the boundary
//   * (that means it's truncated)
//   */


//  int labelTruncated = 0;

//  for ( int i = 0; i < size1; ++i )
//    for ( int j = 0; j < size2; ++j )
//      for ( int k = 0; k < size3; ++k )
//          if ( i==0 || i==size1-1 || j==0 || j==size2-1 || k==0 || k==size3-1 )
//            if ( nucleusMask(i,j,k) != 0 )
//            {
//                labelTruncated = nucleusMask(i,j,k);
//                EVAL(labelTruncated);
//                for (int ii = 0; ii < size1; ++ii)
//                for (int jj = 0; jj < size2; ++jj)
//                  for (int kk = 0; kk < size3; ++kk)
//                     if ( nucleusMask(ii,jj,kk) == labelTruncated )
//                         nucleusMask(ii,jj,kk) = 0;
//            }

//  componentLabelling.apply( nucleusMask );
//  const int numComponents = componentLabelling.getNumLabels();

//  if ( numComponents == 0 )
//  {
//    nucleusMask = originalVoxelMatrix;
//    otsuThresholding.setThreshold( otsuThresholding.getThreshold()*0.8 );
//    otsuThresholding.apply( nucleusMask );
//  }

//  /* check if there are more than one nucleus in the image
//   * regarding to the volumes of the different labels
//   */

//  int numNuclei = 1;

//  if ( numComponents > 1 )
//  {
//    const Vector<unsigned int> histogram = nucleusMask.histogram( 0, 1, numComponents+1 );
//    const Vector<unsigned int> histogram2 = histogram.copy( 1, histogram.getSize()-1 );
//    EVAL(histogram2);

//    Thresholding<float> thresholding;
//    thresholding.setForegroundIsAbove( true );
//    thresholding.setBackground( 0.0 );
//    thresholding.setForeground( 1.0 );

//    nucleusMask.setVoxelCalibration( originalVoxelMatrix.getVoxelCalibration() );

//    VoxelMatrix<float> tempNucleusMask = nucleusMask;
//    for ( int count = 0; count < numComponents; ++count )
//    {
//        if ( histogram2.max()/histogram2[count] < 4 && histogram2.maxPos() != count && histogram2[count]>5000 )
//      //Z  if ( histogram2.max()/histogram2[count] < 4  )
//        {
//            ostringstream iss;
//            iss << numNuclei;
//            for (int k = 0; k < size3; ++k)  thresholding.levelSetMask( tempNucleusMask[k], count+1 );

//  //this was uncommented 12 june
//        //    VoxelMatrix<float> nucleusFill = tempNucleusMask;
//        //    tempNucleusMask.fillIt(nucleusFill);
//            //nucleusMask = nucleusFill;


//            //process to improve nuclei segmentation when the nucleoli touch the envelope
//            VoxelMatrix<float> structElement4;
//            //structElement4.setSize(3,3,3);
//            structElement4.setSize(5,5,5);
//            structElement4.setOnes();

//            VoxelMatrixErosion<float> voxelErosion2;
//            voxelErosion2.setStructElt( structElement4 );
//            voxelErosion2.apply( tempNucleusMask );

//            VoxelMatrix<float> structElement3;
//            structElement3.setSize(3,3,3);
//            structElement3.setOnes();

//            VoxelMatrixDilatation<float> voxelDilatation2;
//            voxelDilatation2.setStructElt( structElement3 );
//            voxelDilatation2.apply( tempNucleusMask );

//            for (int k = 0; k < size3; ++k)  holesFilling.apply( tempNucleusMask[k] );

//            tempNucleusMask.save ( nucleiDir + filename + "-" + iss.str() + ".vm", true );
//            tempNucleusMask = nucleusMask;
//            ++numNuclei;

//        }
//    }

//    for (int k = 0; k < size3; ++k)
//        thresholding.levelSetMask( nucleusMask[k], histogram2.maxPos()+1 );

//    //process to improve nuclei segmentation when the nucleoli touch the envelope
//    VoxelMatrix<float> structElement4;
//    //structElement4.setSize(3,3,3);
//    structElement4.setSize(5,5,5);
//    structElement4.setOnes();

//    VoxelMatrixErosion<float> voxelErosion2;
//    voxelErosion2.setStructElt( structElement4 );
//    voxelErosion2.apply( nucleusMask );

//    VoxelMatrix<float> structElement3;
//    structElement3.setSize(3,3,3);
//    structElement3.setOnes();

//    VoxelMatrixDilatation<float> voxelDilatation2;
//    voxelDilatation2.setStructElt( structElement3 );
//    voxelDilatation2.apply( nucleusMask );

//    for (int k = 0; k < size3; ++k)  holesFilling.apply( nucleusMask[k] );

//    if ( numNuclei != 1 )
//    {
//        ostringstream iss;
//        iss << numNuclei;
//        nucleusMask.save ( nucleiDir + filename + "-" + iss.str() + ".vm", true );
//    }
//    else
//        nucleusMask.save ( nucleiDir + filename + ".vm", true );
//  }

//  else if ( numComponents == 1 )
//  {
//      //process to improve nuclei segmentation when the nucleoli touch the envelope
//      VoxelMatrix<float> structElement4;
//      //structElement4.setSize(3,3,3);
//      structElement4.setSize(5,5,5);
//      structElement4.setOnes();

//      VoxelMatrixErosion<float> voxelErosion2;
//      voxelErosion2.setStructElt( structElement4 );
//      voxelErosion2.apply( nucleusMask );

//      VoxelMatrix<float> structElement3;
//      structElement3.setSize(3,3,3);
//      structElement3.setOnes();

//      VoxelMatrixDilatation<float> voxelDilatation2;
//      voxelDilatation2.setStructElt( structElement3 );
//      voxelDilatation2.apply( nucleusMask );

//      for (int k = 0; k < size3; ++k)  holesFilling.apply( nucleusMask[k] );

//      nucleusMask.save ( nucleiDir + filename + ".vm", true );
//  }

//  //another 3D filling
//  //for (int k = 0; k < size3; ++k)  holesFilling.apply( nucleusMask[k] );


//  return nucleusMask;
//}
