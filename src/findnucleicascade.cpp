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
#include <image.h>
#include <medianfilter.h>
#include <gaussiangradient.h>
#include <sstream>
#include <marchingcubes.h>
#include <trimesh.h>
#include <regionanalysis.h>
#include <watershedtransform.h>
#include <cmath>

#define TRACE
#include <trace.h>
using namespace std;

//void getMaskTriMesh( VoxelMatrix<float> tempMask )
//{
//  MarchingCubes<float> marchingCubes;
//  TriMesh<float> triMesh;
//  Vector<float> bbVertex1, bbVertex2;
//  triMesh = marchingCubes.buildMesh( tempMask, 0.5, true );
//  triMesh.scale( tempMask.getVoxelCalibration().getVoxelSize() );
//  EVAL (triMesh.boundingBox().getVertex1());
//  EVAL (triMesh.boundingBox().getVertex2());
//  bbVertex1 = triMesh.boundingBox().getVertex1();
//  bbVertex2 = triMesh.boundingBox().getVertex2();

//  Vector<float> diff;

//  bbVertex1.operator -=( 1.2, 1.2, 1.2 );
//  bbVertex2.operator +=( diff );
//  EVAL (bbVertex1);
//  EVAL (bbVertex2);
//}

VoxelMatrix <float>  dilateIt( VoxelMatrix<float> tempMask )
{
  EVAL("dilateIt");
  VoxelMatrix<float> structElement, structElement2;
  structElement.setSize(1,1,1);
  structElement.setOnes();

  VoxelMatrixErosion<float> voxelErosion;
  voxelErosion.setStructElt( structElement );
  voxelErosion.apply( tempMask );

  structElement2.setSize(3,3,7);
  structElement2.setOnes();

  VoxelMatrixDilatation<float> voxelDilatation;
  voxelDilatation.setStructElt( structElement2 );
  voxelDilatation.apply( tempMask );

  voxelErosion.apply( tempMask );
  voxelDilatation.apply( tempMask );

  return tempMask;

}

VoxelMatrix <float> labelIt( VoxelMatrix<float>& dilatedNucleusMask, const VoxelMatrix<float>& originalVoxelMatrix )
{
  EVAL("labelIt");
  VoxelMatrix<float> gradientMatrix = originalVoxelMatrix;
  VoxelMatrix<float> nucleusMaskCopy = dilatedNucleusMask;

  //in case the mask is not binarized
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
  gaussianGradient.setSigma( 5 );
  //gaussianGradient.setSigma( 2 );
  //for (int k = 0; k < sizeZ; ++k) gaussianGradient.apply( gradientMatrix [k] );
  gaussianGradient.apply( gradientMatrix );
  //gradientMatrix.save( intermediateProcessesDir + filename + "-gradient.vm", true );
  //gradientMatrix.operator /=( 16 );

  VoxelMatrix <float> regionMatrix = gradientMatrix;
  WatershedTransform<float> watershedTransform;
  watershedTransform.MaskVoxelMatrixProcessing<float>::setMask( nucleusMaskCopy );
  watershedTransform.apply( regionMatrix );

  return regionMatrix;
}

void applyCascade( VoxelMatrix<float>& voxelMatrixMask, const VoxelMatrix<float>& originalVoxelMatrix )
{
  VoxelMatrix<float> tempMask;
  tempMask = dilateIt( voxelMatrixMask );

  VoxelMatrix<float> regionMatrix;
  regionMatrix = labelIt( tempMask, originalVoxelMatrix );

  VoxelMatrix<float> rangeMask, copyVoxelMatrix = originalVoxelMatrix;

  regionMatrix.save( "/home/jarpon/data/new/nucIntEx0.vm", true );

  RegionAnalysis<float> regionAnalysis;
  regionAnalysis.setRegionMatrix( regionMatrix );
  regionAnalysis.run();
  regionAnalysis.mapRegionFeature( rangeMask, REGION_FEATURE_CONTRAST, copyVoxelMatrix );
  rangeMask.save( "/home/jarpon/data/new/nucIntEx.vm", true );

  Vector <float> featureValues = regionAnalysis.computeRegionFeature( REGION_FEATURE_CONTRAST, copyVoxelMatrix );
  //Vector <float> featureValues = regionAnalysis.computeRegionFeature( REGION_FEATURE_AVERAGE_VALUE, rangeMask );

  Vector <unsigned int> histogram = featureValues.histogram( featureValues.min(), 1, floor(featureValues.max())+1 );
  OtsuThresholding<float> otsuThresholding;
  Thresholding<float> thresholding;
  thresholding.setBackground( 0.0 );
  thresholding.setThreshold( otsuThresholding.computeThreshold(histogram) );
  thresholding.apply( copyVoxelMatrix );
  copyVoxelMatrix.save( "/home/jarpon/data/new/nucTh1.vm", true );
  thresholding.apply( voxelMatrixMask );
  voxelMatrixMask.save( "/home/jarpon/data/new/nucTh2.vm", true );

}


VoxelMatrix<float> findNucleusCascadeMethod(const VoxelMatrix<float>& originalVoxelMatrix,
                                            const string& filename, const string& nucleiDir)
{
    VoxelMatrix<float> nucleusMask = originalVoxelMatrix;

    OtsuThresholding<float> otsuThresholding;
    otsuThresholding.setForegroundIsAbove( true );
    otsuThresholding.setForeground( 1.0 );
    otsuThresholding.setBackground( 0.0 );
    otsuThresholding.apply( nucleusMask );

    const int size1 = nucleusMask.getSize1();
    const int size2 = nucleusMask.getSize2();
    const int size3 = nucleusMask.getSize3();

    HolesFillingf holesFilling;

    //to obtain a better filling, it's applied to each 2D slice instead of the complete 3D stack
    for (int k = 0; k < size3; ++k)  holesFilling.apply( nucleusMask[k] );
    holesFilling.apply( nucleusMask );

    VoxelMatrix<float> nucleusFill = nucleusMask;
    nucleusMask.fillIt(nucleusFill);

    /* label each region
     */
    ComponentLabelling<float> componentLabelling;
    componentLabelling.apply( nucleusMask );

    /* check if some label touches the boundary
     * (that means it's truncated)
     */


    int labelTruncated = 0;

    for ( int i = 0; i < size1; ++i )
      for ( int j = 0; j < size2; ++j )
        for ( int k = 0; k < size3; ++k )
            if ( i==0 || i==size1-1 || j==0 || j==size2-1 || k==0 || k==size3-1 )
              if ( nucleusMask(i,j,k) != 0 )
              {
                  labelTruncated = nucleusMask(i,j,k);
                  EVAL(labelTruncated);
                  for (int ii = 0; ii < size1; ++ii)
                  for (int jj = 0; jj < size2; ++jj)
                    for (int kk = 0; kk < size3; ++kk)
                       if ( nucleusMask(ii,jj,kk) == labelTruncated )
                           nucleusMask(ii,jj,kk) = 0;
              }

    componentLabelling.apply( nucleusMask );
    const int numComponents = componentLabelling.getNumLabels();

    if ( numComponents == 0 )
    {
      nucleusMask = originalVoxelMatrix;
      otsuThresholding.setThreshold( otsuThresholding.getThreshold()*0.8 );
      otsuThresholding.apply( nucleusMask );
    }

    /* check if there are more than one nucleus in the image
     * regarding to the volumes of the different labels
     */

    int numNuclei = 1;

    if ( numComponents > 1 )
    {
      const Vector<unsigned int> histogram = nucleusMask.histogram( 0, 1, numComponents+1 );
      const Vector<unsigned int> histogram2 = histogram.copy( 1, histogram.getSize()-1 );
      EVAL(histogram2);

      Thresholding<float> thresholding;
      thresholding.setForegroundIsAbove( true );
      thresholding.setBackground( 0.0 );
      thresholding.setForeground( 1.0 );

      nucleusMask.setVoxelCalibration( originalVoxelMatrix.getVoxelCalibration() );

      VoxelMatrix<float> tempNucleusMask = nucleusMask;
      for ( int count = 0; count < numComponents; ++count )
      {
          if ( histogram2.max()/histogram2[count] < 4 && histogram2.maxPos() != count && histogram2[count]>5000 )
        //Z  if ( histogram2.max()/histogram2[count] < 4  )
          {
              ostringstream iss;
              iss << numNuclei;
              for (int k = 0; k < size3; ++k)  thresholding.levelSetMask( tempNucleusMask[k], count+1 );

              VoxelMatrix<float> nucleusFill = tempNucleusMask;
              tempNucleusMask.fillIt(nucleusFill);
              //nucleusMask = nucleusFill;

              applyCascade( tempNucleusMask, originalVoxelMatrix );



              tempNucleusMask.save ( nucleiDir + filename + "-" + iss.str() + ".vm", true );
              tempNucleusMask = nucleusMask;
              ++numNuclei;

          }
      }

      for (int k = 0; k < size3; ++k)
          thresholding.levelSetMask( nucleusMask[k], histogram2.maxPos()+1 );

      //process to improve nuclei segmentation when the nucleoli touch the envelope
      VoxelMatrix<float> structElement4;
      structElement4.setSize(1,1,1);
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

      for (int k = 0; k < size3; ++k)  holesFilling.apply( nucleusMask[k] );

      if ( numNuclei != 1 )
      {
          ostringstream iss;
          iss << numNuclei;
          nucleusMask.save ( nucleiDir + filename + "-" + iss.str() + ".vm", true );
      }
      else
      {
        applyCascade( nucleusMask, originalVoxelMatrix );
        nucleusMask.save ( nucleiDir + filename + ".vm", true );
      }
    }

    else if ( numComponents == 1 )
    {
        //process to improve nuclei segmentation when the nucleoli touch the envelope
        VoxelMatrix<float> structElement4;
        structElement4.setSize(1,1,1);
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

        for (int k = 0; k < size3; ++k)  holesFilling.apply( nucleusMask[k] );

        applyCascade( nucleusMask, originalVoxelMatrix );

        nucleusMask.save ( nucleiDir + filename + ".vm", true );
    }

    //another 3D filling
    //for (int k = 0; k < size3; ++k)  holesFilling.apply( nucleusMask[k] );


    return nucleusMask;

}
