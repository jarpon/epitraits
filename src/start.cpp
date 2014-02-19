#include <iostream>
#include <findnucleus.h>
#include <gaussiangradient.h>
#include <otsuthresholding.h>
#include <regionanalysis.h>
#include <stack.h>
#include <stackreader.h>
#include <voxelmatrix.h>
#include <watershedtransform.h>

#define TRACE
#include <trace.h>
using namespace std;

const string casename = "20130121_s001_001.lsm";
const string inputDir = "/home/arpon/Desktop/examples/originals/";
const string outputDir = "/home/arpon/Desktop/examples/processed/";




void newProtocol()
{
  //Stack stack;
  //StackReader stackReader;

  //stackReader.read( stack, inputDir+casename + ".sdf" );

 // EVAL(inputDir + casename + ".sdf");

  //stackReader.read( stack, inputDir+casename+ ".lsm" );

    VoxelMatrix<float> voxelMatrix(inputDir + casename);
    //voxelMatrix.load(inputDir + casename);
    //VoxelMatrix<float> findNucleus;
    FindNucleus<float> findNucleus;
    findNucleus.apply( voxelMatrix,  inputDir,  casename,  outputDir );


  //EVAL( voxelMatrix.getSize() );
  //voxelMatrix.save( outputDir+"original", true );

  return;



  //nucleusMask.save( outputDir+casename+".vm", true );

#if 0
  VoxelMatrix<float> gradientMatrix = voxelMatrix;
  GaussianGradient<float> gaussianGradient;
  gaussianGradient.MaskVoxelMatrixProcessing<float>::setMask( nucleusMask );
  gaussianGradient.setSigma( 3.0 );
  gaussianGradient.apply( gradientMatrix );
  gradientMatrix.save( outputDir+casename+"-gradient.vm", true );

  WatershedTransform<float> watershedTransform;
  watershedTransform.MaskVoxelMatrixProcessing<float>::setMask( nucleusMask );
  watershedTransform.apply( gradientMatrix );
  gradientMatrix.save( outputDir+casename+"-watershed.vm", true );

  VoxelMatrix<float> featureMatrix;
  RegionAnalysis<float> regionAnalysis;
  regionAnalysis.setRegionMatrix( gradientMatrix );
  regionAnalysis.run();
  regionAnalysis.mapRegionFeature( featureMatrix, REGION_FEATURE_AVERAGE_VALUE, voxelMatrix );
  featureMatrix.save( outputDir+casename+"-average.vm", true );

  OtsuThresholding<float> otsuThresholding;
  Vector<float> featureValues = regionAnalysis.allRegionValues( featureMatrix );
  Vector<unsigned int> histogram = featureValues.histogram( 0, 1, 256 );
  EVAL( otsuThresholding.computeThreshold( histogram ) );

  Thresholding<float> thresholding;
  thresholding.setBackground( 0.0 );
  thresholding.setForeground( 255.0 );
  thresholding.setThreshold( otsuThresholding.computeThreshold(histogram) );
  thresholding.apply( featureMatrix );
  featureMatrix.save( outputDir+casename+"-cc.vm", true );
#endif
}



int main()
{

    newProtocol();

//    Stack stack;
//    StackReader stackReader;

//    stackReader.read( stack, inputDir+casename+".sdf" );

//    VoxelMatrix<float> voxelMatrix;
//    voxelMatrix = stack.toVoxelMatrix<float>( IMAGE_COLOR_CONTENT_GRAY );

//    VoxelMatrix<float> nucleusMask = findNucleus( voxelMatrix );
//    nucleusMask.save( outputDir+casename+".vm", true );

}

