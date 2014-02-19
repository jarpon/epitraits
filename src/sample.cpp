//#include <gaussiangradient.h>
//#include <otsuthresholding.h>
//#include <regionanalysis.h>
//#include <stack.h>
//#include <stackreader.h>
//#include <watershedtransform.h>

//#define TRACE
//#include <trace.h>

////const string casename = "image5";
////const string inputDir = "/home/andrey/lab/valerieGaudin/temp/stacks/2012.04.20/sdf/";
////const string outputDir = "/home/andrey/lab/valerieGaudin/temp/results/";

//const string casename = "20130129_s018_032";
//const string inputDir = "/home/arpon/Desktop/Test_Javier/dataTest/";
//const string outputDir = "/home/arpon/Desktop/Test_Javier/";

//VoxelMatrix<float> findNucleus(const VoxelMatrix<float>& voxelMatrix)
//{
//  VoxelMatrix<float> nucleusMask = voxelMatrix;
//  OtsuThresholding<float> otsuThresholding;
//  otsuThresholding.setForegroundIsAbove( true );
//  otsuThresholding.setForeground( 1.0 );
//  otsuThresholding.setBackground( 0.0 );
//  otsuThresholding.apply( nucleusMask );
//  return nucleusMask;
//}

//void newProtocol()
//{
//  Stack stack;
//  StackReader stackReader;

//  stackReader.read( stack, inputDir+casename+".sdf" );

//  //stackReader.read( stack, inputDir+casename+ ".lsm" );

//  VoxelMatrix<float> voxelMatrix;
//  voxelMatrix = stack.toVoxelMatrix<float>( IMAGE_COLOR_CONTENT_GRAY );

//  VoxelMatrix<float> nucleusMask = findNucleus( voxelMatrix );
//  nucleusMask.save( outputDir+casename+".vm", true );

//#if 0
//  VoxelMatrix<float> gradientMatrix = voxelMatrix;
//  GaussianGradient<float> gaussianGradient;
//  gaussianGradient.MaskVoxelMatrixProcessing<float>::setMask( nucleusMask );
//  gaussianGradient.setSigma( 3.0 );
//  gaussianGradient.apply( gradientMatrix );
//  gradientMatrix.save( outputDir+casename+"-gradient.vm", true );

//  WatershedTransform<float> watershedTransform;
//  watershedTransform.MaskVoxelMatrixProcessing<float>::setMask( nucleusMask );
//  watershedTransform.apply( gradientMatrix );
//  gradientMatrix.save( outputDir+casename+"-watershed.vm", true );

//  VoxelMatrix<float> featureMatrix;
//  RegionAnalysis<float> regionAnalysis;
//  regionAnalysis.setRegionMatrix( gradientMatrix );
//  regionAnalysis.run();
//  regionAnalysis.mapRegionFeature( featureMatrix, REGION_FEATURE_AVERAGE_VALUE, voxelMatrix );
//  featureMatrix.save( outputDir+casename+"-average.vm", true );

//  OtsuThresholding<float> otsuThresholding;
//  Vector<float> featureValues = regionAnalysis.allRegionValues( featureMatrix );
//  Vector<unsigned int> histogram = featureValues.histogram( 0, 1, 256 );
//  EVAL( otsuThresholding.computeThreshold( histogram ) );

//  Thresholding<float> thresholding;
//  thresholding.setBackground( 0.0 );
//  thresholding.setForeground( 255.0 );
//  thresholding.setThreshold( otsuThresholding.computeThreshold(histogram) );
//  thresholding.apply( featureMatrix );
//  featureMatrix.save( outputDir+casename+"-cc.vm", true );
//#endif
//}
