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
//#include "regionanalysis.h"
#include <regionanalysis3d.h>
#include <watershedtransform.h>
#include <gaussiangradient.h>


#define TRACE
#include <trace.h>
using namespace std;

extern VoxelMatrix<float> applyLabelling(const VoxelMatrix<float>&, const VoxelMatrix<float>&, const string&, const string&);

VoxelMatrix<float> isolateNuclei(const VoxelMatrix<float>& originalVoxelMatrix, const string& filename, const string& parentDir)
{
  ENTER("isolateNuclei");

  VoxelMatrix<float> nucleiMask = originalVoxelMatrix, regionMatrix;

  Thresholding<float> thresholding;
  thresholding.setForegroundIsAbove( true );
  thresholding.setForeground( 1.0 );
  thresholding.setBackground( 0.0 );
  thresholding.apply( nucleiMask );


  HolesFillingf holesFilling;

  int sizeZ = originalVoxelMatrix.getSize3();
  for (int k = 0; k < sizeZ; ++k)  holesFilling.apply( nucleiMask[k] );

  holesFilling.apply( nucleiMask );

  nucleiMask.save( parentDir + "/intermediate_processes/" + filename + "-mask.vm", true );

  regionMatrix = applyLabelling( originalVoxelMatrix, nucleiMask, filename, parentDir + "/intermediate_processes/" );

  EVAL("1");
//  RegionAnalysis3D<float> regionAnalysis;
//  regionAnalysis.setLabelMatrix( regionMatrix );
//  regionAnalysis.setValueMatrix( originalVoxelMatrix );
//  regionAnalysis.run();

//  nucleiMask.setSize( copyVoxelMatrix.getSize() );
//  nucleiMask.setZeros();
//  regionAnalysis.setOutputMatrix( nucleiMask );

  EVAL("1");
  //uncomment - 15 july
  //regionAnalysis.mapRegionFeature( nucleiMask, REGION_FEATURE_CONTRAST, originalVoxelMatrix );

  nucleiMask.save( parentDir + "/intermediate_processes/" + filename + ".vm", true );

  OtsuThresholding<float> otsuThresholding;
  otsuThresholding.setBackground( 0.0 );
  otsuThresholding.applyAlternative( nucleiMask );



//  if ( nucleiMask.max().max().max() > 100 )
//  {
//      VoxelMatrix<float> gradientMatrix = originalVoxelMatrix;

//      GaussianGradient<float> gaussianGradient;
//      gaussianGradient.MaskVoxelMatrixProcessing<float>::setMask( nucleiMask );
//      gaussianGradient.setSigma( 1.4 );
//      gaussianGradient.apply( gradientMatrix );

//      VoxelMatrix <float> regionMatrix = gradientMatrix;
//      WatershedTransform<float> watershedTransform;
//      watershedTransform.MaskVoxelMatrixProcessing<float>::setMask( nucleiMask );
//      watershedTransform.apply( regionMatrix );
//      regionAnalysis.setRegionMatrix( regionMatrix );
//      regionAnalysis.run();
//      regionAnalysis.mapRegionFeature( nucleiMask, REGION_FEATURE_CONTRAST, originalVoxelMatrix );

//      Vector <float> featureValues = regionAnalysis.computeRegionFeature( REGION_FEATURE_CONTRAST, originalVoxelMatrix );
//      Vector <unsigned int> histogram = featureValues.histogram( 0, 1, histogram.getSize()-1 );

//      thresholding.setThreshold( otsuThresholding.computeThreshold(histogram) );
//      thresholding.apply( nucleiMask );
//      //thresholding.applyAlternative( nucleiMask );
//  //    nucleiMask.save( "/home/jarpon/Desktop/test3.vm", true );
//  }

  nucleiMask.setVoxelCalibration( originalVoxelMatrix.getVoxelCalibration() );


  LEAVE();
  return nucleiMask;
}
