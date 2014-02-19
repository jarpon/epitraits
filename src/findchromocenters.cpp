#include <iostream>
#include <componentlabelling.h>
#include <otsuthresholding.h>
//#include "regionanalysis2.h"
#include <regionanalysis.h>
#include <thresholding.h>
#include <voxelmatrix.h>
#include <volumehistogramexpansion.h>
//#include <minfilter.h>

#define TRACE
#include <trace.h>

extern VoxelMatrix<float> applyLabelling(const VoxelMatrix<float>&, const VoxelMatrix<float>&, const string&, const string&);

VoxelMatrix <float> findCCs(const VoxelMatrix<float>& originalVoxelMatrix, const VoxelMatrix<float>& nucleusMask,
                            const string& filename, const string& intermediateProcessesDir )
{
  VoxelMatrix<float> regionMatrix;
  regionMatrix = applyLabelling( originalVoxelMatrix, nucleusMask, filename, intermediateProcessesDir );

  VoxelMatrix<float> ccsMask, copyVoxelMatrix = originalVoxelMatrix;

  VolumeHistogramExpansion<float> normalizerHistogram;
  normalizerHistogram.setMinMax( 0.0, 255.0 );
  normalizerHistogram.apply(copyVoxelMatrix);

  RegionAnalysis<float> regionAnalysis;
  regionAnalysis.setRegionMatrix( regionMatrix );
  regionAnalysis.run();
  regionAnalysis.mapRegionFeature( ccsMask, REGION_FEATURE_AVERAGE_VALUE, copyVoxelMatrix );
  ccsMask.setVoxelCalibration( originalVoxelMatrix.getVoxelCalibration() );
  ccsMask.save( intermediateProcessesDir + filename + "-average.vm", true );

  OtsuThresholding<float> otsuThresholding;
  Vector <float> featureValues = regionAnalysis.allRegionValues( ccsMask );
//  Vector <float> featureValues = regionAnalysis.computeRegionFeature( REGION_FEATURE_AVERAGE_VALUE, copyRegionMatrix );
  Vector <unsigned int> histogram = featureValues.histogram( 0, 1, 256 );

/*
  otsuThresholding.setForegroundIsAbove( true );
  otsuThresholding.setBackground( 0.0 );
  otsuThresholding.setForeground( 255.0 );
  otsuThresholding.computeThreshold(histogram);
  otsuThresholding.apply( featureMatrix );
*/

  Thresholding<float> thresholding;
  thresholding.setBackground( 0.0 );
  thresholding.setForeground( 1.0 );
  thresholding.setThreshold( otsuThresholding.computeThreshold(histogram) );
  EVAL(otsuThresholding.computeThreshold(histogram));
  thresholding.apply( ccsMask );

  //labeling the image
  ComponentLabelling<float> componentLabelling;
  componentLabelling.apply( ccsMask );

  //ccsMask.save( chromocentersDir + fileName + ".vm", true );
  return ccsMask;
}

