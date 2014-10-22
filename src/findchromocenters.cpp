#include <iostream>
#include <componentlabelling.h>
#include <otsuthresholding.h>
//#include "regionanalysis2.h"
#include <regionanalysis.h>
#include <thresholding.h>
#include <voxelmatrix.h>
#include <volumehistogramexpansion.h>
#include <voxelmatrixdilatation.h>
#include <voxelmatrixerosion.h>
#include <holesfilling.h>
#include <medianfilter.h>

#include <cmath>


#define TRACE
#include <trace.h>

extern VoxelMatrix<float> applyLabelling(const VoxelMatrix<float>&, const VoxelMatrix<float>&, const string&, const string&);

VoxelMatrix <float> findCCs(const VoxelMatrix<float>& originalVoxelMatrix, VoxelMatrix<float>& nucleusMask,
                            const string& filename, const string& intermediateProcessesDir )
{
//  VoxelMatrix<float> regionMatrix;
//  regionMatrix = applyLabelling( originalVoxelMatrix, nucleusMask, filename, intermediateProcessesDir );

  VoxelMatrix<float> regionMatrix( intermediateProcessesDir + filename + "-gradient.vm" );
  EVAL(regionMatrix.getSize());

  VoxelMatrix<float> rangeMask, copyVoxelMatrix = originalVoxelMatrix;


  RegionAnalysis<float> regionAnalysis;
  regionAnalysis.setRegionMatrix( regionMatrix );
  regionAnalysis.run();

  EVAL("!");
//  regionAnalysis.mapRegionFeature( rangeMask, REGION_FEATURE_AVERAGE_VALUE, copyVoxelMatrix );
//  rangeMask.save( intermediateProcessesDir + filename + "-average.vm", true );
  regionAnalysis.mapRegionFeature( rangeMask, REGION_FEATURE_CONTRAST, copyVoxelMatrix );
  rangeMask.save( intermediateProcessesDir + filename + "-contrast.vm", true );
//  regionAnalysis.mapRegionFeature( rangeMask, REGION_FEATURE_MAXIMUM_VALUE, copyVoxelMatrix );
//  rangeMask.save( intermediateProcessesDir + filename + "-max.vm", true );


  VoxelMatrix<float> ccsMask = rangeMask;

  OtsuThresholding<float> otsuThresholding;
  Vector <float> featureValues = regionAnalysis.computeRegionFeature( REGION_FEATURE_CONTRAST, copyVoxelMatrix );
//  Vector <float> featureValues = regionAnalysis.computeRegionFeature( REGION_FEATURE_CONTRAST, rangeMask );
  //Vector <float> featureValues = regionAnalysis.computeRegionFeature( REGION_FEATURE_AVERAGE_VALUE, rangeMask );

  //Vector <unsigned int> histogram = featureValues.histogram( 0, 1, featureValues.getSize()-1 );
  //Vector <unsigned int> histogram = featureValues.histogram( 0, 1, featureValues.max() );
  Vector <unsigned int> histogram = featureValues.histogram( featureValues.min(), 1, floor(featureValues.max())+1 );

//  featureValues.sort();
//  EVAL(featureValues);

//  otsuThresholding.setForegroundIsAbove( true );
//  otsuThresholding.setBackground( 0.0 );
//  otsuThresholding.setForeground( 1.0 );
//  otsuThresholding.computeThreshold(histogram);
//  otsuThresholding.applyAlternative( ccsMask );


  Thresholding<float> thresholding;
  thresholding.setBackground( 0.0 );
  //thresholding.setForeground( 1.0 );


//  float factor ;
//  Vector<float> reverseFeatureValues(featureValues.getSize());
//  reverseFeatureValues = featureValues.diff(1);
//  //reverseFeatureValues.sort();
//  reverseFeatureValues.reverse();
//  EVAL(reverseFeatureValues);



  /* Adaptative threshold:
   * (1) IMPLEMENTED: this works with the differences between pairs of the histogram obtained from the contrast map.
   * What this looks for is a stabilization of the constrast values (chromocenters are above (~40 - ~100) of the nuclei values (from negative values to ~40))
   * and after a fall from the brightest contrast through the different chromocenters usually there is this stabilization.
   * To take into account it should be a row of close values that are similar (difference < 2)
   * And it looks for a stabilization found by 6 neighbour values which a difference below 8 between the poles.
   *
   * It works fine but not in all images!!
   */

  /* (2) Other alternative would be to look for a stabilization after a large fall (this works in a lower porcentage.
   */

//  int threshold, counter = 0;

//  for (int i = 0; i < reverseFeatureValues.getSize(); ++i)
//  {
//      if ( (reverseFeatureValues[i] < 2) && ( (i - threshold) == 1 || counter == 0) && featureValues[featureValues.getSize()-i] < featureValues.max()*0.6 )
//      {
//        threshold = i;
//        ++counter;
//      }
//      else counter = 0;

//      if ( counter > 5 && (featureValues[featureValues.getSize()-threshold+counter]-featureValues[featureValues.getSize()-threshold]<8))
//      {
//          EVAL(featureValues[featureValues.getSize()-threshold+counter]);
//          EVAL(featureValues[featureValues.getSize()-threshold]);
//          EVAL(threshold);
//          if ( threshold > 13 )
//            i = reverseFeatureValues.getSize();
//          else
//          {
//              i = threshold-counter+1;
//              threshold = 0;
//              counter = 0;
//          }
//      }
//      else if ( counter > 5 && (featureValues[featureValues.getSize()-threshold+counter]-featureValues[featureValues.getSize()-threshold]>8 ))
//      {
//          EVAL(featureValues[featureValues.getSize()-threshold+counter]);
//          EVAL(featureValues[featureValues.getSize()-threshold]);
//          EVAL(threshold);
//          i = threshold-counter+1;
//          threshold = 0;
//          counter = 0;
//      }
//  }
//  EVAL(threshold)
//  EVAL(counter)

//  EVAL(featureValues[featureValues.getSize() - (threshold-counter) +1 ] + 0.01)

//  if ( counter != 0 ) factor = featureValues[featureValues.getSize() - (threshold-counter) +1 ] - 0.01 ;
//  else factor = otsuThresholding.computeThreshold(histogram);

//  EVAL(otsuThresholding.computeThreshold(histogram));
  //factor = otsuThresholding.computeThreshold(histogram) + abs(featureValues.mean());
  //factor = otsuThresholding.computeThreshold(histogram) * 1.5;
  //factor = featureValues.max() / 1.5;

  //thresholding.setThreshold( factor );

  thresholding.setThreshold( otsuThresholding.computeThreshold(histogram) );
  //thresholding.setThreshold( otsuThresholding.computeThreshold(featureValues) );


//  EVAL(otsuThresholding.computeThreshold(histogram));
//  EVAL(featureValues.max())

  thresholding.applyAlternative( ccsMask );
  //thresholding.apply( ccsMask );
  //ccsMask.save( "/home/jarpon/Desktop/test1.vm", true );

  //to obtain a better filling, it's applied to each 2D slice instead of the complete 3D stack
//  HolesFillingf holesFilling;
//  int sizeZ = ccsMask.getSize3();
//  for (int k = 0; k < sizeZ; ++k)  holesFilling.apply( ccsMask[k] );


//  VoxelMatrix<float> ccsFillMask = ccsMask;
//  ccsMask.fillIt(ccsFillMask);


  VoxelMatrix<float> structElement;
  structElement.setSize(3,3,3);
  structElement.setOnes();

  VoxelMatrixDilatation<float> voxelDilatation;
  voxelDilatation.setStructElt( structElement );
  voxelDilatation.apply( ccsMask );

 //  VoxelMatrix<float> structElement2;
 //  structElement2.setSize(1,1,1);
 //  structElement2.setOnes();

  VoxelMatrixErosion<float> voxelErosion;
  voxelErosion.setStructElt( structElement );
  voxelErosion.apply( ccsMask );


  //labeling the image
  ComponentLabelling<float> componentLabelling;
  componentLabelling.apply( ccsMask );

  ccsMask.setVoxelCalibration( originalVoxelMatrix.getVoxelCalibration() );

  //ccsMask.save( chromocentersDir + filename + ".vm", true );

  return ccsMask;
}
