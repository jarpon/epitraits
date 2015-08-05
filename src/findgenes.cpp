#include <iostream>
#include <componentlabelling.h>
//#include "regionanalysis2.h"
#include <regionanalysis3d.h>
#include "thresholding.h"
#include "voxelmatrix.h"
#include <volumehistogramexpansion.h>
#include <voxelmatrixdilatation.h>
#include <voxelmatrixerosion.h>
#include <holesfilling.h>
#include <medianfilter.h>
#include "otsuthresholding.h"

#include <cmath>


#define TRACE
#include <trace.h>

extern VoxelMatrix<float> applyLabelling(const VoxelMatrix<float>&, const VoxelMatrix<float>&, const string&, const string&);

VoxelMatrix <float> findGenes(const VoxelMatrix<float>& originalVoxelMatrix, VoxelMatrix<float>& nucleusMask,
                            const string& filename, const string& intermediateProcessesDir )
{
  VoxelMatrix<float> regionMatrix = originalVoxelMatrix;
//  nucleusMask.setOnes();



//  OtsuThresholding<float> otsuThresholding2;
//  otsuThresholding2.setForegroundIsAbove( true );
//  otsuThresholding2.setForeground( 1.0 );
//  otsuThresholding2.setBackground( 0.0 );
//  otsuThresholding2.apply( regionMatrix );
//  regionMatrix.save( intermediateProcessesDir + filename + "-a2.vm", true );




  regionMatrix = applyLabelling( originalVoxelMatrix, nucleusMask, filename, intermediateProcessesDir );
  //const int sizeZ = originalVoxelMatrix.getSize3();

  VoxelMatrix<float> rangeMask, copyVoxelMatrix = originalVoxelMatrix;

//  VolumeHistogramExpansion<float> normalizerHistogram;
//  normalizerHistogram.setMinMax( 0.0, 255.0 );
//  normalizerHistogram.apply(copyVoxelMatrix);

//  RegionAnalysis3D<float> regionAnalysis;
//  regionAnalysis.setLabelMatrix( regionMatrix );
//  regionAnalysis.setValueMatrix( copyVoxelMatrix );
//  regionAnalysis.setOutputMatrix( rangeMask );
//  regionAnalysis.run();

//  rangeMask.setSize( copyVoxelMatrix.getSize() );
//  rangeMask.setZeros();


//  regionAnalysis.outputFillRegions( REGION_FEATURE_CONTRACTNESS );
//  rangeMask.save( intermediateProcessesDir + filename + "-contrast.vm", true );

  VoxelMatrix<float> genesMask = rangeMask;

////  int regionsNumber = regionAnalysis.numRegions();
////  Vector<float> regionValues = regionAnalysis.allRegionValues( genesMask );
////  EVAL(regionsNumber);

//  OtsuThresholding<float> otsuThresholding;
//  //Vector <float> featureValues = regionAnalysis.allRegionValues( genesMask );

////  Vector <float> featureValues = regionAnalysis.computeRegionFeature( REGION_FEATURE_CONTRACTNESS );
////  Vector <unsigned int> histogram = featureValues.histogram( 0, 1, histogram.getSize()-1 );

//  //featureValues.sort();
//  EVAL(featureValues);


//  Thresholding<float> thresholding;
//  thresholding.setBackground( 0.0 );
//  //thresholding.setForeground( 1.0 );
//  //thresholding.setThreshold(160);

//  float factor ;
//  Vector<float> reverseFeatureValues(featureValues.getSize());
//  reverseFeatureValues = featureValues.diff(1);
//  reverseFeatureValues.reverse();
//  EVAL(reverseFeatureValues);



//  /* Adaptative threshold:
//   * (1) IMPLEMENTED: this works with the differences between pairs of the histogram obtained from the contrast map.
//   * What this looks for is a stabilization of the constrast values (chromocenters are above (~40 - ~100) of the nuclei values (from negative values to ~40))
//   * and after a fall from the brightest contrast through the different chromocenters usually there is this stabilization.
//   * To take into account it should be a row of close values that are similar (difference < 2)
//   * And it looks for a stabilization found by 6 neighbour values which a difference below 8 between the poles.
//   *
//   * It works fine but not in all images!!
//   */

//  /* (2) Other alternative would be to look for a stabilization after a large fall (this works in a lower porcentage.
//   */

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
////          EVAL(featureValues[featureValues.getSize()-threshold+counter]);
////          EVAL(featureValues[featureValues.getSize()-threshold]);
////          EVAL(threshold);
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
////          EVAL(featureValues[featureValues.getSize()-threshold+counter]);
////          EVAL(featureValues[featureValues.getSize()-threshold]);
////          EVAL(threshold);
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

////  EVAL(otsuThresholding.computeThreshold(histogram));
//  //factor = otsuThresholding.computeThreshold(histogram) + abs(featureValues.mean());
//  //factor = otsuThresholding.computeThreshold(histogram) * 1.5;
//  //factor = featureValues.max() / 1.5;

//  thresholding.setThreshold( factor );
//  //thresholding.setThreshold( 7 );
//  //thresholding.setThreshold( otsuThresholding.computeThreshold(histogram) );


////  EVAL(otsuThresholding.computeThreshold(histogram));
////  EVAL(featureValues.max())

//  thresholding.applyAlternative( genesMask );
//  //thresholding.apply( genesMask );
//  //genesMask.save( "/home/jarpon/Desktop/test1.vm", true );

////  HolesFillingf holesFilling;
////  int sizeZ = originalVoxelMatrix.getSize3();

////  //to obtain a better filling, it's applied to each 2D slice instead of the complete 3D stack
////  for (int k = 0; k < sizeZ; ++k)  holesFilling.apply( genesMask[k] );
////  holesFilling.apply( genesMask );

//  //for (int k = 0; k < sizeZ; ++k)  holesFilling.apply( nucleusMask[k] );
//  //labeling the image

////  HolesFillingf holesFilling;
////  holesFilling.apply( genesMask );

//  ComponentLabelling<float> componentLabelling;
//  componentLabelling.apply( genesMask );

//  VoxelMatrix<float> structElement;
//  structElement.setSize(3,3,3);
//  structElement.setOnes();

//  VoxelMatrixDilatation<float> voxelDilatation;
//  voxelDilatation.setStructElt( structElement );
//  voxelDilatation.apply( genesMask );

////  VoxelMatrix<float> structElement2;
////  structElement2.setSize(3,3,1);
////  structElement2.setOnes();

//  VoxelMatrixErosion<float> voxelErosion;
//  voxelErosion.setStructElt( structElement );
//  voxelErosion.apply( genesMask );

//  //holesFilling.apply( genesMask );

//  genesMask.setVoxelCalibration( originalVoxelMatrix.getVoxelCalibration() );

//  //genesMask.save( chromocentersDir + filename + ".vm", true );

  return genesMask;
}

