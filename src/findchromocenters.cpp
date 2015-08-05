#include <iostream>
#include <componentlabelling.h>
#include "otsuthresholding.h"
//#include "regionanalysis2.h"
//#include "regionanalysis.h"
#include <regionanalysis3d.h>
#include "thresholding.h"
#include "voxelmatrix.h"
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
  VoxelMatrix<float> regionMatrix;
  regionMatrix = applyLabelling( originalVoxelMatrix, nucleusMask, filename, intermediateProcessesDir );

//  VoxelMatrix<float> regionMatrix( intermediateProcessesDir + filename + "-gradient.vm" );
//  EVAL(regionMatrix.getSize());

  VoxelMatrix<float> rangeMask, copyVoxelMatrix = originalVoxelMatrix;


  RegionAnalysis3D<float> regionAnalysis;
  regionAnalysis.setLabelMatrix( regionMatrix );
  regionAnalysis.setValueMatrix( copyVoxelMatrix );
  regionAnalysis.run();

  rangeMask.setSize( copyVoxelMatrix.getSize() );
  rangeMask.setZeros();
  regionAnalysis.setOutputMatrix( rangeMask );

  regionAnalysis.outputFillRegions( REGION_FEATURE_CONTRAST );
  //regionAnalysis.outputFillRegions( REGION_FEATURE_CONTRACTNESS );

  //rangeMask.save( intermediateProcessesDir + filename + "-contrast.vm", true );
  rangeMask.save( intermediateProcessesDir + filename + "-contractness.vm", true );

  OtsuThresholding<float> otsuThresholding;

  Vector <float> featureValues = regionAnalysis.computeRegionFeature( REGION_FEATURE_CONTRAST );
  //Vector <float> featureValues = regionAnalysis.computeRegionFeature( REGION_FEATURE_CONTRACTNESS );
  featureValues.sort();
  Vector <unsigned int> histogram = featureValues.histogram( featureValues.min(), 1, floor(featureValues.max())+1 );

  float threshold = otsuThresholding.computeThreshold( histogram );
  EVAL(featureValues);
  EVAL(histogram);
  EVAL(threshold);
  regionAnalysis.thresholdRegions( featureValues, threshold );
  regionAnalysis.run();
  int num = regionAnalysis.condenseRegionLabels();
  regionAnalysis.run();

  EVAL(num);
  VoxelMatrix<float> ccsMask = regionAnalysis.getLabelMatrix();

  //to obtain a better filling, it's applied to each 2D slice instead of the complete 3D stack
  HolesFillingf holesFilling;
  int sizeZ = ccsMask.getSize3();
  for (int k = 0; k < sizeZ; ++k)  holesFilling.apply( ccsMask[k] );


//  VoxelMatrix<float> ccsFillMask = ccsMask;
//  ccsMask.fillIt(ccsFillMask);


//  VoxelMatrix<float> structElement;
//  structElement.setSize(3,3,3);
//  structElement.setOnes();

//  VoxelMatrixDilatation<float> voxelDilatation;
//  voxelDilatation.setStructElt( structElement );
//  voxelDilatation.apply( ccsMask );

// //  VoxelMatrix<float> structElement2;
// //  structElement2.setSize(1,1,1);
// //  structElement2.setOnes();

//  VoxelMatrixErosion<float> voxelErosion;
//  voxelErosion.setStructElt( structElement );
//  voxelErosion.apply( ccsMask );

//  for (int k = 0; k < sizeZ; ++k)  holesFilling.apply( ccsMask[k] );

  ccsMask.setVoxelCalibration( originalVoxelMatrix.getVoxelCalibration() );

  //ccsMask.save( chromocentersDir + filename + ".vm", true );

  return ccsMask;
}
