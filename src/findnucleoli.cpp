#include <iostream>
#include <componentlabelling.h>
#include <otsuthresholding.h>
//#include "regionanalysis2.h"
#include <regionanalysis3d.h>
#include "thresholding.h"
#include "voxelmatrix.h"
#include <volumehistogramexpansion.h>
#include <voxelmatrixdilatation.h>
#include <voxelmatrixerosion.h>
#include <holesfilling.h>
#include <medianfilter.h>
#include <gaussiangradient.h>
#include <watershedtransform.h>

#include <cmath>


#define TRACE
#include <trace.h>

VoxelMatrix <float> findNucleoli(const VoxelMatrix<float>& originalVoxelMatrix, VoxelMatrix<float>& nucleusMask,
                            const string& filename, const string& intermediateProcessesDir )
{
  EVAL("1");

  VoxelMatrix<float> nucleusMaskCopy = nucleusMask;
  VoxelMatrix<float> gradientMatrix = originalVoxelMatrix;

  EVAL("1");
  //labelling
  GaussianGradient<float> gaussianGradient;
  gaussianGradient.MaskVoxelMatrixProcessing<float>::setMask( nucleusMaskCopy );
  gaussianGradient.setSigma( 2 );
  gaussianGradient.apply( gradientMatrix );
 // gradientMatrix.save( intermediateProcessesDir + filename + "-gradient.vm", true );

  VoxelMatrix <float> regionMatrix = gradientMatrix;
  WatershedTransform<float> watershedTransform;
  watershedTransform.MaskVoxelMatrixProcessing<float>::setMask( nucleusMaskCopy );
  watershedTransform.apply( regionMatrix );
  regionMatrix.save( intermediateProcessesDir + filename + "-watershed.vm", true );

  VoxelMatrix<float> rangeMask, copyVoxelMatrix = originalVoxelMatrix;
  VoxelMatrix<float> nucleoliMask;

  /*

  RegionAnalysis3D<float> regionAnalysis;
  regionAnalysis.setLabelMatrix( regionMatrix );
  regionAnalysis.setValueMatrix( copyVoxelMatrix );
  regionAnalysis.setOutputMatrix( rangeMask );
  regionAnalysis.run();

  regionAnalysis.outputFillRegions( REGION_FEATURE_CONTRACTNESS );
  rangeMask.save( intermediateProcessesDir + filename + "-contrast.vm", true );

  VoxelMatrix<float> nucleoliMask = rangeMask;
  OtsuThresholding<float> otsuThresholding;

  Vector <float> featureValues = regionAnalysis.computeRegionFeature( REGION_FEATURE_CONTRACTNESS );
  Vector <unsigned int> histogram = featureValues.histogram( featureValues.min(), 1, floor(featureValues.max())+1 );

  Thresholding<float> thresholding;
  thresholding.setBackground( 0.0 );
  thresholding.setForeground( 1.0 );
  EVAL("1");
  thresholding.setThreshold( otsuThresholding.computeThreshold(histogram) );

  //thresholding.applyAlternative( nucleoliMask );
  thresholding.apply( nucleoliMask );

  //to obtain a better filling, it's applied to each 2D slice instead of the complete 3D stack
  HolesFillingf holesFilling;
  int sizeZ = nucleoliMask.getSize3();
  for (int k = 0; k < sizeZ; ++k)  holesFilling.apply( nucleoliMask[k] );


//  VoxelMatrix<float> ccsFillMask = nucleoliMask;
//  nucleoliMask.fillIt(ccsFillMask);


  VoxelMatrix<float> structElement;
  structElement.setSize(1,1,1);
  structElement.setOnes();

  VoxelMatrixDilatation<float> voxelDilatation;
  voxelDilatation.setStructElt( structElement );
  voxelDilatation.apply( nucleoliMask );

  VoxelMatrixErosion<float> voxelErosion;
  voxelErosion.setStructElt( structElement );
  voxelErosion.apply( nucleoliMask );

//  VoxelMatrix<float> ccsFillMask = nucleoliMask;
//  nucleoliMask.fillIt(ccsFillMask);

  holesFilling.apply( nucleoliMask );
  for (int k = 0; k < sizeZ; ++k)  holesFilling.apply( nucleoliMask[k] );

  //labeling the image
  ComponentLabelling<float> componentLabelling;
  componentLabelling.apply( nucleoliMask );

  nucleoliMask.setVoxelCalibration( originalVoxelMatrix.getVoxelCalibration() );

  //nucleoliMask.save( chromocentersDir + filename + ".vm", true );
*/
  return nucleoliMask;
}

