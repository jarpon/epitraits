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
//#include <medianfilter.h>
#include <watershedtransform.h>
#include <gaussiangradient.h>

#include <cmath>


#define TRACE
#include <trace.h>

//extern VoxelMatrix<float> applyLabelling(const VoxelMatrix<float>&, const VoxelMatrix<float>&, const string&, const string&);

VoxelMatrix <float> findChromosomes(const VoxelMatrix<float>& originalVoxelMatrix, VoxelMatrix<float>& nucleusMask,
                            const string& filename, const string& intermediateProcessesDir )
{
  //VoxelMatrix<float> regionMatrix;
  //regionMatrix = applyLabelling( originalVoxelMatrix, nucleusMask, filename, intermediateProcessesDir );

  VoxelMatrix<float> regionMatrix = originalVoxelMatrix, rangeMask = originalVoxelMatrix, copyVoxelMatrix = originalVoxelMatrix;

  GaussianGradient<float> gaussianGradient;
  gaussianGradient.MaskVoxelMatrixProcessing<float>::setMask( nucleusMask );
  gaussianGradient.setSigma( 1.4 );
  gaussianGradient.apply( regionMatrix );
  regionMatrix.save( intermediateProcessesDir + filename + "-chromo-gradient.vm", true );

  WatershedTransform<float> watershedTransform;
  watershedTransform.MaskVoxelMatrixProcessing<float>::setMask( nucleusMask );
  watershedTransform.apply( regionMatrix );

  RegionAnalysis<float> regionAnalysis;
  regionAnalysis.setRegionMatrix( regionMatrix );
  regionAnalysis.run();
  regionAnalysis.mapRegionFeature( rangeMask, REGION_FEATURE_CONTRAST, copyVoxelMatrix );
  //regionAnalysis.mapRegionFeature( rangeMask, REGION_FEATURE_MAXIMUM_VALUE, copyVoxelMatrix );
  rangeMask.save( intermediateProcessesDir + filename + "-chromo-max.vm", true );

  VoxelMatrix<float> chromosomeMask = rangeMask;

  OtsuThresholding<float> otsuThresholding;
  //Vector <float> featureValues = regionAnalysis.allRegionValues( rangeMask );
  Vector <float> featureValues = regionAnalysis.computeRegionFeature( REGION_FEATURE_CONTRAST, copyVoxelMatrix );
  //Vector <float> featureValues = regionAnalysis.computeRegionFeature( REGION_FEATURE_MAXIMUM_VALUE, copyVoxelMatrix );
  //Vector <unsigned int> histogram = featureValues.histogram( 0, 1, 255 );
  //featureValues.sort();
  Vector <unsigned int> histogram = featureValues.histogram( featureValues.min(), 1, floor(featureValues.max())+1  );

  //test contrast*shape
  //Vector <float> shapes = regionAnalysis.computeRegionFeature( REGION_FEATURE_SPHERICITY, copyVoxelMatrix );
  //Vector <float> result = featureValues * shapes;
  //result.sort();
  //Vector <unsigned int> histogram = result.histogram( 0, 1, result.getSize()-1 );

  otsuThresholding.setForegroundIsAbove( true );
  otsuThresholding.setBackground( 0.0 );
  otsuThresholding.setForeground( 1.0 );
  double factor = otsuThresholding.computeThreshold(histogram);
  otsuThresholding.setThreshold( factor );
  EVAL(factor);
  EVAL( otsuThresholding.getThreshold() );
  otsuThresholding.applyAlternative( chromosomeMask );

  ComponentLabelling<float> componentLabelling;
  componentLabelling.apply( chromosomeMask );
  EVAL(componentLabelling.getNumLabels());
  EVAL(nucleusMask.sum().sum().sum());

  float threshold = factor;
  bool finish = false;
  do {
    EVAL(chromosomeMask.sum().sum().sum());
    //EVAL(rangeMask.max().max().max());

    if ( chromosomeMask.sum().sum().sum() > nucleusMask.sum().sum().sum()*0.2)
    {
      chromosomeMask = rangeMask;
      threshold *= 1.2;
      EVAL(threshold);
      otsuThresholding.setThreshold( threshold );
      otsuThresholding.applyAlternative( chromosomeMask );
    }
    if ( threshold > rangeMask.max().max().max()*0.9 )
    {
//      WatershedTransform<float> watershedTransform2;
//      watershedTransform2.MaskVoxelMatrixProcessing<float>::setMask( chromosomeMask );
//      watershedTransform2.apply( regionMatrix );

//      RegionAnalysis<float> regionAnalysis2;
//      regionAnalysis2.setRegionMatrix( regionMatrix );
//      regionAnalysis2.run();
//      regionAnalysis2.mapRegionFeature( rangeMask, REGION_FEATURE_MAXIMUM_VALUE, copyVoxelMatrix );
//      rangeMask.save( intermediateProcessesDir + filename + "-chromo-max.vm", true );

//      chromosomeMask = rangeMask;

//      Vector <float> featureValues2 = regionAnalysis.computeRegionFeature( REGION_FEATURE_MAXIMUM_VALUE, copyVoxelMatrix );
//      featureValues2.sort();
//      Vector <unsigned int> histogram2 = featureValues2.histogram( 0, 1, featureValues2.getSize()-1 );

//      otsuThresholding.setThreshold( otsuThresholding.computeThreshold(histogram2) );
//      EVAL( otsuThresholding.getThreshold() );
//      otsuThresholding.apply( chromosomeMask );
      finish = true;
    }
    componentLabelling.apply( chromosomeMask );
    if ( componentLabelling.getNumLabels() > 3 )
    {
      Vector<unsigned int> histogram = chromosomeMask.histogram( 0, 1, componentLabelling.getNumLabels()+1 );
      //Vector<unsigned int> histogram2 = histogram.copy( 1, histogram.getSize()-1 );
     // histogram2.reverse();
      EVAL(histogram);
      //EVAL(histogram2);

      bool twoChromosomes = false;
      int posChromosome = histogram.maxPos()+1;
      int chromosomeArea = histogram[posChromosome];
      EVAL(chromosomeArea);
      int pos2Chromosome;

      for ( int pos = 0; pos < histogram.getSize(); ++pos )
      {
        if ( histogram[pos] > chromosomeArea/2 )
        {
          pos2Chromosome = pos;
          twoChromosomes = true;
          regionMatrix = chromosomeMask;
        }
      }

      int sizeZ = chromosomeMask.getSize3();
      Thresholding<float> thresholding;
      thresholding.setBackground( 0.0 );
      //thresholding.setThreshold( histogram2.maxPos()+1 );
      //thresholding.applyAlternative( chromosomeMask );
      chromosomeMask.save( intermediateProcessesDir + filename + "-chromo-1.vm", true );

      for (int k = 0; k < sizeZ; ++k)  thresholding.levelSetMask( chromosomeMask[k], posChromosome );

      chromosomeMask.save( intermediateProcessesDir + filename + "-chromo-2.vm", true );

      if ( twoChromosomes == true )
      {
        Thresholding<float> secondThresholding;
        secondThresholding.setBackground( 0.0 );
        //secondThresholding.setThreshold( posChromosome );
        //secondThresholding.applyAlternative( regionMatrix );
        for (int k = 0; k < sizeZ; ++k)  secondThresholding.levelSetMask( regionMatrix[k], pos2Chromosome );

        chromosomeMask.operator+=(regionMatrix);

        finish = true;
      }
    }
//    else
//    {
//      EVAL( chromosomeMask.sum().sum().sum()/nucleusMask.sum().sum().sum());
//      VoxelMatrix<float> chromosomeMask = rangeMask;
//      threshold *= 1.2;
//      EVAL(threshold);
//      otsuThresholding.setThreshold( threshold );
//      otsuThresholding.applyAlternative( chromosomeMask );
//    }
  } while ( ( componentLabelling.getNumLabels() > 2 ) && ( finish == false ) && ( componentLabelling.getNumLabels() > 0 || chromosomeMask.sum().sum().sum()/nucleusMask.sum().sum().sum() > 0.2) );
  VoxelMatrix<float> structElement;
  structElement.setSize(5,5,5);
  structElement.setOnes();

  VoxelMatrixDilatation<float> voxelDilatation;
  voxelDilatation.setStructElt( structElement );
  voxelDilatation.apply( chromosomeMask );

  VoxelMatrix<float> structElement2;
  structElement2.setSize(5,5,5);
  structElement2.setOnes();

  VoxelMatrixErosion<float> voxelErosion;
  voxelErosion.setStructElt( structElement2 );
  voxelErosion.apply( chromosomeMask );

  componentLabelling.apply( chromosomeMask );
  chromosomeMask.setVoxelCalibration( originalVoxelMatrix.getVoxelCalibration() );

  return chromosomeMask;
}
