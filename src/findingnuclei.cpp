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


#define TRACE
#include <trace.h>
using namespace std;

VoxelMatrix<float> findNucleus(const VoxelMatrix<float>& originalVoxelMatrix)
{

  /*
  Vector<int> size = originalVoxelMatrix.getSize();
  for (int i = 0; i < size[0]; ++i)
  {
    for (int j = 0; j < size[1]; ++j)
    {
      for (int k = 0; k < size[2]; ++k)
      {
        if ( originalVoxelMatrix(i,j,k) < 0.4 )
        {
          originalVoxelMatrix(i,j,k) = 0;
        }
      }
    }
  }

*/
  VoxelMatrix<float> nucleusMask = originalVoxelMatrix;
//  VolumeHistogramExpansion<float> normalizerHistogram;
//  normalizerHistogram.apply( nucleusMask );


//  VoxelMatrix<float> structElement;
//  structElement.setSize(7,7,7);
//  structElement.setOnes();

//  VoxelMatrixDilatation<float> voxelDilatation;
//  voxelDilatation.setStructElt( structElement );
//  voxelDilatation.apply( nucleusMask );

//  VoxelMatrix<float> structElement2;
//  structElement2.setSize(3,3,3);
//  structElement2.setOnes();

//  VoxelMatrixErosion<float> voxelErosion;
//  voxelErosion.setStructElt( structElement2 );
//  voxelErosion.apply( nucleusMask );



  OtsuThresholding<float> otsuThresholding;
  otsuThresholding.setForegroundIsAbove( true );
  otsuThresholding.setForeground( 1.0 );
  otsuThresholding.setBackground( 0.0 );
  otsuThresholding.apply( nucleusMask );

  HolesFillingf holesFilling;
  int sizeZ = originalVoxelMatrix.getSize3();

  //to obtain a better filling, it's applied to each 2D slice instead of the complete 3D stack
  for (int k = 0; k < sizeZ; ++k)  holesFilling.apply( nucleusMask[k] );
  holesFilling.apply( nucleusMask );
  //for (int k = 0; k < sizeZ; ++k)  holesFilling.apply( nucleusMask[k] );

  ComponentLabelling<float> componentLabelling;
  componentLabelling.apply( nucleusMask );

  const int numComponents = componentLabelling.getNumLabels();
  const Vector<unsigned int> histogram = nucleusMask.histogram( 0, 1, numComponents+1 );
  const Vector<unsigned int> histogram2 = histogram.copy( 1, histogram.getSize()-1 );

  Thresholding<float> thresholding;
  thresholding.setForegroundIsAbove( true );
//  thresholding.setThreshold( 1.0 );  //for images already segmented
  thresholding.setForeground( 1.0 );
  thresholding.setBackground( 0.0 );
//  thresholding.apply( nucleusMask );  //for images already segmented and comment the next threshold

  //a final threshold is applied in each slide to leave only the biggest shape
  //(wich is labeled with 1, the other labeled shapes are set to 0)
  for (int k = 0; k < sizeZ; ++k)  thresholding.levelSetMask( nucleusMask[k], histogram2.maxPos()+1 );

//  //process to improve nuclei segmentation when the nucleoli touch the envelope
  VoxelMatrix<float> structElement3;
  structElement3.setSize(7,7,7);
  structElement3.setOnes();

  VoxelMatrixDilatation<float> voxelDilatation2;
  voxelDilatation2.setStructElt( structElement3 );
  voxelDilatation2.apply( nucleusMask );

  VoxelMatrix<float> structElement4;
  structElement4.setSize(3,3,3);
  structElement4.setOnes();

  VoxelMatrixErosion<float> voxelErosion2;
  voxelErosion2.setStructElt( structElement4 );
  voxelErosion2.apply( nucleusMask );


  //another 3D filling
  holesFilling.apply( nucleusMask );

  //nucleusMask.save( outputDirectory + fileName + ".vm", true );
  return nucleusMask;
}
