#include <voxelmatrix.h>

#define TRACE
#include <trace.h>
using namespace std;


PixelMatrix<float> getProjection(const VoxelMatrix<float>& originalVoxelMatrix)
{
  PixelMatrix<float> zProj( originalVoxelMatrix.getSize1(), originalVoxelMatrix.getSize2() );
  zProj.setZeros();
  int i, j, k;

  for (i = 0; i < originalVoxelMatrix.getSize1(); ++i)
  for (j = 0; j < originalVoxelMatrix.getSize2(); ++j)
  {
    for (k = 1; k < originalVoxelMatrix.getSize3(); ++k)
    {
      if ( originalVoxelMatrix[k][i][j] > zProj(j,i) )
      {
        zProj(j,i) = originalVoxelMatrix[k][i][j];
      }
    }
  }

//  VoxelMatrix<T>& regionMatrix = *_labelMatrix;
//  const int size1 = regionMatrix.getSize1();
//  const int size2 = regionMatrix.getSize2();
//  const int size3 = regionMatrix.getSize3();

//  PixelMatrix<T> labelProjection( size1, size2 );
//  labelProjection.setZeros();

//  for ( int k = 0; k < size3; ++k )
//    for (int i = 0; i < size1; ++i )
//      for ( int j = 0; j < size2; ++j )
//        if ( regionMatrix(i,j,k) == label )
//          labelProjection(i,j) = 1;

  PixelCalibration pixelCalibration;
  pixelCalibration.setPixelHeight(originalVoxelMatrix.getVoxelCalibration().getVoxelHeight() );
  pixelCalibration.setPixelWidth( originalVoxelMatrix.getVoxelCalibration().getVoxelWidth() );
  pixelCalibration.setLengthUnit( originalVoxelMatrix.getVoxelCalibration().getLengthUnit() );

  zProj.setPixelCalibration( pixelCalibration );

  return zProj;
}

