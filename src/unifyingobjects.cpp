#include <voxelmatrix.h>
#include <fileinfo.h>
//#include <regionanalysis.h>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <cmath>

#define TRACE
#include <trace.h>
using namespace std;


VoxelMatrix<float> unifyLabels( const VoxelMatrix<float>& originalVM, const int& oldLabel, const int& newLabel )
{
    const int size1 = originalVM.getSize1();
    const int size2 = originalVM.getSize2();
    const int size3 = originalVM.getSize3();
    int i, j, k;

    VoxelMatrix<float> changedVM = originalVM;

    for (i = 0; i < size1; i++)
      for (j = 0; j < size2; j++)
        for (k = 0; k < size3; k++)
          if ( originalVM(i,j,k) == oldLabel )
            changedVM(i,j,k) = newLabel;

    changedVM.setVoxelCalibration( originalVM.getVoxelCalibration() );

    return changedVM;
}
