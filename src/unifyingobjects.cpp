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


//VoxelMatrix<float> unifyLabels( const VoxelMatrix<float>& originalVM, const int& oldLabel, const int& newLabel )
VoxelMatrix<float> unifyLabels( const VoxelMatrix<float>& originalVM )
{
    const int size1 = originalVM.getSize1();
    const int size2 = originalVM.getSize2();
    const int size3 = originalVM.getSize3();
    int i, j, k;

    int oldLabel;
    string input1 = "";
    while (true) {
     cout << "Please enter the label you want to replace:\n>";
     getline(cin, input1);

     // This code converts from string to number safely.
     stringstream myStream(input1);
     if (myStream >> oldLabel)
       break;
     cout << "Invalid number, please try again" << endl;
    }

    int newLabel;
    string input2 = "";
    while (true) {
     cout << "Please enter the new label you want to have:\n>";
     getline(cin, input2);

     // This code converts from string to number safely.
     stringstream myStream(input2);
     if (myStream >> newLabel)
       break;
     cout << "Invalid number, please try again" << endl;
    }


    cout << "You are going to replace: " << oldLabel << " with " << newLabel << endl << endl;


    VoxelMatrix<float> changedVM = originalVM;

    for (i = 0; i < size1; i++)
      for (j = 0; j < size2; j++)
        for (k = 0; k < size3; k++)
          if ( changedVM(i,j,k) == oldLabel )
            changedVM(i,j,k) = newLabel;

    changedVM.setVoxelCalibration( originalVM.getVoxelCalibration() );

    return changedVM;
}
