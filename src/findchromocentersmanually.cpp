#include <iostream>
#include <componentlabelling.h>
#include <regionanalysis.h>
#include <thresholding.h>
#include <voxelmatrix.h>
#include <voxelmatrixdilatation.h>
#include <voxelmatrixerosion.h>
#include <holesfilling.h>

#include <cmath>
#include <string>
#include <sstream>

using namespace std;

#define TRACE
#include <trace.h>

extern VoxelMatrix<float> applyLabelling(const VoxelMatrix<float>&, const VoxelMatrix<float>&, const string&, const string&);

void findCCsManually(const VoxelMatrix<float>& originalVoxelMatrix, VoxelMatrix<float>& nucleusMask,
                            const string& filename, const string& intermediateProcessesDir , const string& chromocentersDir)
{
  VoxelMatrix<float> regionMatrix;
  regionMatrix = applyLabelling( originalVoxelMatrix, nucleusMask, filename, intermediateProcessesDir );

  VoxelMatrix<float> rangeMask, copyVoxelMatrix = originalVoxelMatrix;

  RegionAnalysis<float> regionAnalysis;
  regionAnalysis.setRegionMatrix( regionMatrix );
  regionAnalysis.run();
  regionAnalysis.mapRegionFeature( rangeMask, REGION_FEATURE_CONTRAST, copyVoxelMatrix );

  VoxelMatrix<float> ccsMask = rangeMask;

  float inputThreshold = 0;
  string input = "";
  while (true) {
   cout << "Please enter the required value for the correct segmentation:\n>";
   getline(cin, input);

   // This code converts from string to number safely.
   stringstream myStream(input);
   if (myStream >> inputThreshold)
     break;
   cout << "Invalid number, please try again" << endl;
  }
  cout << "You are going to use: " << inputThreshold << " as threshold" << endl << endl;

  EVAL(inputThreshold);

  if ( inputThreshold == -1 )
    return;
  else
  {

    Thresholding<float> thresholding;
    thresholding.setBackground( 0.0 );
    thresholding.setThreshold( inputThreshold );
    thresholding.applyAlternative( ccsMask );

  //  HolesFillingf holesFilling;
  //  int sizeZ = originalVoxelMatrix.getSize3();

    //to obtain a better filling, it's applied to each 2D slice instead of the complete 3D stack
  //  for (int k = 0; k < sizeZ; ++k)  holesFilling.apply( ccsMask[k] );
  //  holesFilling.apply( ccsMask );

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

    ccsMask.save ( chromocentersDir + filename + ".vm", true );

  }

  //return ccsMask;
  return;
}
