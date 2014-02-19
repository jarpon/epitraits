#include <iostream>
#include <dataset.h>
#include <fileinfo.h>
//#include <math.h>
//#include "regionanalysis2.h"
#include <regionanalysis.h>
#include <voxelmatrix.h>
#include <componentlabelling.h>
#include <watershedtransform.h>


#define TRACE
#include <trace.h>
using namespace std;

extern VoxelMatrix<float> findNucleus(const VoxelMatrix<float>&);
extern void nucleusAnalysis(const VoxelMatrix<float>&, VoxelMatrix<float>&, const string&, const int&, DataSet&);
extern VoxelMatrix<float> findCCs( const VoxelMatrix<float>&, const VoxelMatrix<float>&,
                                   const string&, const string&);
extern void chromocentersAnalysis(const VoxelMatrix<float>&, VoxelMatrix<float>&, VoxelMatrix<float>&,
                                  const string&, const int&, int&, DataSet&, DataSet&);

void analyseSample(
    const string& filepath, const int& numNucleus, int& totalNumCCs,
    DataSet& nucleiDataset, DataSet& chromocentersDataset)
{


/*//////////////////////////            OPEN FILE            ///////////////////////////*/

  FileInfo fileInfo( filepath );
  const string casename = fileInfo.baseName();
  const string inputDir = fileInfo.dirName();
  const string outputDirSegmentation = inputDir + "nuclei/";


  //the image is initialized (optional to open directly the nucleus from its own file)
  VoxelMatrix<float> originalVoxelMatrix( filepath );//open an original vm file
//  VoxelMatrix<float> nucleusMask( outputDir + "/nucleus/" + casename + ".vm" ); //open the nucleus segmented file, if it is called, remember to comment the segmentation below

  /*! Nucleus segmentation.
  ****************************************************************/
  VoxelMatrix<float> nucleusMask;
  //nucleusMask = findNucleus( filepath, casename, outputDirSegmentation );//get the nucleus mask

  /*! Nucleus analysis.
  ****************************************************************/
  //nucleusAnalysis( casename, numNucleus, outputDirSegmentation );//get the nucleus mask

  /*! Nucleus analysis.
  ****************************************************************/
  VoxelMatrix<float> ccsMask;
  //ccsMask = findCCs( filepath, regionMatrix, casename, outputDir );//get the ccs mask



  return;
}


