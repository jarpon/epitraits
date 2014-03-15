#include <iostream>
#include <dataset.h>
#include <dirent.h>
#include <errno.h>
#include <fileinfo.h>
//#include <QString>
#include <unistd.h>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <voxelmatrix.h>
#include <trimesh.h>
#include <sstream>
#include <programerror.h>

#define TRACE
#include <trace.h>
using namespace std;

extern VoxelMatrix<float> findNucleus(const VoxelMatrix<float>&);
extern void nucleusAnalysis(const VoxelMatrix<float>&, VoxelMatrix<float>&,
                            const string&, const string&, const int&, DataSet&);
extern VoxelMatrix<float> findCCs( const VoxelMatrix<float>&, const VoxelMatrix<float>&,
                                   const string&, const string&);
extern void chromocentersAnalysis(VoxelMatrix<float>&, const string&, const string&,
                                  const int&, int&, DataSet&, DataSet&, DataSet&);
extern void spatialModelAnalysis(TriMesh<float>&,
                          const string&, const string&, const int& );
extern void spatialModelEvaluator(const string&, const string&, const string&, const int,
                          DataSet&, RandomGenerator&);
extern void realDataEvaluator(const string&, const string&, const string&, const int&,
                          DataSet&, RandomGenerator&);
extern void uniformTest( const string&, const string&, DataSet&);

//extern void analyseSample(const string&, const int&, int&, DataSet&, DataSet&);

int main(int argc, char* argv[])
{
  RandomGenerator randomGenerator;
  string filepath, filename, parentDir, originalDir, nucleiDir, chromocentersDir, intermediateProcessesDir, analysisDir, shapesDir;


  if ( (argc == 1) || ( argv[1] == std::string("-h") ) || ( argv[1] == std::string("--help") ) )
  {
    // Check the value of argc. If not enough parameters have been passed, inform user and exit.
    // Inform the user of how to use the program
    cout << "Usage: epitraits [OPTION process to apply] FILES... " << endl;
    cout << "Run a global or specific process of segmentation and/or quantification of nuclei" << endl;
    cout << "            without any option, it will run all processes" << endl;
    cout << "                                                         " << endl;
    cout << "       -p   <process>, sets the process to execute if not all process are needed" << endl;
    cout << "                     , previous processed are assumed already done" << endl;
    cout << "           '1' to segment the nucleus" << endl;
    cout << "           '2' to analyse and quantify the nucleus" << endl;
    cout << "           '3' to segment the chromocenters" << endl;
    cout << "           '4' to analyse and quantify chromocenters" << endl;
    cout << "           '5' to generate spatial models taking into account:" << endl;
    cout << "                 - real distances to the border" << endl;
    cout << "                 - real equivalent radius" << endl;
    cout << "               >> this function must call files inside segmented_chromocenters folder" << endl;
    cout << "                                                          " << endl;
    cout << "       to test model methods themself:" << endl;
    cout << "           '6' and after choose descriptor and constraints" << endl;
    cout << "               '1' to use function F" << endl;
    cout << "               '2' to use function G" << endl;
    cout << "               '3' to use function H" << endl;
    cout << "               '4' to use distance to the border descriptor" << endl;
    cout << "                  '0' to not use constraints" << endl;
    cout << "                  '1' to use sized constraints" << endl;
    cout << "                  '2' to use distances to the border constraints" << endl;
    cout << "                  '3' to use both constraints" << endl;
    cout << "                                                         " << endl;
    cout << "       to test real data:" << endl;
    cout << "           '7' and after choose descriptor" << endl;
    cout << "               '1' to use function F" << endl;
    cout << "               '2' to use function G" << endl;
    cout << "               '3' to use function H" << endl;
    cout << "               '4' to use distance to the border descriptor" << endl;
    cout << "                                                         " << endl;
    cout << "       FILES vm files" << endl;
    cout << "**********************************************************************************" << endl;
    cout << "   The following folders must exist:" << endl;
    cout << "         segmented_nuclei" << endl;
    cout << "         segmented_chromocenters" << endl;
    cout << "         intermediate_processes" << endl;
    cout << "         analysis" << endl;
    cout << "         shapes" << endl;
  }

  /*! Run all the process from 1 to 4 on each image of the folder.
   *  This means: nucleus's segmentation & quantification
   *              chromocenters' segmentation & quantification
  ****************************************************************/
  else if ( argv[1] != std::string("-p") )
  {
    ENTER("Complete nuclei and chromocenters segmentation and analysis");
    // we set datafiles and file string
    DataSet nucleiDataset;
    DataSet chromocentersDataset;

    //these are just counters to use at the datafiles
    int numNucleus = 0;
    int totalNumCCs = 0;

    //if we got enough parameters but no options
    for (int i = 1; i < argc; i++)
    {
      filename = argv[i];
      DataSet individualChromocentersDataset;

      // only vm files are processed
      if(filename.substr(filename.find_last_of(".") + 1) == "vm")
      {
        FileInfo fileInfo (filename);
        if ( fileInfo.isAbsolutePath() == false )
        {
          originalDir = getcwd( argv[i], 2048 );
          originalDir = originalDir + "/";
          filepath = originalDir + filename;
          FileInfo fileInfo (filepath);
          parentDir = originalDir.substr(0,originalDir.find_last_of("/\\"));
          parentDir = parentDir.substr(0,parentDir.find_last_of("/\\")) + "/";
          filename = fileInfo.baseName();
        }

        else
        {
          filepath = filename;
          filename = fileInfo.baseName();
          originalDir = fileInfo.dirName();
          parentDir = originalDir.substr(0,originalDir.find_last_of("/\\"));
          parentDir = parentDir.substr(0,parentDir.find_last_of("/\\"));
          //originalDir = parentDir + "/originals_vm/";
        }

        EVAL(filename);

        nucleiDir = parentDir + "/segmented_nuclei/";
        intermediateProcessesDir = parentDir + "/intermediate_processes/";
        chromocentersDir = parentDir + "/segmented_chromocenters/";
        analysisDir = parentDir + "/analysis/";

        VoxelMatrix<float> originalVoxelMatrix( filepath );
        VoxelMatrix<float> nucleusMask;
        nucleusMask = findNucleus( originalVoxelMatrix );//get the nucleus mask
        nucleusMask.save ( nucleiDir + filename + ".vm", true );
        nucleusAnalysis( originalVoxelMatrix, nucleusMask, filename, parentDir, numNucleus, nucleiDataset );
        VoxelMatrix<float> ccsMask;
        ccsMask = findCCs( originalVoxelMatrix, nucleusMask, filename, intermediateProcessesDir);
        ccsMask.save ( chromocentersDir + filename + ".vm", true );
        chromocentersAnalysis( ccsMask, filename, parentDir, numNucleus, totalNumCCs, nucleiDataset, chromocentersDataset, individualChromocentersDataset );
        individualChromocentersDataset.save(analysisDir + filename + "_chromocenters.csv", true );
        ++numNucleus;
      }
      else cout << "Error opening the image, it must be a VoxelMatrix image " << endl;
    }

    // dataSet is saved into a file at the end of the function
    nucleiDataset.save(analysisDir + "nuclei.csv", true );
    chromocentersDataset.save(analysisDir + "chromocenters.csv", true );
    LEAVE();
  }

  /*! Run a process from 1 to 4 on each image of the folder.
  ****************************************************************/
  else if ( ( argv[1] == std::string("-p") ) &&
            ( argv[2] == std::string("1") || argv[2] == std::string("2") || argv[2] == std::string("3") || argv[2] == std::string("4") )
            && ( argc > 3 ) )
  {
    DataSet nucleiDataset;
    DataSet chromocentersDataset;
    int numNucleus = 0;
    int totalNumCCs = 0;

    const string process = argv[2];
    EVAL(process);

    for (int i = 3; i < argc; i++)
    {
      filename = argv[i];

      // only vm files are processed
      if(filename.substr(filename.find_last_of(".") + 1) == "vm")
      {
        DataSet individualChromocentersDataset;

        EVAL(filename);
        FileInfo fileInfo (filename);

        if ( fileInfo.isAbsolutePath() == false )
        {
          originalDir = getcwd( argv[i], 2048 );
          originalDir = originalDir + "/";
          filepath = originalDir + filename;
          FileInfo fileInfo (filepath);
          parentDir = originalDir.substr(0,originalDir.find_last_of("/\\"));
          parentDir = parentDir.substr(0,parentDir.find_last_of("/\\")) + "/";
          filename = fileInfo.baseName();
        }

        else
        {
          filepath = filename;
          filename = fileInfo.baseName();
          originalDir = fileInfo.dirName();
          parentDir = originalDir.substr(0,originalDir.find_last_of("/\\"));
          parentDir = parentDir.substr(0,parentDir.find_last_of("/\\"));
          //originalDir = parentDir + "/originals_vm/";
        }

        nucleiDir = parentDir + "/segmented_nuclei/";
        VoxelMatrix<float> originalVoxelMatrix( filepath );


        if ( process == "1" )
        {
          VoxelMatrix<float> nucleusMask;
          nucleusMask = findNucleus( originalVoxelMatrix );
          nucleusMask.save ( nucleiDir + filename + ".vm", true );
        }

        else if ( process == "2" )
        {
          VoxelMatrix<float> nucleusMask ( nucleiDir + filename + ".vm" );
          nucleusAnalysis( originalVoxelMatrix, nucleusMask, filename, parentDir, numNucleus, nucleiDataset );
          ++numNucleus;
        }

        else if ( process == "3" )
        {
          intermediateProcessesDir = parentDir + "/intermediate_processes/";
          chromocentersDir = parentDir + "/segmented_chromocenters/";
          analysisDir = parentDir + "/analysis/";
          VoxelMatrix<float> nucleusMask ( nucleiDir + filename + ".vm" );
          VoxelMatrix<float> ccsMask;
          ccsMask = findCCs( originalVoxelMatrix, nucleusMask, filename, intermediateProcessesDir);
          ccsMask.save ( chromocentersDir + filename + ".vm", true );
        }

        else if ( process == "4" )
        {
          chromocentersDir = parentDir + "/segmented_chromocenters/";
          analysisDir = parentDir + "/analysis/";
          VoxelMatrix<float> nucleusMask ( nucleiDir + filename + ".vm" );
          VoxelMatrix<float> ccsMask ( chromocentersDir + filename + ".vm" );
          chromocentersAnalysis( ccsMask, filename, parentDir, numNucleus, totalNumCCs, nucleiDataset, chromocentersDataset, individualChromocentersDataset );
          individualChromocentersDataset.save(analysisDir + filename + "_chromocenters.csv", true );
          ++numNucleus;
        }

      }
      else cout << "Error opening the image, it must be a VoxelMatrix image " << endl;
    }

    if ( argv[2] == std::string("2") )
      nucleiDataset.save(analysisDir + "nuclei.csv", true );

    else if ( argv[2] == std::string("4") )
    {
      nucleiDataset.save(analysisDir + "nuclei_extended.csv", true );
      chromocentersDataset.save(analysisDir + "chromocenters.csv", true );
    }

  }

  /*! Generates spatial models
  ****************************************************************/
  else if ( ( argv[1] == std::string("-p") ) && ( argv[2] == std::string("5") ) && ( argc > 3 ) )
  {
    ENTER("Spatial model generator");

    //introduce HERE the number of patterns to generate for each nucleus;
    const int numPatterns = 50;

    // if we got enough parameters and options...
    for (int i = 3; i < argc; i++)
    {
      filename = argv[i];

      // only vm files are processed
      if(filename.substr(filename.find_last_of(".") + 1) == "vm")
      {
        EVAL(filename);
        FileInfo fileInfo (filename);

        if ( fileInfo.isAbsolutePath() == false )
        {
          originalDir = getcwd( argv[i], 2048 );
          originalDir = originalDir + "/";
          filepath = originalDir + filename;
          FileInfo fileInfo (filepath);
          parentDir = originalDir.substr(0,originalDir.find_last_of("/\\"));
          parentDir = parentDir.substr(0,parentDir.find_last_of("/\\")) + "/";
          filename = fileInfo.baseName();
        }

        else
        {
          filepath = filename;
          filename = fileInfo.baseName();
          originalDir = fileInfo.dirName();
          parentDir = originalDir.substr(0,originalDir.find_last_of("/\\"));
          parentDir = parentDir.substr(0,parentDir.find_last_of("/\\"));
        }

        chromocentersDir = parentDir + "/segmented_chromocenters/";
        analysisDir = parentDir + "/analysis/";
        shapesDir = parentDir + "/shapes/";

        EVAL(filename);


        //VoxelMatrix<float> ccsMask ( chromocentersDir + filename + ".vm" );
        TriMesh<float> nucleusTriMesh ( shapesDir + filename + "_nucleus.tm" );

        spatialModelAnalysis( nucleusTriMesh, filename, parentDir, numPatterns );

      }
      else cout << "Error" << endl;
    }
    LEAVE();
  }


  /*! Evaluates spatial descriptors
  ****************************************************************/
  else if ( ( argv[1] == std::string("-p") ) &&
            ( argv[2] == std::string("6") || argv[2] == std::string("7") ) &&
            ( argv[3] == std::string("1") || argv[3] == std::string("2") || argv[3] == std::string("3") || argv[3] == std::string("4") ) &&
//            ( argv[4] == std::string("0") || argv[4] == std::string("1") || argv[4] == std::string("2") || argv[4] == std::string("3") ) &&
            ( argc > 5 ) )
  {
    ENTER("Spatial model analysis");

    DataSet dataSet;

    //introduce HERE the number of patterns to generate for each nucleus;
//    const int numPatterns = 99;
    string test, function;
    int constraints;

    if ( argv[2] == std::string("6") )       test = "model";
    else if ( argv[2] == std::string("7") )  test = "data";

    if ( argv[3] == std::string("2") )       function = "G";
    else if ( argv[3] == std::string("3") )  function = "H";
    else if ( argv[3] == std::string("4") )  function = "B";
    else //if ( argv[3] == std::string("1") )
                                             function = "F";

    if ( argv[4] == std::string("0") )       constraints = 0;
    else if ( argv[4] == std::string("1") )  constraints = 1;
    else if ( argv[4] == std::string("2") )  constraints = 2;
    else if ( argv[4] == std::string("3") )  constraints = 3;

    // if we got enough parameters and options...
    for ( int i = 5 ; i < argc; ++i)
    {

      filename = argv[i];

      // only vm files are processed
      if(filename.substr(filename.find_last_of(".") + 1) == "vm")
      {
        EVAL(filename);
        FileInfo fileInfo (filename);

        if ( fileInfo.isAbsolutePath() == false )
        {
          originalDir = getcwd( argv[i], 2048 );
          originalDir = originalDir + "/";
          filepath = originalDir + filename;
          FileInfo fileInfo (filepath);
          parentDir = originalDir.substr(0,originalDir.find_last_of("/\\"));
          parentDir = parentDir.substr(0,parentDir.find_last_of("/\\")) + "/";
          filename = fileInfo.baseName();
        }

        else
        {
          filepath = filename;
          filename = fileInfo.baseName();
          originalDir = fileInfo.dirName();
          parentDir = originalDir.substr(0,originalDir.find_last_of("/\\"));
          parentDir = parentDir.substr(0,parentDir.find_last_of("/\\"));
        }

        analysisDir = parentDir + "/analysis/";
        shapesDir = parentDir + "/shapes/";

        EVAL(filename);

        if  ( test == "model" )
        {
          ostringstream oss;
          oss << constraints;
          spatialModelEvaluator( filename, parentDir, function, constraints, dataSet, randomGenerator );
          dataSet.save(analysisDir + oss.str() + "/" + function + "/" + filename + "_model.csv", true );
        }
        else if ( test == "data" )
        {
          //realDataEvaluator( filename, parentDir, function, dataSet, randomGenerator );
          realDataEvaluator( filename, parentDir, function, constraints, dataSet, randomGenerator );
//          dataSet.save(analysisDir + filename + "_" + function + ".csv", true );
        }
        else
        {
          ProgramError error;
          error.setWhat( "Inputs are wrong" );
          error.setWhat( "Call epitraits help for more information" );
        }
      }
      else
        cout << "Error" << endl;
    }

    if ( test == "data" )
    {
      ostringstream oss;
      oss << constraints;
      dataSet.save(analysisDir + "pValues_" + oss.str() + function + oss.str() + ".csv", true );
    }
    LEAVE();
  }

  /*! Gets an error!
  ****************************************************************/
  else
  {
    ProgramError error;
    error.setWhat( "Error calling the program" );
    error.setWhat( "Call epitraits to see the help" );
  }

  return 0;
}
