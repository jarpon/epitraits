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
#include <stopwatch.h>

#define TRACE
#include <trace.h>
using namespace std;

extern void testsStatisticalTests();
extern VoxelMatrix<float> findNucleus(const VoxelMatrix<float>&);
//extern VoxelMatrix<float> findNucleusCascadeMethod(const VoxelMatrix<float>&,
//                                                   const string&, const string&);
extern VoxelMatrix<float> findMoreNuclei(const VoxelMatrix<float>&,
                                         const string&, const string&);
extern VoxelMatrix<float> findNucleusAlternative(const VoxelMatrix<float>&);
extern void nucleusAnalysis(const VoxelMatrix<float>&, VoxelMatrix<float>&,
                            const string&, const string&, const int&, DataSet&);
extern VoxelMatrix<float> findCCs( const VoxelMatrix<float>&, VoxelMatrix<float>&,
                                   const string&, const string&);
extern VoxelMatrix<float> findNucleoli( const VoxelMatrix<float>&, VoxelMatrix<float>&,
                                   const string&, const string&);
extern VoxelMatrix<float> findCCs16bits( const VoxelMatrix<float>&, VoxelMatrix<float>&,
                                   const string&, const string&);
extern void findCCsManually( const VoxelMatrix<float>&, VoxelMatrix<float>&,
                                   const string&, const string&, const string&);
extern VoxelMatrix<float> findChromosomes( const VoxelMatrix<float>&, VoxelMatrix<float>&,
                                   const string&, const string&);
extern VoxelMatrix<float> findGenes(const VoxelMatrix<float>&, VoxelMatrix<float>&,
                                   const string&, const string&);
extern void chromocentersAnalysis(VoxelMatrix<float>&, const string&, const string&,
                                  const int&, int&, DataSet&, DataSet&, DataSet&);
extern void chromocentersInterdistances(const string&,
                           const string&, DataSet&);
extern void chromosomesAnalysis(VoxelMatrix<float>&, const string&, const string&,
                                  const int&, int&, DataSet&, DataSet&, DataSet&);
//extern void nucleoliAnalysis(VoxelMatrix<float>&, const string&, const string&,
//                                  int&, DataSet&, DataSet&);
//extern void spatialModelAnalysis(TriMesh<float>&,
//                          const string&, const string&, const int& );
//extern void spatialModelEvaluator(const string&, const string&, const string&, const int,
//                          DataSet&, RandomGenerator&);
extern void realDataEvaluator(const string&, const string&, const string&, const int&,
                          DataSet&, RandomGenerator&);
extern void twoCompartmentsEvaluator(const string&, const string&, const string&, const int&,
                          DataSet&, RandomGenerator&);
//extern void nucleoliEvaluator(const string&, const string&, const string&, const int&,
//                          DataSet&, RandomGenerator&);
extern void uniformTest(const string&, const string&, DataSet&);
extern VoxelMatrix<float> isolateNuclei(const VoxelMatrix<float>&,
                                        const string&, const string&);
//extern VoxelMatrix<float> unifyLabels( const VoxelMatrix<float>&,
//                                           const int&, const int&);
extern VoxelMatrix<float> unifyLabels( const VoxelMatrix<float>&);
extern void doIt(const string&, const string&, RandomGenerator&);
extern void doIt2(const string&, const string&);
extern void test2Distributions(const string&, const string&, const string& = 0);

//extern void analyzeSample(const string&, const int&, int&, DataSet&, DataSet&);
int main(int argc, char* argv[])
{
//  Stopwatch stopWatch;
//  stopWatch.start( "Global process" );

  RandomGenerator randomGenerator;
  string filepath, filename, parentDir, originalDir, originalVMDir, nucleiDir, chromocentersDir, intermediateProcessesDir, analysisDir, shapesDir;


  if ( (argc == 1) || ( argv[1] == std::string("-h") ) || ( argv[1] == std::string("--help") ) )
  {
    // Check the value of argc. If not enough parameters have been passed, inform user and exit.
    // Inform the user of how to use the program
    cout << "                                                         " << endl;
    cout << "Usage: epitraits [OPTION process to apply] FILES... " << endl;
    cout << "Run a global or specific process of segmentation and/or quantification of nuclei" << endl;
    cout << "            without any option, it will run all processes" << endl;
    cout << "                                                         " << endl;
    cout << "       -p   <process>, sets the process to execute if not all process are needed" << endl;
    cout << "                     , previous processed are assumed already done" << endl;
    cout << "                                                         " << endl;
    cout << "        segmentation and quantification              " << endl;
    cout << "           '1' to segment the nucleus" << endl;
    cout << "           '1+' to segment more than one nucleus per image" << endl;
    cout << "           '1c' to segment the nucleus with a cascade method" << endl;
    cout << "           '2' to analyze and quantify the nucleus" << endl;
    cout << "           '3' to segment the chromocenters" << endl;
    cout << "           '3m' to segment manually the chromocenters; introduce -1 to discard the stack" << endl;
    cout << "           '4' to analyze and quantify chromocenters" << endl;
    cout << "           '4-unify' to unify labels (i.e. chromocenter which appears as 2 different labels" << endl;
    cout << "           '4-interdistances' to quantify interdistances among chromocenters" << endl;
    cout << "           '5' to generate spatial models taking into account:" << endl;
    cout << "                 - real distances to the border" << endl;
    cout << "                 - real equivalent radius" << endl;
    cout << "               >> this function must call files inside segmented_chromocenters folder" << endl;
    cout << "                                                         " << endl;
    cout << "           '8' to find and segment genes" << endl;
    cout << "           '9' to find and segment the chromosome" << endl;
    cout << "          '10' to analyze chromosomes" << endl;
    cout << "          '11' to segment the nucleoli" << endl;
    cout << "          '12' to analyze the nucleoli" << endl;
    cout << "                                                          " << endl;
    cout << "         spatial analysis (test and process real data)              " << endl;
    cout << "         to test model methods themself:" << endl;
    cout << "           '6' and after choose descriptor and constraints" << endl;
    cout << "               '1' to use F-function" << endl;
    cout << "               '2' to use G-function" << endl;
    cout << "               '3' to use H-function" << endl;
    cout << "               '4' to use B-function" << endl;
    cout << "               '5' to use C-function" << endl;
    cout << "                  '0' to not use de complete random model" << endl;
    cout << "                  '1' to use hardcore distances to constrain the model" << endl;
    cout << "                  '2' to use distances to the border constraints" << endl;
    cout << "                  '3' to use the hardcore and the boundary constraints" << endl;
    cout << "                  '4' to use maximal repulsion model" << endl;
    cout << "                                                         " << endl;
    cout << "          to study real data:" << endl;
    cout << "           '7' (to study chromocenters organization) and after choose descriptor and model" << endl;
    cout << "           '7-2' (to study two kinds of organization) and after choose descriptor and model" << endl;
    cout << "           '7-nucleoli' (to study nucleoli organization) and after choose descriptor and model" << endl;
    cout << "               '1' to use function F" << endl;
    cout << "               '2' to use function G" << endl;
    cout << "               '3' to use function H" << endl;
    cout << "               '4' to use function B - distance to the border" << endl;
    cout << "               '5' to use function C - distance to the centroid" << endl;
    cout << "               'all' to use all of them" << endl;
    cout << "                  '0' to not use constraints" << endl;
    cout << "                  '1' to use sized constraints ~ random" << endl;
    cout << "                  '2' to use distances to the border constraints" << endl;
    cout << "                  '3' to use both constraints" << endl;
    cout << "                  '4' to use maximal repulsion" << endl;
    cout << "                  '5' to use maximal repulsion with distance to the border constraint" << endl;
    cout << "                                                         " << endl;

    cout << "       FILES vm files" << endl;
    cout << "**********************************************************************************" << endl;
    cout << "   The following folders must exist:" << endl;
    cout << "         originals_vm" << endl;
    cout << "         segmented_nuclei" << endl;
    cout << "         segmented_chromocenters / segmented_chromosomes / segmented_genes / segmented_nucleoli" << endl;
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
    Stopwatch stopWatch;
    stopWatch.start( "Global process" );

    ENTER("Complete nuclei and chromocenters segmentation and analysis");
    // we set datafiles and file string
    DataSet nucleiDataset;
    DataSet chromocentersDataset;
    DataSet nucleoliDataset;

    //these are just counters to use at the datafiles
    int numNucleus = 0;
    int totalNumCCs = 0;
    int totalNumNucleoli = 0;

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
        nucleiDataset.save(analysisDir + "nuclei.csv", true );
        ++numNucleus;
        //EVAL("Nucleus analysis done!");
        VoxelMatrix<float> ccsMask;
        ccsMask = findCCs( originalVoxelMatrix, nucleusMask, filename, intermediateProcessesDir);
        ccsMask.save ( chromocentersDir + filename + ".vm", true );
        chromocentersAnalysis( ccsMask, filename, parentDir, numNucleus, totalNumCCs, nucleiDataset, chromocentersDataset, individualChromocentersDataset );
        individualChromocentersDataset.save(analysisDir + filename + "_chromocenters.csv", true );
        //++numNucleus;
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
            ( argv[2] == std::string("1") || argv[2] == std::string("1a") || argv[2] == std::string("1+") ||  argv[2] == std::string("1c")
              || argv[2] == std::string("2")
              || argv[2] == std::string("3") || argv[2] == std::string("3_16b") || argv[2] == std::string("3m")
              || argv[2] == std::string("4") || argv[2] == std::string("4-interdistances") || argv[2] == std::string("4-unify")
              || argv[2] == std::string("8") || argv[2] == std::string("9") || argv[2] == std::string("10")  || argv[2] == std::string("11") || argv[2] == std::string("12") )
            && ( argc > 3 ) )
  {

    Stopwatch stopWatch;
    stopWatch.start( "Global process" );

    DataSet nucleiDataset;
    DataSet chromocentersDataset;
    DataSet chromosomesDataset;
    DataSet nucleoliDataset;

    int numNucleus = 0;
    int totalNumCCs = 0;
    int totalNumChromosomes = 0;
    int totalNumNucleoli = 0;

    const string process = argv[2];
    EVAL(process);

    for (int i = 3; i < argc; i++)
    {
      filename = argv[i];

      // only vm files are processed
      if(filename.substr(filename.find_last_of(".") + 1) == "vm")
      {
//        DataSet individualChromocentersDataset;


        EVAL(filename);
        FileInfo fileInfo (filename);

        if ( fileInfo.isAbsolutePath() == false )
        {
          originalDir = getcwd( argv[i], 4048 );
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
        originalVMDir = parentDir + "/originals_vm/";
        intermediateProcessesDir = parentDir + "/intermediate_processes/";
        chromocentersDir = parentDir + "/segmented_chromocenters/";
        analysisDir = parentDir + "/analysis/";

        EVAL(filename);
        EVAL(originalVMDir);

        if ( process == "1" )
        {
          ENTER("Nuclei segmentation");
          VoxelMatrix<float> originalVoxelMatrix( originalVMDir + filename + ".vm" );
          VoxelMatrix<float> nucleusMask;
          nucleusMask = findNucleus( originalVoxelMatrix );
          nucleusMask.save ( nucleiDir + filename + ".vm", true );
          LEAVE();
        }

        if ( process == "1+" )
        {
          ENTER("Nuclei segmentation, looking for more than 1 nucleus");
          VoxelMatrix<float> originalVoxelMatrix( originalVMDir + filename + ".vm" );
          VoxelMatrix<float> nucleusMask;
          nucleusMask = findMoreNuclei( originalVoxelMatrix, filename, nucleiDir );
          //nucleusMask.save ( nucleiDir + filename + ".vm", true );
          LEAVE();
        }

        if ( process == "1a" )
        {
          ENTER("Nuclei alternative segmentation");
          VoxelMatrix<float> originalVoxelMatrix( originalVMDir + filename + ".vm" );
          VoxelMatrix<float> nucleusMask;
          nucleusMask = findNucleusAlternative( originalVoxelMatrix );
          nucleusMask.save ( nucleiDir + filename + ".vm", true );
          LEAVE();
        }

        if ( process == "1c" )
        {
          ENTER("Nuclei segmentation with cascade method");
          VoxelMatrix<float> originalVoxelMatrix( originalVMDir + filename + ".vm" );
          VoxelMatrix<float> nucleusMask;
          //nucleusMask = findNucleusCascadeMethod( originalVoxelMatrix, filename, nucleiDir );
          nucleusMask.save ( nucleiDir + filename + ".vm", true );
          LEAVE();
        }

        else if ( process == "2" )
        {
          ENTER("Nuclei quantification");
          string originalName = filename.substr( 0,filename.find_last_of("-")  );
          EVAL(originalName);
          VoxelMatrix<float> originalVoxelMatrix( originalVMDir + originalName + ".vm" );
          VoxelMatrix<float> nucleusMask ( nucleiDir + filename + ".vm" );
          nucleusAnalysis( originalVoxelMatrix, nucleusMask, filename, parentDir, numNucleus, nucleiDataset );
          nucleiDataset.save(analysisDir + "nuclei.csv", true );
          ++numNucleus;
          LEAVE();
        }

        else if ( process == "3" )
        {
          ENTER("Chromocenters segmentation");
          string originalName = filename.substr( 0,filename.find_last_of("-")  );
          VoxelMatrix<float> originalVoxelMatrix( originalVMDir + originalName + ".vm" );
          //VoxelMatrix<float> originalVoxelMatrix( originalVMDir + filename + ".vm" );
          VoxelMatrix<float> nucleusMask ( nucleiDir + filename + ".vm" );
          VoxelMatrix<float> ccsMask;
          ccsMask = findCCs( originalVoxelMatrix, nucleusMask, filename, intermediateProcessesDir);
          ccsMask.save ( chromocentersDir + filename + ".vm", true );
          LEAVE();
        }

        else if ( process == "3_16b" )
        {
          ENTER("Chromocenters segmentation (16bits images)");
          string originalName = filename.substr( 0,filename.find_last_of("-")  );
          VoxelMatrix<float> originalVoxelMatrix( originalVMDir + originalName + ".vm" );
          //VoxelMatrix<float> originalVoxelMatrix( originalVMDir + filename + ".vm" );
          VoxelMatrix<float> nucleusMask ( nucleiDir + filename + ".vm" );
          VoxelMatrix<float> ccsMask;
          ccsMask = findCCs16bits( originalVoxelMatrix, nucleusMask, filename, intermediateProcessesDir);
          ccsMask.save ( chromocentersDir + filename + ".vm", true );
          LEAVE();
        }

        else if ( process == "3m" )
        {
          ENTER("Chromocenters manual segmentation");
          string originalName = filename.substr( 0,filename.find_last_of("-")  );
          VoxelMatrix<float> originalVoxelMatrix( originalVMDir + originalName + ".vm" );
          //VoxelMatrix<float> originalVoxelMatrix( originalVMDir + filename + ".vm" );
          VoxelMatrix<float> nucleusMask ( nucleiDir + filename + ".vm" );
          //VoxelMatrix<float> ccsMask;
          findCCsManually( originalVoxelMatrix, nucleusMask, filename, intermediateProcessesDir, chromocentersDir);
          //ccsMask = findCCsManually( originalVoxelMatrix, nucleusMask, filename, intermediateProcessesDir);
          //ccsMask.save ( chromocentersDir + filename + ".vm", true );
          LEAVE();
        }

        else if ( process == "4" )
        {
          ENTER("Chromocenters quantification");
          DataSet individualChromocentersDataset;
          VoxelMatrix<float> ccsMask ( chromocentersDir + filename + ".vm" );
          chromocentersAnalysis( ccsMask, filename, parentDir, numNucleus, totalNumCCs, nucleiDataset, chromocentersDataset, individualChromocentersDataset );
          individualChromocentersDataset.save(analysisDir + filename + "_chromocenters.csv", true );
          ++numNucleus;
          LEAVE();
        }

        else if ( process == "4-interdistances" )
        {
          ENTER("Chromocenters interdistances quantification");
          DataSet individualChromocentersDataset;
          //VoxelMatrix<float> ccsMask ( chromocentersDir + filename + ".vm" );
          chromocentersInterdistances( filename, analysisDir, individualChromocentersDataset );
          individualChromocentersDataset.save(analysisDir + filename + "_chromocentersInterdistances.csv", true );
          LEAVE();
        }

        else if ( process == "4-unify" )
        {
          ENTER("Unifying labeling of two different regions");
//          const int oldLabel = argv[3];
//          const int newLabel = argv[4];
//          stringstream oldRegionLabel(argv[3]);
//          int oldLabel;
//          oldRegionLabel >> oldLabel;
//          stringstream newRegionLabel(argv[4]);
//          int newLabel;
//          newRegionLabel >> newLabel;
          const VoxelMatrix<float> ccsMask ( chromocentersDir + filename + ".vm" );
          //VoxelMatrix<float> changedVM = unifyLabels( ccsMask, oldLabel, newLabel );
          VoxelMatrix<float> changedVM = unifyLabels( ccsMask );
          //changedVM.save( chromocentersDir + filename + "_unified.vm" );
          changedVM.save( chromocentersDir + filename + ".vm", true );
          LEAVE();
        }

        else if ( process == "8" )
        {
          ENTER("Genes' search");
          string genesDir = parentDir + "/segmented_genes/";
          string originalName = filename.substr( 0,filename.find_last_of("-")  );
          VoxelMatrix<float> originalVoxelMatrix( originalVMDir + originalName + ".vm" );
          VoxelMatrix<float> nucleusMask ( nucleiDir + filename + ".vm" );
          VoxelMatrix<float> genes;
          genes = findGenes( originalVoxelMatrix, nucleusMask, filename, intermediateProcessesDir );
          //genes = findGenes( originalVoxelMatrix, originalVoxelMatrix, filename, intermediateProcessesDir );
          genes.save ( genesDir + filename + ".vm", true );
          LEAVE();
        }

        else if ( process == "9" )
        {
          ENTER("Chromosome segmentation");
          string chromosomesDir = parentDir + "/segmented_chromosomes/";
          string originalName = filename.substr( 0,filename.find_last_of("-")  );
          VoxelMatrix<float> originalVoxelMatrix( originalVMDir + originalName + ".vm" );
          VoxelMatrix<float> nucleusMask ( nucleiDir + filename + ".vm" );
          VoxelMatrix<float> chromosomeMask;
          chromosomeMask = findChromosomes( originalVoxelMatrix, nucleusMask, filename, intermediateProcessesDir);
          chromosomeMask.save ( chromosomesDir + filename + ".vm", true );
          LEAVE();
        }

        else if ( process == "10" )
        {
          ENTER("Chromosomes quantification");
          DataSet individualChromosomesDataset;
          string chromosomesDir = parentDir + "/segmented_chromosomes/";
          VoxelMatrix<float> chromosomesMask ( chromosomesDir + filename + ".vm" );
          EVAL(chromosomesMask.max().max().max());
          chromosomesAnalysis( chromosomesMask, filename, parentDir, numNucleus, totalNumChromosomes, nucleiDataset, chromosomesDataset, individualChromosomesDataset );
          individualChromosomesDataset.save(analysisDir + filename + "_chromosomes.csv", true );
          ++numNucleus;
          LEAVE();
        }

        else if ( process == "11" )
        {
          ENTER("Nucleoli segmentation");
          string nucleoliDir = parentDir + "/segmented_nucleoli/";
          string originalName = filename.substr( 0,filename.find_last_of("-")  );
          VoxelMatrix<float> originalVoxelMatrix( originalVMDir + originalName + ".vm" );
          VoxelMatrix<float> nucleusMask ( nucleiDir + filename + ".vm" );
          VoxelMatrix<float> nucleoliMask;
          nucleoliMask = findNucleoli( originalVoxelMatrix, nucleusMask, filename, intermediateProcessesDir);
          nucleoliMask.save ( nucleoliDir + filename + ".vm", true );
          LEAVE();
        }

        else if ( process == "12" )
        {
          ENTER("Nucleoli quantification");
          string nucleoliDir = parentDir + "/segmented_nucleoli/";
          DataSet individualNucleoliDataset;
          string originalName = filename.substr( 0,filename.find_last_of("-")  );
          VoxelMatrix<float> nucleoliMask ( nucleoliDir + filename + ".vm" );
          //nucleoliAnalysis( nucleoliMask, filename, parentDir, totalNumNucleoli, nucleoliDataset, individualNucleoliDataset );
          individualNucleoliDataset.save(analysisDir + filename + "_nucleoli.csv", true );
          LEAVE();
        }

      }
      else cout << "Error opening the image, it must be a VoxelMatrix image " << endl;

      if ( argv[2] == std::string("4") )
      {
        nucleiDataset.save(analysisDir + "nuclei_extended.csv", true );
        chromocentersDataset.save(analysisDir + "ccs.csv", true );
      }
    }

    if ( argv[2] == std::string("2") )
      nucleiDataset.save(analysisDir + "nuclei.csv", true );

    else if ( argv[2] == std::string("4") )
    {
      nucleiDataset.save(analysisDir + "nuclei_extended.csv", true );
      chromocentersDataset.save(analysisDir + "ccs.csv", true );
    }

    else if ( argv[2] == std::string("10") )
    {
      nucleiDataset.save(analysisDir + "nuclei_extended.csv", true );
      chromosomesDataset.save(analysisDir + "ccs.csv", true );
    }

    else if ( argv[2] == std::string("12") )
      nucleoliDataset.save(analysisDir + "nucleoli.csv", true );

    stopWatch.stop( "Global process" );
    stopWatch.print();


  }

  /*! Generates spatial models
  ****************************************************************/
  else if ( ( argv[1] == std::string("-p") ) && ( argv[2] == std::string("5") ) && ( argc > 3 ) )
  {

    Stopwatch stopWatch;
    stopWatch.start( "Global process" );

    ENTER("Spatial model generator");

    //introduce HERE the number of patterns to generate for each nucleus;
    const int numPatterns = 10;

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
        //TriMesh<float> nucleusTriMesh ( shapesDir + filename + "_nucleus.tm" );
        TriMesh<float> nucleusTriMesh ( shapesDir + filename + ".tm" );

//        spatialModelAnalysis( nucleusTriMesh, filename, parentDir, numPatterns );

      }
      else cout << "Error" << endl;
    }
    LEAVE();

    stopWatch.stop( "Global process" );
    stopWatch.print();

  }


  /*! Evaluates spatial descriptors
  ****************************************************************/
  else if ( ( argv[1] == std::string("-p") ) &&
            ( argv[2] == std::string("6") || argv[2] == std::string("7") || argv[2] == std::string("7-2") || argv[2] == std::string("7-nucleoli") || argv[2] == std::string("7-2distributions") ) &&
            ( argv[3] == std::string("1") || argv[3] == std::string("2") || argv[3] == std::string("3") || argv[3] == std::string("4") || argv[3] == std::string("5") || argv[3] == std::string("all") ) &&
//            ( argv[4] == std::string("0") || argv[4] == std::string("1") || argv[4] == std::string("2") || argv[4] == std::string("3") ) || argv[4] == std::string("4") ) &&
            ( argc > 5 ) )
  {

    Stopwatch stopWatch;
    stopWatch.start( "Global process" );

    ENTER("Spatial model analysis");

    DataSet dataSet, errorDataSet;
    int numErrors;

    //introduce HERE the number of patterns to generate for each nucleus;
//    const int numPatterns = 99;
    string test, function;
    int constraints;

    //if ( argv[2] == std::string("6") )       test = "model";
    if ( argv[2] == std::string("6") )       test = "2compartments";
    else if ( argv[2] == std::string("7") )  test = "data";
    else if ( argv[2] == std::string("7-nucleoli") )  test = "nucleoli";
    else if ( argv[2] == std::string("7-2") )  test = "2distributions";

    if ( argv[3] == std::string("2") )       function = "G";
    else if ( argv[3] == std::string("3") )  function = "H";
    else if ( argv[3] == std::string("4") )  function = "B";
    else if ( argv[3] == std::string("5") )  function = "C";
    else if ( argv[3] == std::string("all") )  function = "all";
    else //if ( argv[3] == std::string("1") )
                                             function = "F";

    if ( argv[4] == std::string("0") )       constraints = 0;
    else if ( argv[4] == std::string("1") )  constraints = 1;
    else if ( argv[4] == std::string("2") )  constraints = 2;
    else if ( argv[4] == std::string("3") )  constraints = 3;
    else if ( argv[4] == std::string("4") )  constraints = 4;
    else if ( argv[4] == std::string("5") )  constraints = 5;



    // if we got enough parameters and options...
    for ( int i = 5 ; i < argc; ++i)
    {

      Stopwatch stopWatchNucleus;
      stopWatchNucleus.start( "Nucleus processed" );

      try
      {

        filename = argv[i];

        // only vm files are processed
        if(filename.substr(filename.find_last_of(".") + 1) == "vm")
        {
          EVAL(filename);
          FileInfo fileInfo (filename);

          if ( fileInfo.isAbsolutePath() == false )
          {
            originalDir = getcwd( argv[i], 4048 );
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

//          if  ( test == "model" )
//          {
//            ostringstream oss;
//            oss << constraints;
//            spatialModelEvaluator( filename, parentDir, function, constraints, dataSet, randomGenerator );
////            dataSet.save(analysisDir + oss.str() + "/" + function + "/" + filename + "_model.csv", true );
////            dataSet.save(analysisDir + "indexes_" + oss.str() + function + oss.str() + ".csv", true );
//          }
          if  ( test == "2distributions" )
          {
            twoCompartmentsEvaluator( filename, parentDir, function, constraints, dataSet, randomGenerator );
            ostringstream oss;
            oss << constraints;
            dataSet.save(analysisDir + "indexes_" + function + oss.str() + ".csv", true );
          }
          else if ( test == "data" )
          {
            //string originalName = filename.substr( 0,filename.find_last_of("-")  );
            //realDataEvaluator( filename, parentDir, function, dataSet, randomGenerator );
            realDataEvaluator( filename, parentDir, function, constraints, dataSet, randomGenerator );

            //realDataEvaluator( originalName, parentDir, function, constraints, dataSet, randomGenerator );

            ostringstream oss;
            oss << constraints;
            dataSet.save(analysisDir + "indexes_" + function + oss.str() + ".csv", true );
//            dataSet.save(analysisDir + "indexes_" + function + oss.str() + "_random.csv", true );
          }
          else if ( test == "nucleoli" )
          {
            //string originalName = filename.substr( 0,filename.find_last_of("-")  );
            //realDataEvaluator( filename, parentDir, function, dataSet, randomGenerator );
//            nucleoliEvaluator( filename, parentDir, function, constraints, dataSet, randomGenerator );

            //realDataEvaluator( originalName, parentDir, function, constraints, dataSet, randomGenerator );

            ostringstream oss;
            oss << constraints;
            dataSet.save(analysisDir + "indexes_" + oss.str() + function + oss.str() + ".csv", true );
          }
          else if ( test == "2distributions" )
          {

            test2Distributions( filename, filename );

          }
          else
          {
            ProgramError error;
            error.setWhat( "Inputs are wrong" );
            error.setWhat( "Call epitraits help for more information" );
            throw error;
          }
        }
        else
          cout << "Error" << endl;

        if ( test == "data" )
        {
          ostringstream oss;
          oss << constraints;
         // dataSet.save(analysisDir + "indexes_" + oss.str() + function + oss.str() + ".csv", true );

        }

      }
      catch(Exception exception)
      {
        EVAL(exception.getWhat());
        ostringstream oss;
        oss << constraints;
        errorDataSet.setValue("id",numErrors,numErrors+1);
        errorDataSet.setValue("filename",numErrors,filename);
        //errorDataSet.save(analysisDir + "failedNuclei" + function + oss.str() + ".csv", true );
        errorDataSet.save(analysisDir + "failedNuclei" + function + oss.str() + "_random.csv", true );
        ++numErrors;

      }

      stopWatchNucleus.stop( "Nucleus processed" );
      stopWatchNucleus.print();
    }

    LEAVE();

    stopWatch.stop( "Global process" );
    stopWatch.print();

  }

  /*! Temporal
  ****************************************************************/
  else if ( ( argv[1] == std::string("-p") ) && ( argv[2] == std::string("99") ) && ( argc > 3 ) )
  {

    Stopwatch stopWatch;
    stopWatch.start( "Global process" );

    ENTER("Temporal To Do Things!");

    // if we got enough parameters and options...
    for (int i = 3; i < argc; i++)
    {
      filename = argv[i];

      // only vm files are processed

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

        //const VoxelMatrix<float> originalVoxelMatrix( filepath );

        //analysisDir = parentDir + "/analysis/";
        EVAL(filename);
        EVAL(parentDir);

        //test2Distributions( filepath, filename, filename );

        testsStatisticalTests();
        //doIt2( filename, parentDir );
//        VoxelMatrix<float>nucleiMask = isolateNuclei( originalVoxelMatrix );
//        nucleiMask.save( "/home/jarpon/Desktop/" + filename + "-nucleus.vm", true );

        //doIt( filename, parentDir, randomGenerator);


    }
    LEAVE();

    stopWatch.stop( "Global process" );
    stopWatch.print();

  }

  /*! Crio sections
  ****************************************************************/
  else if ( ( argv[1] == std::string("-p") ) && ( argv[2] == std::string("extractnuclei") ) && ( argc > 3 ) )
  {

    Stopwatch stopWatch;
    stopWatch.start( "Global process" );

    ENTER("Temporal To Do Things!");

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

        nucleiDir = parentDir + "/segmented_nuclei/";

        const VoxelMatrix<float> originalVoxelMatrix( filepath );
        VoxelMatrix<float> nucleiMask;
        nucleiMask = isolateNuclei( originalVoxelMatrix, filename, parentDir );
        nucleiMask.save ( nucleiDir + filename + ".vm", true );
      }
      else cout << "Error" << endl;
    }

    LEAVE();

    stopWatch.stop( "Global process" );
    stopWatch.print();

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
