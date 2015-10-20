#include <spatialmodelevaluator.h>
//#include "spatialmodelevaluator2.h"
#include <spatialdescriptorfunctionf.h>
#include <spatialdescriptorfunctiong.h>
#include <spatialdescriptorfunctionh.h>
#include <spatialdescriptorfunctionb.h>
#include <spatialdescriptorfunctionc.h>
#include <spatialmodelmaximalrepulsion3d.h>
//#include "spatialdescriptorborder.h"
//#include "spatialdescriptorcentroid.h"
//#include "maximalrepulsion.h"
#include "maxrepulsionwithdistances.h"
#include <spatialmodelcompleterandomness3d.h>
#include <spatialmodelborderdistance3d.h>
#include <spatialmodelhardcoreborderdistance3d.h>
#include <spatialmodelhardcoredistance3d.h>
#include <spatialmodel.h>
//#include "trimeshspatialmodel.h"
#include <trimesh.h>
#include <voxelmatrix.h>
#include <fileinfo.h>
//#include <regionanalysis.h>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <cmath>
#include <iostream>
#include <string>

#include <stringtools.h>
#include <fstream>
#include <sstream>
#include <unistd.h>
 #include <cstdlib>

#define TRACE
#include <trace.h>
using namespace std;


Vector<string> getNamesList(std::string inputDir)
{
  std::string fileTmp = inputDir + "nuclei.csv";

  EVAL(fileTmp);
  char sbuffer[200];
  char fileIn[200];
  char currentLine[250];
  std::ifstream inputFile;
  chdir(inputDir.c_str());
  //sprintf( sbuffer, "ls -1 *.csv > %s",  fileTmp.c_str() );
  //system( sbuffer );

  // ouvre le tmp crée et lit par ligne
  sprintf( fileIn, "%s", fileTmp.c_str());
  inputFile.open ( fileIn, std::ifstream::in );
  while ( inputFile.good() )
  {
    inputFile.getline( currentLine, 80 );
    if (inputFile.gcount () == 0)
    {
      PRINT("fichier vide")
          break;
    }
    PRINT("Dossier  en traitement : " );
    EVAL(currentLine);
  }
  Vector<string> list;
  return list;
}

Vector<string> loadFromCSV2( const int& rows, const int& columns )
{
  EVAL(rows);
  EVAL(columns);
    float data[rows][columns];
    std::ifstream file( "/home/jarpon/data/wt/all/analysis/nuclei.csv");

    for(int row = 0; row < rows; ++row)
    {
        std::string line;
        std::getline(file, line);
        if ( !file.good() )
            break;

        std::stringstream iss(line);
        //EVAL(line); //good

        for (int col = 0; col < columns; ++col)
        {
            std::string val;
            std::getline(iss, val, '\t');
            //EVAL(val);
            if ( !iss.good() )
            {
                break;
                  EVAL('here!');
            }

            std::stringstream convertor(val);
            convertor >> data[row][col];
            //EVAL(data[row][col])
    EVAL('here!');
        }
    }
    EVAL('here!');
    Vector<string> listNames;
    listNames.setSize( rows - 1 );
    for( int i=0; i<int(rows); i++ )
      listNames[i] = data[i+1][0];

    return listNames;
}


Vector<string> loadFromCSV( const std::string& filename )
{
  EVAL(filename);

    std::ifstream       file( filename.c_str() );
    std::vector< std::vector<std::string> >   matrix;
    std::vector<std::string>   row;
    std::string                line;
    std::string                cell;

    while( file )
    {
        std::getline(file,line);
        std::stringstream lineStream(line);
        row.clear();

        while( std::getline( lineStream, cell, ',' ) )
            row.push_back( cell );

        if( !row.empty() )
            matrix.push_back( row );
    }

    for( int i=0; i<int(matrix.size()); i++ )
    {
        for( int j=0; j<int(matrix[i].size()); j++ )
        {
            std::cout << matrix[i][j] << " ";
            EVAL(matrix[i][j])
        }

        std::cout << std::endl;
    }

    Vector<string> listNames;
    for( int i=0; i<int(matrix.size()); i++ )
      listNames[i] = matrix[0][i];

    return listNames;
}

void evaluator(
  const TriMesh<float>& nucleusTriMesh,
  TriMeshSpatialModel<float>& triMeshSpatialModel,
  const string& filename, const string& parentDir,
  const string& function, const int& constraints,
  DataSet& dataSet, RandomGenerator& randomGenerator)
{
  const string analysisDir = parentDir + "/analysis/";
  string classif = parentDir;
  classif = classif.substr( classif.find_last_of("/\\")+1, classif.length() );

  //new data
  const DataSet datasetNucleus( analysisDir + filename + "_chromocenters.csv" );
//  const DataSet datasetNucleus( analysisDir + filename + "_nucleoli.csv" );
  const DataSet globalAnalysis( analysisDir + "nuclei.csv" );
  DataSet copyGlobalAnalysis = globalAnalysis;



  Vector<string> nucleiNames;
  //nucleiNames.setSize( copyGlobalAnalysis.numRows() );
  nucleiNames = loadFromCSV2 ( copyGlobalAnalysis.numRows()+1, copyGlobalAnalysis.size()[1]+1);




  const int numPoints = datasetNucleus.size()[0];

  const int numPatterns = 99;

  SpatialModelEvaluator<float,float> modelEvaluator;
  modelEvaluator.setModel( triMeshSpatialModel );
  modelEvaluator.setNumMonteCarloSamples( numPatterns ); //to check uniformity
  modelEvaluator.setPrecision( 0.05 );

  SpatialDescriptor<float>* spatialDescriptor;

  DataSet saveTest;

  if ( function == "all" )
  {
    PRINT("all functions");
    SpatialModelCompleteRandomness3D<float> tempTriMeshSpatialModel;
    tempTriMeshSpatialModel.setRandomGenerator( randomGenerator );
    tempTriMeshSpatialModel.setTriMesh( nucleusTriMesh );
    tempTriMeshSpatialModel.initialize();
    Vertices<float> evaluationPositions = tempTriMeshSpatialModel.drawSample( 10000 );

    SpatialDescriptorFunctionF<float>* spatialDescriptorFunctionF;
    spatialDescriptorFunctionF = new SpatialDescriptorFunctionF<float>();
    spatialDescriptorFunctionF->setEvaluationPositions( evaluationPositions );
    spatialDescriptor = spatialDescriptorFunctionF;
    modelEvaluator.addDescriptor( *spatialDescriptor );

    spatialDescriptor = new SpatialDescriptorFunctionG<float>();
    modelEvaluator.addDescriptor( *spatialDescriptor );

    spatialDescriptor = new SpatialDescriptorFunctionH<float>();
    modelEvaluator.addDescriptor( *spatialDescriptor );

    SpatialDescriptorFunctionB<float>* spatialDescriptorFunctionB;
    spatialDescriptorFunctionB = new SpatialDescriptorFunctionB<float>();
    spatialDescriptorFunctionB->setTriMesh( nucleusTriMesh );
    spatialDescriptor = spatialDescriptorFunctionB;
    modelEvaluator.addDescriptor( *spatialDescriptor );

    SpatialDescriptorFunctionC<float>* spatialDescriptorFunctionC;
    spatialDescriptorFunctionC = new SpatialDescriptorFunctionC<float>();
    spatialDescriptorFunctionC->setCenter( nucleusTriMesh.cog() );
    spatialDescriptor = spatialDescriptorFunctionC;
    modelEvaluator.addDescriptor( *spatialDescriptor );

  }
  else if ( function == "G" )
  {
    PRINT("G");
    spatialDescriptor = new SpatialDescriptorFunctionG<float>();
  }
  else if ( function == "H" )
  {
    PRINT("H");
    spatialDescriptor = new SpatialDescriptorFunctionH<float>();
  }
  else if ( function == "B" )
  {
    PRINT("B");
    SpatialDescriptorFunctionB<float>* spatialDescriptorFunctionB;
    spatialDescriptorFunctionB = new SpatialDescriptorFunctionB<float>();
    spatialDescriptorFunctionB->setTriMesh( nucleusTriMesh );
    spatialDescriptor = spatialDescriptorFunctionB;
  }
  else if ( function == "C" )
  {
    PRINT("C");
    SpatialDescriptorFunctionC<float>* spatialDescriptorFunctionC;
    spatialDescriptorFunctionC = new SpatialDescriptorFunctionC<float>();
    spatialDescriptorFunctionC->setCenter( nucleusTriMesh.cog() );
    spatialDescriptor = spatialDescriptorFunctionC;
  }
  else //if ( function == "F" )
  {
    PRINT("F");
    SpatialModelCompleteRandomness3D<float> tempTriMeshSpatialModel;
    tempTriMeshSpatialModel.setRandomGenerator( randomGenerator );
    tempTriMeshSpatialModel.setTriMesh( nucleusTriMesh );
    tempTriMeshSpatialModel.initialize();
    Vertices<float> evaluationPositions = tempTriMeshSpatialModel.drawSample( 10000 );
    //evaluationPositions.save( parentDir + "/spatial_models/" + filename + "_Fpattern-" + ".vx", true);
    SpatialDescriptorFunctionF<float>* spatialDescriptorFunctionF;
    spatialDescriptorFunctionF = new SpatialDescriptorFunctionF<float>();
    spatialDescriptorFunctionF->setEvaluationPositions( evaluationPositions );
    spatialDescriptor = spatialDescriptorFunctionF;
  }


  if ( function != "all" )
  {
    modelEvaluator.setDescriptor( *spatialDescriptor );

    Vertices<float> vertices ( 3, numPoints, 0, 0 );
    for ( int i = 0; i < numPoints; ++i )
    {
      vertices[i][0] = datasetNucleus.getValue<float>( "centroidCoordX", i );
      vertices[i][1] = datasetNucleus.getValue<float>( "centroidCoordY", i );
      vertices[i][2] = datasetNucleus.getValue<float>( "centroidCoordZ", i );
      EVAL(vertices[i]);
    }


    ostringstream iss; //we suppose as much 99 labels
    iss << constraints;

    float sdi = modelEvaluator.eval( vertices, &saveTest );

//    Vector<float> output = modelEvaluator.evalSDIandMaxDiff( vertices, &saveTest );
//    EVAL( output[0] );
//    EVAL( output[1] );


   // saveTest.save( analysisDir + iss.str() + "/" + function + "/" + filename + ".csv", true );
  //  saveTest.save( analysisDir + iss.str() + "/" + function + "/" + filename + "_random.csv", true );
    const int row = dataSet.size()[0];
    EVAL(row);
    dataSet.setValue( "nucleus", row, filename );
    dataSet.setValue( "class", row, classif );//classification: mutant, tissue, etc.
    dataSet.setValue( "descriptor", row, function );//spatial descriptor
    dataSet.setValue( "index", row, sdi );
    //dataSet.setValue( "index", row, output[0] );
    //dataSet.setValue( "signedMaxDiff", row, output[1] );



    EVAL( copyGlobalAnalysis.columnType(0) );

//    for ( int i = 0; i < copyGlobalAnalysis.numRows(); ++i )
//    {
//      int tempNames = copyGlobalAnalysis.getValue<int>( copyGlobalAnalysis.variableNames()[0], i );
//      nucleiNames[i] = StringTools::toString( tempNames );
//    EVAL(nucleiNames[i]);
//    }

    for ( int i = 0; i < copyGlobalAnalysis.numRows(); ++i )
    {
      nucleiNames[i] = StringTools::toString(
            copyGlobalAnalysis.getValue<int>( copyGlobalAnalysis.variableNames()[0], i ) );
      EVAL(nucleiNames[i]);
    }

//    if ( copyGlobalAnalysis.columnType(0) == "%s" )
//      nucleiNames = copyGlobalAnalysis.getValues<string>( copyGlobalAnalysis.variableNames()[0] );
//    else if ( copyGlobalAnalysis.columnType(0) == "%d" )
//    {
//      for ( int i = 0; i < copyGlobalAnalysis.numRows(); ++i )
//      {
//        int tempNames = copyGlobalAnalysis.getValue<int>( copyGlobalAnalysis.variableNames()[0], i );
//        ostringstream tempName;
//        tempName << tempNames;
//        nucleiNames[i] = StringTools::toString( tempNames );
//        EVAL(nucleiNames[i]);
////        EVAL(tempNames);
//        EVAL(filename);
//      }
//    }
//    else if ( copyGlobalAnalysis.columnType(0) == "%f" )
//    {
//      for ( int i = 0; i < copyGlobalAnalysis.numRows(); ++i )
//      {
//        float tempNames = copyGlobalAnalysis.getValue<float>( copyGlobalAnalysis.variableNames()[0], i );
//        ostringstream tempName;
//        tempName << tempNames;
//        nucleiNames[i] = tempName.str();
//        EVAL(nucleiNames[i]);
//        EVAL(tempNames);
//        EVAL(filename);
//      }
//    }

    int numCurrentNucleus = 0;

    EVAL(filename);
    for ( int j = 0; j < nucleiNames.getSize(); ++j )
      if ( nucleiNames[j] == filename )
      {
        EVAL( nucleiNames[j] )
        numCurrentNucleus = j;
      }

    EVAL (numCurrentNucleus);
    string na = "NA";

    if ( ( sdi <= 1 ) && ( numCurrentNucleus != -1 ) )
      copyGlobalAnalysis.setValue( function, numCurrentNucleus, sdi );
    else
      copyGlobalAnalysis.setValue( function, numCurrentNucleus, na );

  }
  else
  {
    Vertices<float> vertices ( 3, numPoints, 0, 0 );
    for ( int i = 0; i < numPoints; ++i )
    {
      vertices[i][0] = datasetNucleus.getValue<float>( "centroidCoordX", i );
      vertices[i][1] = datasetNucleus.getValue<float>( "centroidCoordY", i );
      vertices[i][2] = datasetNucleus.getValue<float>( "centroidCoordZ", i );
      EVAL(vertices[i]);
    }

    const int row = dataSet.numRows();
    vector<float> sdis;
    vector<int> ranks;
//    vector<float> maxDiff;

//    modelEvaluator.evalSDIandMaxDiff( vertices, pValues, ranks, maxDiff);
    modelEvaluator.eval( vertices, sdis, ranks);

    const Vector<string>& nucleiNames = globalAnalysis.getValues<string>( "name" );
    int numCurrentNucleus = 0;

    for ( int j = 0; j < nucleiNames.getSize(); ++j )
      if ( nucleiNames[j] == filename )
        numCurrentNucleus = j;

    EVAL(numCurrentNucleus);

    string na = "NA";
    EVAL(sdis.size());

    for ( int jj = 0; jj < sdis.size(); ++jj )
    {
      if ( ( sdis[jj] < 0 || sdis[jj] > 1 )  && ( jj = 0 ) )
        copyGlobalAnalysis.setValue( "F", numCurrentNucleus, na );
      else
        copyGlobalAnalysis.setValue( "F", numCurrentNucleus, sdis[jj] );
      if ( ( sdis[jj] < 0 || sdis[jj] > 1 )  && ( jj = 1 ) )
        copyGlobalAnalysis.setValue( "G", numCurrentNucleus, na );
      else
        copyGlobalAnalysis.setValue( "G", numCurrentNucleus, sdis[jj] );
      if ( ( sdis[jj] < 0 || sdis[jj] > 1 )  && ( jj = 2 ) )
        copyGlobalAnalysis.setValue( "H", numCurrentNucleus, na );
      else
        copyGlobalAnalysis.setValue( "H", numCurrentNucleus, sdis[jj] );
      if ( ( sdis[jj] < 0 || sdis[jj] > 1 )  && ( jj = 3 ) )
        copyGlobalAnalysis.setValue( "B", numCurrentNucleus, na );
      else
        copyGlobalAnalysis.setValue( "B", numCurrentNucleus, sdis[jj] );
      if ( ( sdis[jj] < 0 || sdis[jj] > 1 )  && ( jj = 4 ) )
        copyGlobalAnalysis.setValue( "C", numCurrentNucleus, na );
      else
        copyGlobalAnalysis.setValue( "C", numCurrentNucleus, sdis[jj] );
    }

    dataSet.setValue( "nucleus", row, filename );
    dataSet.setValue( "class", row, classif );//classification: mutant, tissue, etc.
    dataSet.setValue( "F-SDI", row, sdis[0] );
//    dataSet.setValue( "F-maxDiff", row, maxDiff[0] );
    dataSet.setValue( "G-SDI", row, sdis[1] );
//    dataSet.setValue( "G-maxDiff", row, maxDiff[1] );
    dataSet.setValue( "H-SDI", row, sdis[2] );
//    dataSet.setValue( "H-maxDiff", row, maxDiff[2] );
    dataSet.setValue( "B-SDI", row, sdis[3] );
//    dataSet.setValue( "B-maxDiff", row, maxDiff[3] );
    dataSet.setValue( "C-SDI", row, sdis[4] );
//    dataSet.setValue( "C-maxDiff", row, maxDiff[4] );

    saveTest.setValue( "nucleus", 1, filename );
    saveTest.setValue( "class", 1, classif );//classification: mutant, tissue, etc.
    saveTest.setValue( "F-SDI", 1, sdis[0] );
//    saveTest.setValue( "F-maxDiff", 1, maxDiff[0] );
    saveTest.setValue( "G-SDI", 1, sdis[1] );
//    saveTest.setValue( "G-maxDiff", 1, maxDiff[1] );
    saveTest.setValue( "H-SDI", 1, sdis[2] );
//    saveTest.setValue( "H-maxDiff", 1, maxDiff[2] );
    saveTest.setValue( "B-SDI", 1, sdis[3] );
//    saveTest.setValue( "B-maxDiff", 1, maxDiff[3] );
    saveTest.setValue( "C-SDI", 1, sdis[4] );
//    saveTest.setValue( "C-maxDiff", 1, maxDiff[4] );

    ostringstream iss; //we suppose as much 99 labels
    iss << constraints;

    saveTest.save( analysisDir + iss.str() + "/" + filename + "_all.csv", true );
    copyGlobalAnalysis.save( analysisDir + "/" + "nuclei_complete.csv", true );

  }



}


/*! Preparing constraints for the model
****************************************************************/
void evaluator_completeSpatialRandomness(
  const string& filename, const string& parentDir,
  const string& function, const int& constraints,
  DataSet& dataSet,
  RandomGenerator& randomGenerator)
{
  PRINT("spatialModelEvaluator_completeSpatialRandomness");

  //const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + "_nucleus.tm" );
  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + ".tm" );
  SpatialModelCompleteRandomness3D<float> triMeshSpatialModel;
  triMeshSpatialModel.setRandomGenerator( randomGenerator );
  triMeshSpatialModel.setTriMesh( nucleusTriMesh );
  triMeshSpatialModel.initialize();

  evaluator(
    nucleusTriMesh,
    triMeshSpatialModel,
    filename, parentDir,
    function, constraints,
    dataSet, randomGenerator );
}

void evaluator_sizeConstrained(
  const string& filename, const string& parentDir,
  const string& function, const int& constraints,
  DataSet& dataSet,
  RandomGenerator& randomGenerator)
{
  PRINT("spatialModelEvaluator_sizeConstrained");

  const string analysisDir = parentDir + "/analysis/";

  //new data
  //const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + "_nucleus.tm" );
  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + ".tm" );
  //const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/nuclei/" + filename + ".tm" );
  //old data
  //const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + "-nucleus.tm" );

  //new data
  const DataSet datasetNucleus( analysisDir + filename + "_chromocenters.csv" );
  //const DataSet datasetNucleus( analysisDir + filename + "_nucleoli.csv" );
  //old data
  //const DataSet datasetNucleus( analysisDir + filename + ".csv" );

  //new data
  const Vector<float> eqRadii = datasetNucleus.getValues<float>( "equivalentRadius_tm" );
  //old data
  //const Vector<float> eqRadii = datasetNucleus.getValues<float>( "chromocenterRadius" );
  EVAL(eqRadii);

  SpatialModelHardcoreDistance3D<float> triMeshSpatialModel;
  triMeshSpatialModel.setRandomGenerator( randomGenerator );
  triMeshSpatialModel.setTriMesh( nucleusTriMesh );
  triMeshSpatialModel.setHardcoreDistances( eqRadii );
  triMeshSpatialModel.initialize();

  evaluator(
    nucleusTriMesh,
    triMeshSpatialModel,
    filename, parentDir,
    function, constraints,
    dataSet, randomGenerator );
}

void evaluator_distanceConstrained(
  const string& filename, const string& parentDir,
  const string& function, const int& constraints,
  DataSet& dataSet,
  RandomGenerator& randomGenerator)
{
  PRINT("spatialModelEvaluator_distanceConstrained");

  //const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + "_nucleus.tm" );
  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + ".tm" );
  const string analysisDir = parentDir + "/analysis/";
  const DataSet datasetNucleus( analysisDir + filename + "_chromocenters.csv" );
  //const DataSet datasetNucleus( analysisDir + filename + ".csv" );
  const Vector<float> distancesToBorder = datasetNucleus.getValues<float>( "distanceToTheBorder" );

  SpatialModelBorderDistance3D<float> triMeshSpatialModel;
  triMeshSpatialModel.setRandomGenerator( randomGenerator );
  triMeshSpatialModel.setTriMesh( nucleusTriMesh );
  triMeshSpatialModel.setDistancesToBorder( distancesToBorder );
  triMeshSpatialModel.initialize();

  evaluator(
    nucleusTriMesh,
    triMeshSpatialModel,
    filename, parentDir,
    function, constraints,
    dataSet, randomGenerator );
}

void evaluator_sizeAndDistanceConstrained(
  const string& filename, const string& parentDir,
  const string& function, const int& constraints,
  DataSet& dataSet,
  RandomGenerator& randomGenerator)
{
  PRINT("spatialModelEvaluator_sizeAndDistanceConstrained");

  //randomGenerator.init(11);
  const string analysisDir = parentDir + "/analysis/";

  //new data
  //const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + "_nucleus.tm" );
  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + ".tm" );
  //old data
  //const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + "-nucleus.tm" );

  //new data
  const DataSet datasetNucleus( analysisDir + filename + "_chromocenters.csv" );
  //old data
  //const DataSet datasetNucleus( analysisDir + filename + ".csv" );

  //new data
  const Vector<float> eqRadii = datasetNucleus.getValues<float>( "equivalentRadius_tm" );
  //old data
  //const Vector<float> eqRadii = datasetNucleus.getValues<float>( "chromocenterRadius" );
  EVAL(eqRadii);
  const Vector<float> distancesToBorder = datasetNucleus.getValues<float>( "distanceToTheBorder" );

  SpatialModelHardcoreBorderDistance3D<float> triMeshSpatialModel;
  triMeshSpatialModel.setRandomGenerator( randomGenerator );
  triMeshSpatialModel.setTriMesh( nucleusTriMesh );
  triMeshSpatialModel.setDistancesToBorder( distancesToBorder );
  triMeshSpatialModel.setHardcoreDistances( eqRadii );
  triMeshSpatialModel.initialize();

  evaluator(
    nucleusTriMesh, triMeshSpatialModel, filename, parentDir,
    function, constraints, dataSet, randomGenerator );
}

void evaluator_MaximalRepulsionConstrained(
  const string& filename, const string& parentDir,
  const string& function, const int& constraints,
  DataSet& dataSet,
  RandomGenerator& randomGenerator)
{
  PRINT("spatialModelEvaluator_MaximalRepulsionConstrained");

  const string analysisDir = parentDir + "/analysis/";

  //new data
  //const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + "_nucleus.tm" );
  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + ".tm" );
  //old data
  //const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + "-nucleus.tm" );

  //new data
  const DataSet datasetNucleus( analysisDir + filename + "_chromocenters.csv" );
  //old data
  //const DataSet datasetNucleus( analysisDir + filename + ".csv" );

  //new data
  const Vector<float> eqRadii = datasetNucleus.getValues<float>( "equivalentRadius_tm" );
  //old data
  //const Vector<float> eqRadii = datasetNucleus.getValues<float>( "chromocenterRadius" );
  EVAL(eqRadii);

  SpatialModelMaximalRepulsion3D<float> triMeshSpatialModel;
  triMeshSpatialModel.setRandomGenerator( randomGenerator );
  triMeshSpatialModel.setTriMesh( nucleusTriMesh );
  triMeshSpatialModel.setHardcoreDistances( eqRadii );
  //triMeshSpatialModel.setNumMonteCarloCycles( 1000 );
  triMeshSpatialModel.initialize();

  evaluator(
    nucleusTriMesh, triMeshSpatialModel, filename, parentDir,
    function, constraints, dataSet, randomGenerator );
}


//void evaluator_MaxRepulsionWithDistanceToTheBorderConstrained(
//  const string& filename, const string& parentDir,
//  const string& function, const int& constraints,
//  DataSet& dataSet,
//  RandomGenerator& randomGenerator)
//{
//  PRINT("spatialModelEvaluator_MaxRepulsionWithDistanceToTheBorderConstrained");

//  const string analysisDir = parentDir + "/analysis/";

//  //new data
//  //const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + "_nucleus.tm" );
//  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + ".tm" );
//  //old data
//  //const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + "-nucleus.tm" );

//  //new data
//  const DataSet datasetNucleus( analysisDir + filename + "_chromocenters.csv" );
//  //old data
//  //const DataSet datasetNucleus( analysisDir + filename + ".csv" );
//  const Vector<float> distancesToBorder = datasetNucleus.getValues<float>( "distanceToTheBorder" );

//  //new data
//  const Vector<float> eqRadii = datasetNucleus.getValues<float>( "equivalentRadius_tm" );
//  //old data
//  //const Vector<float> eqRadii = datasetNucleus.getValues<float>( "chromocenterRadius" );
//  EVAL(eqRadii);

//  MaxRepulsionWithDistancesTriMeshSpatialModel<float> triMeshSpatialModel;
//  triMeshSpatialModel.setRandomGenerator( randomGenerator );
//  triMeshSpatialModel.setTriMesh( nucleusTriMesh );
//  triMeshSpatialModel.setDistanceToBorder( distancesToBorder );
//  triMeshSpatialModel.setHardcoreDistance( eqRadii );
//  triMeshSpatialModel.initialize();

//  evaluator(
//    nucleusTriMesh, triMeshSpatialModel, filename, parentDir,
//    function, constraints, dataSet, randomGenerator );
//}


/*! Chooses case depending on the constraints
****************************************************************/
void realDataEvaluator(
  const string& filename, const string& parentDir,
  const string& function, const int& constraints,
  DataSet& dataSet,
  RandomGenerator& randomGenerator)
{
  switch( constraints )
  {
    case 0:
      evaluator_completeSpatialRandomness(
        filename, parentDir,
        function, constraints,
        dataSet, randomGenerator );
      break;
    case 1:
      evaluator_sizeConstrained(
        filename, parentDir,
        function, constraints,
        dataSet, randomGenerator );
      break;
    case 2:
      evaluator_distanceConstrained(
        filename, parentDir,
        function, constraints,
        dataSet, randomGenerator );
      break;
    case 3:
      evaluator_sizeAndDistanceConstrained(
        filename, parentDir,
        function, constraints,
        dataSet, randomGenerator );
      break;
    case 4:
      evaluator_MaximalRepulsionConstrained(
        filename, parentDir,
        function, constraints,
        dataSet, randomGenerator );
      break;
//    case 5:
//    evaluator_MaxRepulsionWithDistanceToTheBorderConstrained(
//      filename, parentDir,
//      function, constraints,
//      dataSet, randomGenerator );
//    break;
  }
}
