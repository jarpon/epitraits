#include <spatialmodelevaluatorbase.h>
//#include "spatialmodelevaluator2.h"
#include <spatialdescriptorfunctionf.h>
#include <spatialdescriptorfunctiong.h>
#include <spatialdescriptorfunctionh.h>
#include <spatialdescriptorfunctionb.h>
#include <spatialdescriptorfunctionc.h>
#include "spatialdescriptorfunctionz.h"
#include <spatialmodelmaximalrepulsion3d.h>
#include "maxrepulsionwithdistances.h"
#include <spatialmodelcompleterandomness3d.h>
#include <spatialmodelborderdistance3d.h>
#include <spatialmodelhardcoreborderdistance3d.h>
#include <spatialmodelhardcoredistance3d.h>
#include <spatialmodel.h>
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
#include <programerror.h>
#include <vertexstack.h>

#include <stringtools.h>
#include <fstream>
#include <sstream>
#include <unistd.h>
 #include <cstdlib>

#define TRACE
#include <trace.h>
using namespace std;


void evaluator(
  const string& filename, const string& parentDir, const string& patternsDir,
  const string& function, const int& constraints, DataSet& dataSet, RandomGenerator& randomGenerator )
{
  ENTER("Analyzing data using existing spatial patterns");
  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + ".tm" );

  string classif = parentDir;
  classif = classif.substr( classif.find_last_of("/\\")+1, classif.length() );
  //open data info
  const string analysisDir = parentDir + "/analysis/";
//  const DataSet datasetNucleus( analysisDir + filename + "_chromocenters.csv" );
//  const DataSet datasetNucleus( analysisDir + filename + "_nucleoli.csv" );
  //DataSet globalAnalysis( analysisDir + "nuclei.csv" );
  DataSet globalAnalysis( analysisDir + "nuclei_extended.csv" );
  const DataSet ccsInfo( analysisDir + "ccs.csv" );

  Vector<string> tempFileNames;
  tempFileNames = ccsInfo.getValues<string>( ccsInfo.variableNames()[0] );

  int lastPos, numCCS = 0;

  for ( int j = 0; j < tempFileNames.getSize(); ++j )
    if ( tempFileNames[j] == filename )
    {
      lastPos = j;
      ++ numCCS;
    }

  if ( numCCS == 0 )
  {
    ProgramError error;
    error.setWhat( "There is a problem with the chromocenters dataset" );
    throw error;
    return;
  }

  VertexStack<float> spatialPatterns;
  spatialPatterns.load( patternsDir + filename + ".vs" );
  EVAL(spatialPatterns.getSize());

  ShapeSet<float> modelPatterns, comparisonPatterns;

  const int numTotalPatterns = spatialPatterns.getSize();
  Vertices<float> temp( 3, numCCS, 0, 0 );
  for ( int i = 0; i < numTotalPatterns/2; ++i )
  {
//    temp = spatialPatterns[i];
//    modelPatterns.insert( i, temp );
    modelPatterns.addShape( &spatialPatterns[i] );
  }
  for ( int j = numTotalPatterns/2; j < numTotalPatterns; ++j )
  {
//    temp = spatialPatterns[j];
//    comparisonPatterns.insert( j, temp );
    comparisonPatterns.addShape( &spatialPatterns[j] );
  }

  Vertices<float> vertices( 3, numCCS, 0, 0 );
  int k = 0;
  for ( int j = lastPos - numCCS + 1 ; j < lastPos + 1; ++j, ++k )
  {
    vertices[k][0] = ccsInfo.getValue<float>( "centroidCoordX", j );
    vertices[k][1] = ccsInfo.getValue<float>( "centroidCoordY", j );
    vertices[k][2] = ccsInfo.getValue<float>( "centroidCoordZ", j );
    EVAL(vertices[k]);
  }

  Vector<string> nucleiNames ;
  nucleiNames = globalAnalysis.getValues<string>( globalAnalysis.variableNames()[0] );

  //unifying datasets
  int numCurrentNucleus;

  for ( int j = 0; j < nucleiNames.getSize(); ++j )
    if ( nucleiNames[j] == filename )
      numCurrentNucleus = j;

  EVAL (numCurrentNucleus);

  //const int numPoints = vertices.getNumVertices();
  const int numPatterns = 99;

  SpatialModelEvaluatorBase<float,float> modelEvaluator;
  //modelEvaluator.setModel( triMeshSpatialModel );
  //modelEvaluator.setNumMonteCarloSamples( numPatterns ); //to check uniformity
  modelEvaluator.setPrecision( 0.05 );

  SpatialDescriptor<float>* spatialDescriptor;


  //setting function parameters
  if ( function == "all" )
  {
    PRINT("all functions");

    SpatialModelCompleteRandomness3D<float> tempTriMeshSpatialModel;
    tempTriMeshSpatialModel.setRandomGenerator( randomGenerator );
    tempTriMeshSpatialModel.setTriMesh( nucleusTriMesh );
    tempTriMeshSpatialModel.initialize();
    Vertices<float> evaluationPositions = tempTriMeshSpatialModel.drawSample( 10000 );
    //evaluationPositions.save( parentDir + "/" + filename + "_Fpattern-" + ".vx", true); //to check uniformity of the F patterns

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

    spatialDescriptor = new SpatialDescriptorFunctionZ<float>();
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
  else if ( function == "Z" )
  {
    PRINT("Z");
    spatialDescriptor = new SpatialDescriptorFunctionZ<float>();
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


  //processing data
  if ( function != "all" )
  {
    string newVariable;
    switch ( constraints )
    {
      case 0:
        newVariable = ( "SpatialModelCompleteRandomness3D_" + function + "-SDI" );
        break;
      case 1:
        newVariable = ( "SpatialModelHardcoreDistance3D_" + function + "-SDI" );
        break;
      case 2:
        newVariable = ( "spatialModelBorderDistance3D_" + function + "-SDI" );
        break;
      case 3:
        newVariable = ( "SpatialModelBorderHardcoreDistance3D_" + function + "-SDI" );
        break;
      case 4:
        newVariable = ( "SpatialModelMaximalRepulsion3D_" + function + "-SDI" );
        break;
    }

    //globalAnalysis.setValues<float>( newVariable, -1 );
    float sdi;
    int row, rank;
    vector<float> sdis;
    vector<int> ranks;

    modelEvaluator.setDescriptor( *spatialDescriptor );

    ostringstream iss; //we suppose as much 99 labels
    iss << constraints;
    DataSet saveTest;
    try {
      //sdi = modelEvaluator.eval( vertices, &saveTest );
      modelEvaluator.eval( vertices, modelPatterns, comparisonPatterns, sdis, ranks, &saveTest );

  //    Vector<float> output = modelEvaluator.evalSDIandMaxDiff( vertices, &saveTest );
  //    EVAL( output[0] );
  //    EVAL( output[1] );


      saveTest.save( analysisDir + iss.str() + "/" + function + "/" + filename + ".csv", true );
    //  saveTest.save( analysisDir + iss.str() + "/" + function + "/" + filename + "_random.csv", true );
      row = dataSet.size()[0];

      dataSet.setValue( "nucleus", row, filename );
      dataSet.setValue( "class", row, classif );//classification: mutant, tissue, etc.
      dataSet.setValue( "descriptor", row, function );//spatial descriptor
      dataSet.setValue( "index", row, sdi );
      //dataSet.setValue( "index", row, output[0] );
      //dataSet.setValue( "signedMaxDiff", row, output[1] );

      globalAnalysis.setValue( newVariable, numCurrentNucleus, sdi );
    }
    catch( Exception exception ) {
      EVAL( exception.getWhat() );
      EVAL( exception.getWhere() );

      const int row = dataSet.size()[0];

      dataSet.setValue( "nucleus", row, filename );
      dataSet.setValue( "class", row, classif );//classification: mutant, tissue, etc.
      dataSet.setValue( "descriptor", row, function );//spatial descriptor
      dataSet.setValue( "index", row, sqrt(-1) );


    }

  }
  else
  {
//    Vertices<float> vertices ( 3, numPoints, 0, 0 );
//    for ( int i = 0; i < numPoints; ++i )
//    {
//      vertices[i][0] = datasetNucleus.getValue<float>( "centroidCoordX", i );
//      vertices[i][1] = datasetNucleus.getValue<float>( "centroidCoordY", i );
//      vertices[i][2] = datasetNucleus.getValue<float>( "centroidCoordZ", i );
//      EVAL(vertices[i]);
//    }

    const int row = dataSet.numRows();
    vector<float> sdis;
    vector<int> ranks;
    DataSet saveTest;

//    vector<float> maxDiff;

    try {
//    modelEvaluator.evalSDIandMaxDiff( vertices, pValues, ranks, maxDiff);

   // modelEvaluator.eval( vertices, sdis, ranks, &saveTest );

    modelEvaluator.eval( vertices, modelPatterns, comparisonPatterns, sdis, ranks, &saveTest );

    EVAL( sdis[0] );
    EVAL( sdis[1] );
    EVAL( sdis[2] );
    EVAL( sdis[3] );
    EVAL( sdis[4] );
    EVAL( sdis[5] );

    string newVariable;
    switch ( constraints )
    {
      case 0: newVariable = "SpatialModelCompleteRandomness3D_";
      case 1: newVariable = "SpatialModelHardcoreDistance3D_";
      case 2: newVariable = "spatialModelBorderDistance3D_";
      case 3: newVariable = "SpatialModelBorderHardcoreDistance3D_";
      case 4: newVariable = "SpatialModelMaximalRepulsion3D_";
    }

    //unifying datasets
    for ( int jj = 0; jj < sdis.size(); ++jj )
    {
      if ( ( sdis[jj] < 0 || sdis[jj] > 1 )  && ( jj = 0 ) )
        globalAnalysis.setValue( newVariable + "F-SDI", numCurrentNucleus, sqrt(-1) );
      else
        globalAnalysis.setValue( newVariable + "F-SDI", numCurrentNucleus, sdis[jj] );
      if ( ( sdis[jj] < 0 || sdis[jj] > 1 )  && ( jj = 1 ) )
        globalAnalysis.setValue( newVariable + "G-SDI", numCurrentNucleus, sqrt(-1) );
      else
        globalAnalysis.setValue( newVariable + "G-SDI", numCurrentNucleus, sdis[jj] );
      if ( ( sdis[jj] < 0 || sdis[jj] > 1 )  && ( jj = 2 ) )
        globalAnalysis.setValue( newVariable + "H-SDI", numCurrentNucleus, sqrt(-1) );
      else
        globalAnalysis.setValue( newVariable + "H-SDI", numCurrentNucleus, sdis[jj] );
      if ( ( sdis[jj] < 0 || sdis[jj] > 1 )  && ( jj = 3 ) )
        globalAnalysis.setValue( newVariable + "B-SDI", numCurrentNucleus, sqrt(-1) );
      else
        globalAnalysis.setValue( newVariable + "B-SDI", numCurrentNucleus, sdis[jj] );
      if ( ( sdis[jj] < 0 || sdis[jj] > 1 )  && ( jj = 4 ) )
        globalAnalysis.setValue( newVariable + "C-SDI", numCurrentNucleus, sqrt(-1) );
      else
        globalAnalysis.setValue( newVariable + "C-SDI", numCurrentNucleus, sdis[jj] );
      if ( ( sdis[jj] < 0 || sdis[jj] > 1 )  && ( jj = 5 ) )
        globalAnalysis.setValue( newVariable + "Z-SDI", numCurrentNucleus, sqrt(-1) );
      else
        globalAnalysis.setValue( newVariable + "Z-SDI", numCurrentNucleus, sdis[jj] );
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
    dataSet.setValue( "Z-SDI", row, sdis[4] );
//    dataSet.setValue( "Z-maxDiff", row, maxDiff[5] );

    ostringstream iss; //we have 4 constraints
    iss << constraints;
    saveTest.save( analysisDir + iss.str() + "/" + filename + ".csv", true );

    }
    catch( Exception exception ) {
      EVAL( exception.getWhat() );
      EVAL( exception.getWhere() );
      const int row = dataSet.size()[0];

      dataSet.setValue( "nucleus", row, filename );
      dataSet.setValue( "class", row, classif );//classification: mutant, tissue, etc.
      dataSet.setValue( "descriptor", row, function );//spatial descriptor
      dataSet.setValue( "index", row, sqrt(-1) );
      //dataSet.setValue( "index", row, output[0] );
      //dataSet.setValue( "signedMaxDiff", row, output[1] );

    }

  }

  globalAnalysis.save( analysisDir + "/" + "nuclei_complete2.csv", true );

  LEAVE();
}

/*! Chooses folder input of the patterns depending on the spatial model
****************************************************************/
void realDataEvaluatorExternalPatterns(
  const string& filename, const string& parentDir,
  const string& function, const int& constraints,
    DataSet& dataSet, RandomGenerator& randomGenerator )
{
  switch( constraints )
  {
    case 0:
      const string patternsDir = parentDir + "/patterns/SpatialModelCompleteRandomness3D/";
      evaluator(
        filename, parentDir, patternsDir,
            function, constraints, dataSet, randomGenerator );
        break;
//    case 1:
//    const string patternsDir = parentDir + "/patterns/SpatialModelHardcoreDistance3D/";
//    evaluator(
//      filename, parentDir, patternsDir,
//          function, constraints, dataSet, randomGenerator );
//      break;
//    case 2:
//    const string patternsDir = parentDir + "/patterns/SpatialModelBorderDistance3D/";
//    evaluator(
//      filename, parentDir, patternsDir,
//          function, dataSet, randomGenerator );
//      break;
//    case 3:
//    const string patternsDir = parentDir + "/patterns/SpatialModelHardcoreBorderDistance3D/";
//    evaluator(
//      filename, parentDir, patternsDir,
//          function, dataSet, randomGenerator );
//      break;
//    case 4:
//    const string patternsDir = parentDir + "/patterns/SpatialModelMaximalRepulsion3D/";
//    evaluator(
//      filename, parentDir, patternsDir,
//          function, dataSet, randomGenerator );
//      break;
//    case 5:
//    evaluator(
//      filename, parentDir, patternsDir,
//      function, dataSet );
//    break;
  }
}
