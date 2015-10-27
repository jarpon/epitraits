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


void evaluator(
  const TriMesh<float>& nucleusTriMesh,
  TriMeshSpatialModel<float>& triMeshSpatialModel,
  const string& filename, const string& parentDir,
  const string& function, const int& constraints,
  DataSet& dataSet, RandomGenerator& randomGenerator)
{
  string classif = parentDir;
  classif = classif.substr( classif.find_last_of("/\\")+1, classif.length() );

  //open data info
  const string analysisDir = parentDir + "/analysis/";
  const DataSet datasetNucleus( analysisDir + filename + "_chromocenters.csv" );
//  const DataSet datasetNucleus( analysisDir + filename + "_nucleoli.csv" );
  DataSet globalAnalysis( analysisDir + "nuclei.csv" );
  const DataSet ccsInfo( analysisDir + "ccs.csv" );

  Vertices<float> vertices;
  Vector<string> tempFileNames;

  tempFileNames = ccsInfo.getValues<string>( ccsInfo.variableNames()[0] );

  Vector<int> posCCSlist, indCCSlist;
  indCCSlist.setSize( 1 );

  for ( int j = 0; j < tempFileNames.getSize(); ++j )
    if ( tempFileNames[j] == filename )
    {
      EVAL(j);

      indCCSlist[0] = j;
      posCCSlist.append( indCCSlist );
    }
  EVAL( posCCSlist.getSize() );
  EVAL (posCCSlist);

  for ( int j = posCCSlist[0]; j < posCCSlist[posCCSlist.getSize()-1]; ++j )
  {
    EVAL(j);
    vertices[j][0] = ccsInfo.getValue<float>( "centroidCoordX", j );
    vertices[j][1] = ccsInfo.getValue<float>( "centroidCoordY", j );
    vertices[j][2] = ccsInfo.getValue<float>( "centroidCoordZ", j );
    EVAL(vertices[j]);
  }

  Vector<string> nucleiNames ;
  nucleiNames = globalAnalysis.getValues<string>( globalAnalysis.variableNames()[0] );

  //unifying datasets
  int numCurrentNucleus;

  for ( int j = 0; j < nucleiNames.getSize(); ++j )
    if ( nucleiNames[j] == filename )
      numCurrentNucleus = j;

//  EVAL (numCurrentNucleus);

  //const int numPoints = vertices.getNumVertices();
  const int numPatterns = 99;

  SpatialModelEvaluator<float,float> modelEvaluator;
  modelEvaluator.setModel( triMeshSpatialModel );
  modelEvaluator.setNumMonteCarloSamples( numPatterns ); //to check uniformity
  modelEvaluator.setPrecision( 0.05 );

  SpatialDescriptor<float>* spatialDescriptor;

  DataSet saveTest;

  //setting function parameters
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
    int row;

    modelEvaluator.setDescriptor( *spatialDescriptor );

//    Vertices<float> vertices ( 3, numPoints, 0, 0 );
//    for ( int i = 0; i < numPoints; ++i )
//    {
//      vertices[i][0] = datasetNucleus.getValue<float>( "centroidCoordX", i );
//      vertices[i][1] = datasetNucleus.getValue<float>( "centroidCoordY", i );
//      vertices[i][2] = datasetNucleus.getValue<float>( "centroidCoordZ", i );
//      EVAL(vertices[i]);
//    }


    ostringstream iss; //we suppose as much 99 labels
    iss << constraints;

    try {
      sdi = modelEvaluator.eval( vertices, &saveTest );

  //    Vector<float> output = modelEvaluator.evalSDIandMaxDiff( vertices, &saveTest );
  //    EVAL( output[0] );
  //    EVAL( output[1] );


     // saveTest.save( analysisDir + iss.str() + "/" + function + "/" + filename + ".csv", true );
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
    catch( Exception() ) {
      const int row = dataSet.size()[0];

      dataSet.setValue( "nucleus", row, filename );
      dataSet.setValue( "class", row, classif );//classification: mutant, tissue, etc.
      dataSet.setValue( "descriptor", row, function );//spatial descriptor
      dataSet.setValue( "index", row, sqrt(-1) );
      //dataSet.setValue( "index", row, output[0] );
      //dataSet.setValue( "signedMaxDiff", row, output[1] );

      globalAnalysis.setValue( newVariable + "G-SDI", numCurrentNucleus, sqrt(-1) );
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
//    vector<float> maxDiff;

//    modelEvaluator.evalSDIandMaxDiff( vertices, pValues, ranks, maxDiff);
    modelEvaluator.eval( vertices, sdis, ranks);

    int numCurrentNucleus;

    for ( int j = 0; j < nucleiNames.getSize(); ++j )
      if ( nucleiNames[j] == filename )
        numCurrentNucleus = j;

    EVAL(numCurrentNucleus);


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
    globalAnalysis.save( analysisDir + "/" + "nuclei_complete.csv", true );

  }

  globalAnalysis.save( analysisDir + "/" + "nuclei.csv", true );


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

  //open data info
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

  //open data info
  const string analysisDir = parentDir + "/analysis/";
  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + ".tm" );
  const DataSet datasetNucleus( analysisDir + filename + "_chromocenters.csv" );
  const Vector<float> eqRadii = datasetNucleus.getValues<float>( "equivalentRadius_tm" );
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

  //open data info
  const string analysisDir = parentDir + "/analysis/";
  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + ".tm" );
  const DataSet datasetNucleus( analysisDir + filename + "_chromocenters.csv" );
  const Vector<float> distancesToBorder = datasetNucleus.getValues<float>( "distanceToTheBorder" );
  EVAL(distancesToBorder);

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

  //open data info
  const string analysisDir = parentDir + "/analysis/";
  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + ".tm" );
  const DataSet datasetNucleus( analysisDir + filename + "_chromocenters.csv" );
  const Vector<float> eqRadii = datasetNucleus.getValues<float>( "equivalentRadius_tm" );
  const Vector<float> distancesToBorder = datasetNucleus.getValues<float>( "distanceToTheBorder" );
  EVAL(eqRadii);
  EVAL(distancesToBorder);

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

  //open data info
  const string analysisDir = parentDir + "/analysis/";
  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + ".tm" );
  const DataSet datasetNucleus( analysisDir + filename + "_chromocenters.csv" );
  const Vector<float> eqRadii = datasetNucleus.getValues<float>( "equivalentRadius_tm" );
  EVAL(eqRadii);

  SpatialModelMaximalRepulsion3D<float> triMeshSpatialModel;
  triMeshSpatialModel.setRandomGenerator( randomGenerator );
  triMeshSpatialModel.setTriMesh( nucleusTriMesh );
  //triMeshSpatialModel.setNumMonteCarloCycles( 10 );
  triMeshSpatialModel.setHardcoreDistances( eqRadii );
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
