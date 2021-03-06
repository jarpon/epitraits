#include <spatialmodelevaluator.h>
//#include "spatialmodelevaluator2.h"
#include <spatialdescriptorfunctionf.h>
#include <spatialdescriptorfunctiong.h>
#include <spatialdescriptorfunctionh.h>
#include <spatialdescriptorfunctionb.h>
#include <spatialdescriptorfunctionc.h>
#include "spatialdescriptorfunctionz.h"
#include "spatialdescriptorfunctionlrd.h"
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
  const string& filename, const string& parentDir, const string& spatialModel,
  const string& function, const int& constraints,
  DataSet& dataSet, RandomGenerator& randomGenerator)
{
  string classif = parentDir;
  classif = classif.substr( classif.find_last_of("/\\")+1, classif.length() );
  //open data info
  const string analysisDir = parentDir + "/analysis/";
//  const DataSet datasetNucleus( analysisDir + filename + "_chromocenters.data" );
//  const DataSet datasetNucleus( analysisDir + filename + "_nucleoli.data" );
  //DataSet globalAnalysis( analysisDir + "nuclei.data" );
  DataSet globalAnalysis( analysisDir + "nuclei_extended.data" );
  const DataSet ccsInfo( analysisDir + "ccs.data" );

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
    EVAL("Nucleus not found");
    return;
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

  SpatialModelEvaluator<float,float> modelEvaluator;
  modelEvaluator.setModel( triMeshSpatialModel );
  modelEvaluator.setNumMonteCarloSamples( numPatterns ); //to check uniformity
  //modelEvaluator.setPrecision( 0.05 );
  modelEvaluator.setPrecision( 0.01 );

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

    spatialDescriptor = new SpatialDescriptorFunctionLRD<float>();
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
  else if ( function == "LRD" )
  {
    PRINT("LRD");
    spatialDescriptor = new SpatialDescriptorFunctionLRD<float>();
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
    //globalAnalysis.setValues<float>( newVariable, -1 );
    float sdi;
    int row;

    modelEvaluator.setDescriptor( *spatialDescriptor );

    ostringstream iss; //we suppose as much 99 labels
    iss << constraints;
    DataSet saveTest;
    try {
      sdi = modelEvaluator.eval( vertices, &saveTest );
      EVAL(sdi);

  //    Vector<float> output = modelEvaluator.evalSDIandMaxDiff( vertices, &saveTest );
  //    EVAL( output[0] );
  //    EVAL( output[1] );


      saveTest.save( analysisDir + "/" + spatialModel + "/" + function + "/" + filename + ".data", true );
    //  saveTest.save( analysisDir + iss.str() + "/" + function + "/" + filename + "_random.data", true );
      row = dataSet.size()[0];

      dataSet.setValue( "nucleus", row, filename );
      dataSet.setValue( "class", row, classif );//classification: mutant, tissue, etc.
      dataSet.setValue( "descriptor", row, function );//spatial descriptor
      dataSet.setValue( "index", row, sdi );
      //dataSet.setValue( "index", row, output[0] );
      //dataSet.setValue( "signedMaxDiff", row, output[1] );

      //globalAnalysis.setValue( newVariable, numCurrentNucleus, sdi );
      globalAnalysis.setValue( spatialModel + "_" + function + "-SDI", numCurrentNucleus, sdi );
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

    modelEvaluator.eval( vertices, sdis, ranks, &saveTest );

    EVAL( sdis[0] );
    EVAL( sdis[1] );
    EVAL( sdis[2] );
    EVAL( sdis[3] );
    EVAL( sdis[4] );
    EVAL( sdis[5] );
    EVAL( sdis[6] );

    //unifying datasets
    for ( int jj = 0; jj < sdis.size(); ++jj )
    {
      if ( ( sdis[jj] < 0 || sdis[jj] > 1 )  && ( jj = 0 ) )
        globalAnalysis.setValue( spatialModel + "_F-SDI", numCurrentNucleus, sqrt(-1) );
      else
        globalAnalysis.setValue( spatialModel + "_F-SDI", numCurrentNucleus, sdis[jj] );
      if ( ( sdis[jj] < 0 || sdis[jj] > 1 )  && ( jj = 1 ) )
        globalAnalysis.setValue( spatialModel + "_G-SDI", numCurrentNucleus, sqrt(-1) );
      else
        globalAnalysis.setValue( spatialModel + "_G-SDI", numCurrentNucleus, sdis[jj] );
      if ( ( sdis[jj] < 0 || sdis[jj] > 1 )  && ( jj = 2 ) )
        globalAnalysis.setValue( spatialModel + "_H-SDI", numCurrentNucleus, sqrt(-1) );
      else
        globalAnalysis.setValue( spatialModel + "_H-SDI", numCurrentNucleus, sdis[jj] );
      if ( ( sdis[jj] < 0 || sdis[jj] > 1 )  && ( jj = 3 ) )
        globalAnalysis.setValue( spatialModel + "_B-SDI", numCurrentNucleus, sqrt(-1) );
      else
        globalAnalysis.setValue( spatialModel + "_B-SDI", numCurrentNucleus, sdis[jj] );
      if ( ( sdis[jj] < 0 || sdis[jj] > 1 )  && ( jj = 4 ) )
        globalAnalysis.setValue( spatialModel + "_C-SDI", numCurrentNucleus, sqrt(-1) );
      else
        globalAnalysis.setValue( spatialModel + "_C-SDI", numCurrentNucleus, sdis[jj] );
      if ( ( sdis[jj] < 0 || sdis[jj] > 1 )  && ( jj = 5 ) )
        globalAnalysis.setValue( spatialModel + "_Z-SDI", numCurrentNucleus, sqrt(-1) );
      else
        globalAnalysis.setValue( spatialModel + "_Z-SDI", numCurrentNucleus, sdis[jj] );
      if ( ( sdis[jj] < 0 || sdis[jj] > 1 )  && ( jj = 6 ) )
        globalAnalysis.setValue( spatialModel + "_LRD-SDI", numCurrentNucleus, sqrt(-1) );
      else
        globalAnalysis.setValue( spatialModel + "_LRD-SDI", numCurrentNucleus, sdis[jj] );
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
    dataSet.setValue( "Z-SDI", row, sdis[5] );
//    dataSet.setValue( "Z-maxDiff", row, maxDiff[5] );
    dataSet.setValue( "LRD-SDI", row, sdis[6] );
//    dataSet.setValue( "LRD-maxDiff", row, maxDiff[5] );

    ostringstream iss; //we have 4 constraints
    iss << constraints;
    saveTest.save( analysisDir + iss.str() + "/" + filename + ".data", true );

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

  globalAnalysis.save( analysisDir + "/" + "nuclei_complete2.data", true );


}


/*! Preparing constraints for the model
****************************************************************/
void evaluator_completeSpatialRandomness(
  const string& filename, const string& parentDir, const string& spatialModel,
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
    filename, parentDir, spatialModel,
    function, constraints,
    dataSet, randomGenerator );
}

void evaluator_sizeConstrained(
  const string& filename, const string& parentDir, const string& spatialModel,
  const string& function, const int& constraints,
  DataSet& dataSet,
  RandomGenerator& randomGenerator)
{
  PRINT("spatialModelEvaluator_sizeConstrained");

  //open data info
  const string analysisDir = parentDir + "/analysis/";
  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + ".tm" );

  const DataSet ccsInfo( analysisDir + "ccs.data" );

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
    EVAL("Nucleus not found");
    return;
  }

  Vector<float> eqRadiiTemp( numCCS );
  Vector<float> distancesToBorder( numCCS );
  int k = 0;
  for ( int j = lastPos - numCCS + 1 ; j < lastPos + 1; ++j, ++k )
  {
    eqRadiiTemp[k] = ccsInfo.getValue<float>( "equivalentRadius_ZprojCorrection", j );
//    eqRadiiTemp[k] = ccsInfo.getValue<float>( "equivalentRadius_PSFVolCorrection", j );
    distancesToBorder[k] = ccsInfo.getValue<float>( "distanceToTheBorder", j );
  }


  for ( int i = 0; i < eqRadiiTemp.getSize(); ++i )
    if ( eqRadiiTemp[i] > distancesToBorder[i] )
      eqRadiiTemp[i] = distancesToBorder[i];

  const Vector<float> eqRadii = eqRadiiTemp;
  EVAL(eqRadii);

  SpatialModelHardcoreDistance3D<float> triMeshSpatialModel;
  triMeshSpatialModel.setRandomGenerator( randomGenerator );
  triMeshSpatialModel.setTriMesh( nucleusTriMesh );
  triMeshSpatialModel.setHardcoreDistances( eqRadii );
  triMeshSpatialModel.initialize();

  evaluator(
    nucleusTriMesh,
    triMeshSpatialModel,
    filename, parentDir, spatialModel,
    function, constraints,
    dataSet, randomGenerator );
}

void evaluator_distanceConstrained(
  const string& filename, const string& parentDir, const string& spatialModel,
  const string& function, const int& constraints,
  DataSet& dataSet,
  RandomGenerator& randomGenerator)
{
  PRINT("spatialModelEvaluator_distanceConstrained");

  //open data info
  const string analysisDir = parentDir + "/analysis/";
  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + ".tm" );
  const DataSet datasetNucleus( analysisDir + filename + "_chromocenters.data" );
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
    filename, parentDir, spatialModel,
    function, constraints,
    dataSet, randomGenerator );
}

void evaluator_sizeAndDistanceConstrained(
  const string& filename, const string& parentDir, const string& spatialModel,
  const string& function, const int& constraints,
  DataSet& dataSet,
  RandomGenerator& randomGenerator)
{
  PRINT("spatialModelEvaluator_sizeAndDistanceConstrained");

  //open data info
  const string analysisDir = parentDir + "/analysis/";
  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + ".tm" );

  const DataSet ccsInfo( analysisDir + "ccs.data" );

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
    EVAL("Nucleus not found");
    return;
  }

  Vector<float> eqRadiiTemp( numCCS );
  Vector<float> distancesToBorder( numCCS );
  int k = 0;
  for ( int j = lastPos - numCCS + 1 ; j < lastPos + 1; ++j, ++k )
  {
    eqRadiiTemp[k] = ccsInfo.getValue<float>( "equivalentRadius_ZprojCorrection", j );
//    eqRadiiTemp[k] = ccsInfo.getValue<float>( "equivalentRadius_PSFVolCorrection", j );
    distancesToBorder[k] = ccsInfo.getValue<float>( "distanceToTheBorder", j );
  }


  for ( int i = 0; i < eqRadiiTemp.getSize(); ++i )
    if ( eqRadiiTemp[i] > distancesToBorder[i] )
      eqRadiiTemp[i] = distancesToBorder[i];

  const Vector<float> eqRadii = eqRadiiTemp;
  EVAL(eqRadii);
  EVAL(distancesToBorder);

  SpatialModelHardcoreBorderDistance3D<float> triMeshSpatialModel;
  triMeshSpatialModel.setRandomGenerator( randomGenerator );
  triMeshSpatialModel.setTriMesh( nucleusTriMesh );
  triMeshSpatialModel.setDistancesToBorder( distancesToBorder );
  triMeshSpatialModel.setHardcoreDistances( eqRadii );
  triMeshSpatialModel.initialize();

  evaluator(
    nucleusTriMesh, triMeshSpatialModel, filename, parentDir, spatialModel,
    function, constraints, dataSet, randomGenerator );
}

void evaluator_MaximalRepulsionConstrained(
  const string& filename, const string& parentDir, const string& spatialModel,
  const string& function, const int& constraints,
  DataSet& dataSet,
  RandomGenerator& randomGenerator)
{
  PRINT("spatialModelEvaluator_MaximalRepulsionConstrained");

  //open data info
  const string analysisDir = parentDir + "/analysis/";
  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + ".tm" );

  const DataSet ccsInfo( analysisDir + "ccs.data" );

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
    EVAL("Nucleus not found");
    return;
  }

  Vector<float> eqRadiiTemp( numCCS );
  Vector<float> distancesToBorder( numCCS );
  int k = 0;
  for ( int j = lastPos - numCCS + 1 ; j < lastPos + 1; ++j, ++k )
  {
    eqRadiiTemp[k] = ccsInfo.getValue<float>( "equivalentRadius_ZprojCorrection", j );
//    eqRadiiTemp[k] = ccsInfo.getValue<float>( "equivalentRadius_PSFVolCorrection", j );
    distancesToBorder[k] = ccsInfo.getValue<float>( "distanceToTheBorder", j );
  }


  for ( int i = 0; i < eqRadiiTemp.getSize(); ++i )
    if ( eqRadiiTemp[i] > distancesToBorder[i] )
      eqRadiiTemp[i] = distancesToBorder[i];

  const Vector<float> eqRadii = eqRadiiTemp;
  EVAL(eqRadii);
  EVAL(distancesToBorder);


  SpatialModelMaximalRepulsion3D<float> triMeshSpatialModel;
  triMeshSpatialModel.setRandomGenerator( randomGenerator );
  triMeshSpatialModel.setTriMesh( nucleusTriMesh );
  triMeshSpatialModel.setNumMonteCarloCycles( 2000 );
  triMeshSpatialModel.setHardcoreDistances( eqRadii );
  triMeshSpatialModel.initialize();
  triMeshSpatialModel.initializeBeta( numCCS );

  evaluator(
    nucleusTriMesh, triMeshSpatialModel, filename, parentDir, spatialModel,
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
//  const DataSet datasetNucleus( analysisDir + filename + "_chromocenters.data" );
//  //old data
//  //const DataSet datasetNucleus( analysisDir + filename + ".data" );
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
  string spatialModel;
  switch( constraints )
  {
    case 0:
      spatialModel = "SpatialModelCompleteRandomness3D";
      evaluator_completeSpatialRandomness(
        filename, parentDir, spatialModel,
        function, constraints,
        dataSet, randomGenerator );
      break;
    case 1:
      spatialModel = "SpatialModelHardcoreDistance3D";
      evaluator_sizeConstrained(
        filename, parentDir, spatialModel,
        function, constraints,
        dataSet, randomGenerator );
      break;
    case 2:
      spatialModel = "SpatialModelBorderDistance3D";
      evaluator_distanceConstrained(
        filename, parentDir, spatialModel,
        function, constraints,
        dataSet, randomGenerator );
      break;
    case 3:
      spatialModel = "SpatialModelHardcoreBorderDistance3D";
      evaluator_sizeAndDistanceConstrained(
        filename, parentDir, spatialModel,
        function, constraints,
        dataSet, randomGenerator );
      break;
    case 4:
      spatialModel = "SpatialModelMaximalRepulsion3D";
      evaluator_MaximalRepulsionConstrained(
        filename, parentDir, spatialModel,
        function, constraints,
        dataSet, randomGenerator );
      break;
//    case 5:
//    evaluator_MaxRepulsionWithDistanceToTheBorderConstrained(
//      filename, parentDir, spatialModel,
//      function, constraints,
//      dataSet, randomGenerator );
//    break;
  }
}
