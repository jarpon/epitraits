//#include <spatialmodelevaluator.h>
#include "spatialmodelevaluator2.h"
#include <spatialdescriptorfunctionf.h>
#include "spatialdescriptorfunctiongg.h"
#include <spatialdescriptorfunctionh.h>
#include "spatialdescriptorborder.h"
#include "spatialdescriptorcentroid.h"
#include "maximalrepulsion.h"
#include "maxrepulsionwithdistances.h"
#include "spatialdescriptormaxima.h"
#include "trimeshspatialmodel2.h"
#include <trimesh.h>
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


void evaluator(
  const TriMesh<float>& nucleusTriMesh,
  TriMeshSpatialModel2<float>& triMeshSpatialModel,
  const string& filename, const string& parentDir,
  const string& function, const int& constraints,
  DataSet& dataSet, RandomGenerator& randomGenerator)
{
  const string analysisDir = parentDir + "/analysis/";
  string classif = parentDir;
  classif = classif.substr(classif.find_last_of("/\\")+1,classif.length());

  //new data
  const DataSet datasetNucleus( analysisDir + filename + "_chromocenters.csv" );
  const DataSet datasetKinda2( analysisDir + filename + "_kinda2.csv" );

  const Vector<float> eqRadii = datasetNucleus.getValues<float>( "equivalentRadius_tm" );
  const Vector<float> eqRadiiKinda2 = datasetKinda2.getValues<float>( "equivalentRadius_tm" );

  Vertices<float> vertices1 ( 3, datasetNucleus.size()[0], 0, 0 );
  for ( int i = 0; i < datasetNucleus.size()[0]; ++i )
  {
    vertices1[i][0] = datasetNucleus.getValue<float>( "centroidCoordX", i );
    vertices1[i][1] = datasetNucleus.getValue<float>( "centroidCoordY", i );
    vertices1[i][2] = datasetNucleus.getValue<float>( "centroidCoordZ", i );
    EVAL(vertices1[i]);
  }

  Vertices<float> vertices2 ( 3, datasetKinda2.size()[0], 0, 0 );
  for ( int i = 0; i < datasetKinda2.size()[0]; ++i )
  {
    vertices2[i][0] = datasetKinda2.getValue<float>( "centroidCoordX", i );
    vertices2[i][1] = datasetKinda2.getValue<float>( "centroidCoordY", i );
    vertices2[i][2] = datasetKinda2.getValue<float>( "centroidCoordZ", i );
    EVAL(vertices2[i]);
  }

  const int numPoints = datasetNucleus.size()[0];

  //const int numPatterns = 99;
  const int numPatterns = 1;

  SpatialModelEvaluator<float,float> modelEvaluator;
  modelEvaluator.setModel( triMeshSpatialModel );
  modelEvaluator.setNumMonteCarloSamples( numPatterns ); //to check uniformity
  modelEvaluator.setPrecision( 0.05 );

  SpatialDescriptor<float>* spatialDescriptor;

  DataSet saveTest;

  if ( function == "all" )
  {
    PRINT("all functions");
    TriMeshSpatialModel2<float> tempTriMeshSpatialModel2;
    tempTriMeshSpatialModel2.setRandomGenerator( randomGenerator );
    tempTriMeshSpatialModel2.setTriMesh( nucleusTriMesh );
    tempTriMeshSpatialModel2.initialize();
    Vertices<float> evaluationPositions = tempTriMeshSpatialModel2.drawSample( 10000 );

    SpatialDescriptorFunctionF<float>* spatialDescriptorFunctionF;
    spatialDescriptorFunctionF = new SpatialDescriptorFunctionF<float>();
    spatialDescriptorFunctionF->setEvaluationPositions( evaluationPositions );
    spatialDescriptor = spatialDescriptorFunctionF;
    modelEvaluator.addDescriptor( *spatialDescriptor );

    spatialDescriptor = new SpatialDescriptorFunctionGG<float>();
    modelEvaluator.addDescriptor( *spatialDescriptor );

    spatialDescriptor = new SpatialDescriptorFunctionH<float>();
    modelEvaluator.addDescriptor( *spatialDescriptor );

    SpatialDescriptorDistanceToBorder<float>* spatialDescriptorDistanceToBorder;
    spatialDescriptorDistanceToBorder = new SpatialDescriptorDistanceToBorder<float>();
    spatialDescriptorDistanceToBorder->setTriMesh( nucleusTriMesh );
    spatialDescriptor = spatialDescriptorDistanceToBorder;
    modelEvaluator.addDescriptor( *spatialDescriptor );

    SpatialDescriptorDistanceToCentroid<float>* spatialDescriptorDistanceToCentroid;
    spatialDescriptorDistanceToCentroid = new SpatialDescriptorDistanceToCentroid<float>();
    spatialDescriptorDistanceToCentroid->setTriMesh( nucleusTriMesh );
    spatialDescriptor = spatialDescriptorDistanceToCentroid;
    modelEvaluator.addDescriptor( *spatialDescriptor );

  }
  else if ( function == "G" )
  {
    PRINT("G'");
    SpatialDescriptorFunctionGG<float>* spatialDescriptorGG;
    spatialDescriptorGG = new SpatialDescriptorFunctionGG<float>();
//    spatialDescriptorGG->setVerticesKind1( vertices1 );
//    spatialDescriptorGG->setVerticesKind2( vertices2 );
    spatialDescriptorGG->setVertices( vertices1, vertices2 );
    spatialDescriptor = spatialDescriptorGG;
  }
  else if ( function == "H" )
  {
    PRINT("H");
    spatialDescriptor = new SpatialDescriptorFunctionH<float>();
  }
  else if ( function == "B" )
  {
    PRINT("B");
    SpatialDescriptorDistanceToBorder<float>* spatialDescriptorDistanceToBorder;
    spatialDescriptorDistanceToBorder = new SpatialDescriptorDistanceToBorder<float>();
    spatialDescriptorDistanceToBorder->setTriMesh( nucleusTriMesh );
    spatialDescriptor = spatialDescriptorDistanceToBorder;
  }
  else if ( function == "C" )
  {
    PRINT("C");
    SpatialDescriptorDistanceToCentroid<float>* spatialDescriptorDistanceToCentroid;
    spatialDescriptorDistanceToCentroid = new SpatialDescriptorDistanceToCentroid<float>();
    spatialDescriptorDistanceToCentroid->setTriMesh( nucleusTriMesh );
    spatialDescriptor = spatialDescriptorDistanceToCentroid;
  }
  else if ( function == "FMod" )
  {
    PRINT("FMod");
    TriMeshSpatialModel2<float> tempTriMeshSpatialModel2;
    tempTriMeshSpatialModel2.setRandomGenerator( randomGenerator );
    tempTriMeshSpatialModel2.setTriMesh( nucleusTriMesh );
    const DataSet datasetNucleus( analysisDir + filename + ".csv" );
    const Vector<float> eqRadii = datasetNucleus.getValues<float>( "chromocenterRadius" );
    tempTriMeshSpatialModel2.setDistanceToBorderRange( eqRadii.min() , eqRadii.max() +10  );
    tempTriMeshSpatialModel2.initialize();
    Vertices<float> evaluationPositions = tempTriMeshSpatialModel2.drawSample( 10000 );
    evaluationPositions.save( parentDir + "/spatial_models/" + filename + "_FModpattern.vx", true);
    SpatialDescriptorFunctionF<float>* spatialDescriptorFunctionF;
    spatialDescriptorFunctionF = new SpatialDescriptorFunctionF<float>();
    spatialDescriptorFunctionF->setEvaluationPositions( evaluationPositions );
    spatialDescriptor = spatialDescriptorFunctionF;
  }
  else //if ( function == "F" )
  {
    PRINT("F");
    TriMeshSpatialModel2<float> tempTriMeshSpatialModel2;
    tempTriMeshSpatialModel2.setRandomGenerator( randomGenerator );
    tempTriMeshSpatialModel2.setTriMesh( nucleusTriMesh );
    tempTriMeshSpatialModel2.initialize();
    Vertices<float> evaluationPositions = tempTriMeshSpatialModel2.drawSample( 10000 );
    //evaluationPositions.save( parentDir + "/spatial_models/" + filename + "_Fpattern-" + ".vx", true);
    SpatialDescriptorFunctionF<float>* spatialDescriptorFunctionF;
    spatialDescriptorFunctionF = new SpatialDescriptorFunctionF<float>();
    spatialDescriptorFunctionF->setEvaluationPositions( evaluationPositions );
    spatialDescriptor = spatialDescriptorFunctionF;
  }


  if ( function != "all" )
  {
    modelEvaluator.setDescriptor( *spatialDescriptor );

    ostringstream iss; //we suppose as much 99 labels
    iss << constraints;

  //  float pValue = modelEvaluator.eval( vertices, &saveTest );
  //  EVAL( pValue );

    Vertices<float> verticesAll ( 3, datasetNucleus.size()[0]+datasetKinda2.size()[0], 0, 0 );
    verticesAll.append( vertices1 );
    verticesAll.append( vertices2 );

    EVAL("here");
//    Vector<float> output = modelEvaluator.evalSDIandMaxDiff( verticesAll, &saveTest );
//    EVAL( output[0] );
//    EVAL( output[1] );

    vector<float> pValues;
    vector<int> ranks;
    vector<float> maxDiff;

    modelEvaluator.evalSDIandMaxDiff( vertices1, pValues, ranks, maxDiff);

    const int row = dataSet.numRows();

    saveTest.setValue( "nucleus", row, filename );
    saveTest.setValue( "class", row, classif );//classification: mutant, tissue, etc.
    saveTest.setValue( "G-SDI", row, pValues[0] );
    saveTest.setValue( "G-maxDiff", row, maxDiff[0] );

    saveTest.save( analysisDir + iss.str() + "/" + function + "/" + filename + "_2kind.csv", true );
  //  saveTest.save( analysisDir + iss.str() + "/" + function + "/" + filename + "_random.csv", true );
    dataSet.setValue( "nucleus", row, filename );
    dataSet.setValue( "class", row, classif );//classification: mutant, tissue, etc.
    dataSet.setValue( "descriptor", row, function );//spatial descriptor
    //dataSet.setValue( "index", row, pValue );

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

    vector<float> pValues;
    vector<int> ranks;
    vector<float> maxDiff;

    modelEvaluator.evalSDIandMaxDiff( vertices1, pValues, ranks, maxDiff);
//    modelEvaluator.eval( vertices, pValues, ranks);

    dataSet.setValue( "nucleus", row, filename );
    dataSet.setValue( "class", row, classif );//classification: mutant, tissue, etc.
    dataSet.setValue( "F-SDI", row, pValues[0] );
    dataSet.setValue( "F-maxDiff", row, maxDiff[0] );
    dataSet.setValue( "G-SDI", row, pValues[1] );
    dataSet.setValue( "G-maxDiff", row, maxDiff[1] );
    dataSet.setValue( "H-SDI", row, pValues[2] );
    dataSet.setValue( "H-maxDiff", row, maxDiff[2] );
    dataSet.setValue( "B-SDI", row, pValues[3] );
    dataSet.setValue( "B-maxDiff", row, maxDiff[3] );
    dataSet.setValue( "C-SDI", row, pValues[4] );
    dataSet.setValue( "C-maxDiff", row, maxDiff[4] );

    saveTest.setValue( "nucleus", 1, filename );
    saveTest.setValue( "class", 1, classif );//classification: mutant, tissue, etc.
    saveTest.setValue( "F-SDI", 1, pValues[0] );
    saveTest.setValue( "F-maxDiff", 1, maxDiff[0] );
    saveTest.setValue( "G-SDI", 1, pValues[1] );
    saveTest.setValue( "G-maxDiff", 1, maxDiff[1] );
    saveTest.setValue( "H-SDI", 1, pValues[2] );
    saveTest.setValue( "H-maxDiff", 1, maxDiff[2] );
    saveTest.setValue( "B-SDI", 1, pValues[3] );
    saveTest.setValue( "B-maxDiff", 1, maxDiff[3] );
    saveTest.setValue( "C-SDI", 1, pValues[4] );
    saveTest.setValue( "C-maxDiff", 1, maxDiff[4] );

    ostringstream iss; //we suppose as much 99 labels
    iss << constraints;

    saveTest.save( analysisDir + iss.str() + "/" + filename + "_all.csv", true );
  }



}


/*! Preparing constraints for the model
****************************************************************/
void completeSpatialRandomness(
  const string& filename, const string& parentDir,
  const string& function, const int& constraints,
  DataSet& dataSet,
  RandomGenerator& randomGenerator)
{
  PRINT("spatialModelcompleteSpatialRandomness");

  //const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + "_nucleus.tm" );
  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + ".tm" );
  TriMeshSpatialModel2<float> triMeshSpatialModel;
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

void sizeConstrained(
  const string& filename, const string& parentDir,
  const string& function, const int& constraints,
  DataSet& dataSet,
  RandomGenerator& randomGenerator)
{
  PRINT("spatialModelsizeConstrained");

  const string analysisDir = parentDir + "/analysis/";

  //new data
  //const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + "_nucleus.tm" );
  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + ".tm" );
  //const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/nuclei/" + filename + ".tm" );
  //old data
  //const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + "-nucleus.tm" );

  //new data
  const DataSet datasetCCs( analysisDir + filename + "_chromocenters.csv" );
  const DataSet datasetKinda2( analysisDir + filename + "_kinda2.csv" );

  //new data
  const Vector<float> eqRadii = datasetCCs.getValues<float>( "equivalentRadius_tm" );
  const Vector<float> eqRadiiKinda2 = datasetKinda2.getValues<float>( "equivalentRadius_tm" );
  //old data
  //const Vector<float> eqRadii = datasetNucleus.getValues<float>( "chromocenterRadius" );

  TriMeshSpatialModel2<float> triMeshSpatialModel;
  triMeshSpatialModel.setRandomGenerator( randomGenerator );
  triMeshSpatialModel.setTriMesh( nucleusTriMesh );
  triMeshSpatialModel.addDistribution( datasetCCs.numRows(), eqRadii );
  EVAL(eqRadii);
  triMeshSpatialModel.addDistribution( datasetKinda2.numRows(), eqRadiiKinda2 );
  EVAL(eqRadiiKinda2);
  triMeshSpatialModel.initialize();

  evaluator(
    nucleusTriMesh,
    triMeshSpatialModel,
    filename, parentDir,
    function, constraints,
    dataSet, randomGenerator );
}

//void distanceConstrained(
//  const string& filename, const string& parentDir,
//  const string& function, const int& constraints,
//  DataSet& dataSet,
//  RandomGenerator& randomGenerator)
//{
//  PRINT("spatialModeldistanceConstrained");

//  //const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + "_nucleus.tm" );
//  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + ".tm" );
//  const string analysisDir = parentDir + "/analysis/";
//  //const DataSet datasetNucleus( analysisDir + filename + "_chromocenters.csv" );
//  const DataSet datasetNucleus( analysisDir + filename + ".csv" );
//  const Vector<float> distancesToBorder = datasetNucleus.getValues<float>( "distanceToTheBorder" );

//  TriMeshSpatialModel2<float> triMeshSpatialModel;
//  triMeshSpatialModel.setRandomGenerator( randomGenerator );
//  triMeshSpatialModel.setTriMesh( nucleusTriMesh );
//  triMeshSpatialModel.setDistanceToBorder( distancesToBorder );
//  triMeshSpatialModel.initialize();

//  evaluator(
//    nucleusTriMesh,
//    triMeshSpatialModel,
//    filename, parentDir,
//    function, constraints,
//    dataSet, randomGenerator );
//}

//void sizeAndDistanceConstrained(
//  const string& filename, const string& parentDir,
//  const string& function, const int& constraints,
//  DataSet& dataSet,
//  RandomGenerator& randomGenerator)
//{
//  PRINT("spatialModelsizeAndDistanceConstrained");

//  //randomGenerator.init(11);
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

//  //new data
//  const Vector<float> eqRadii = datasetNucleus.getValues<float>( "equivalentRadius_tm" );
//  //old data
//  //const Vector<float> eqRadii = datasetNucleus.getValues<float>( "chromocenterRadius" );
//  EVAL(eqRadii);
//  const Vector<float> distancesToBorder = datasetNucleus.getValues<float>( "distanceToTheBorder" );

//  TriMeshSpatialModel2<float> triMeshSpatialModel;
//  triMeshSpatialModel.setRandomGenerator( randomGenerator );
//  triMeshSpatialModel.setTriMesh( nucleusTriMesh );
//  triMeshSpatialModel.setDistanceToBorder( distancesToBorder );
//  triMeshSpatialModel.setHardcoreDistances( eqRadii );
//  triMeshSpatialModel.initialize();

//  evaluator(
//    nucleusTriMesh, triMeshSpatialModel, filename, parentDir,
//    function, constraints, dataSet, randomGenerator );
//}

//void MaximalRepulsionConstrained(
//  const string& filename, const string& parentDir,
//  const string& function, const int& constraints,
//  DataSet& dataSet,
//  RandomGenerator& randomGenerator)
//{
//  PRINT("spatialModelMaximalRepulsionConstrained");

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

//  //new data
//  const Vector<float> eqRadii = datasetNucleus.getValues<float>( "equivalentRadius_tm" );
//  //old data
//  //const Vector<float> eqRadii = datasetNucleus.getValues<float>( "chromocenterRadius" );
//  EVAL(eqRadii);

//  MaximalRepulsionTriMeshSpatialModel2<float> triMeshSpatialModel;
//  triMeshSpatialModel.setRandomGenerator( randomGenerator );
//  triMeshSpatialModel.setTriMesh( nucleusTriMesh );
//  triMeshSpatialModel.setHardcoreDistance( eqRadii );
//  triMeshSpatialModel.initialize();

//  evaluator(
//    nucleusTriMesh, triMeshSpatialModel, filename, parentDir,
//    function, constraints, dataSet, randomGenerator );
//}


//void MaxRepulsionWithDistanceToTheBorderConstrained(
//  const string& filename, const string& parentDir,
//  const string& function, const int& constraints,
//  DataSet& dataSet,
//  RandomGenerator& randomGenerator)
//{
//  PRINT("spatialModelMaxRepulsionWithDistanceToTheBorderConstrained");

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

//  MaxRepulsionWithDistancesTriMeshSpatialModel2<float> triMeshSpatialModel;
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
void twoCompartmentsEvaluator(
  const string& filename, const string& parentDir,
  const string& function, const int& constraints,
  DataSet& dataSet,
  RandomGenerator& randomGenerator)
{
  switch( constraints )
  {
    case 0:
      completeSpatialRandomness(
        filename, parentDir,
        function, constraints,
        dataSet, randomGenerator );
      break;
    case 1:
      sizeConstrained(
        filename, parentDir,
        function, constraints,
        dataSet, randomGenerator );
      break;
//    case 2:
//      distanceConstrained(
//        filename, parentDir,
//        function, constraints,
//        dataSet, randomGenerator );
//      break;
//    case 3:
//      sizeAndDistanceConstrained(
//        filename, parentDir,
//        function, constraints,
//        dataSet, randomGenerator );
//      break;
//    case 4:
//      MaximalRepulsionConstrained(
//        filename, parentDir,
//        function, constraints,
//        dataSet, randomGenerator );
//      break;
//    case 5:
//    MaxRepulsionWithDistanceToTheBorderConstrained(
//      filename, parentDir,
//      function, constraints,
//      dataSet, randomGenerator );
//    break;
  }
}

