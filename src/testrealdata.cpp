#include <spatialmodelevaluator.h>
#include <spatialdescriptorfunctionf.h>
#include <spatialdescriptorfunctiong.h>
#include <spatialdescriptorfunctionh.h>
#include "spatialdescriptorborder.h"
#include "maximalrepulsion.h"
#include "spatialdescriptormaxima.h"
#include "trimeshspatialmodel.h"
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
  TriMeshSpatialModel<float>& triMeshSpatialModel,
  const string& filename, const string& parentDir,
  const string& function, const int& constraints,
  DataSet& dataSet, RandomGenerator& randomGenerator)
{
  const string analysisDir = parentDir + "/analysis/";
  string classif = parentDir;
  classif = classif.substr(classif.find_last_of("/\\")+1,classif.length());

  //new data
  const DataSet datasetNucleus( analysisDir + filename + "_chromocenters.csv" );

  const int numPoints = datasetNucleus.size()[0];

  const int numPatterns = 99;

  SpatialModelEvaluator<float,float> modelEvaluator;
  modelEvaluator.setModel( triMeshSpatialModel );
  modelEvaluator.setNumRandomSamples( numPatterns ); //to check uniformity
  modelEvaluator.setPrecision( 0.05 );

  SpatialDescriptor<float>* spatialDescriptor;

  DataSet saveTest;

  if ( function == "G" )
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
    SpatialDescriptorDistanceToBorder<float>* spatialDescriptorDistanceToBorder;
    spatialDescriptorDistanceToBorder = new SpatialDescriptorDistanceToBorder<float>();
    spatialDescriptorDistanceToBorder->setTriMesh( nucleusTriMesh );
    spatialDescriptor = spatialDescriptorDistanceToBorder;
  }
  else if ( function == "FMod" )
  {
    PRINT("FMod");
    TriMeshSpatialModel<float> tempTriMeshSpatialModel;
    tempTriMeshSpatialModel.setRandomGenerator( randomGenerator );
    tempTriMeshSpatialModel.setTriMesh( nucleusTriMesh );
    const DataSet datasetNucleus( analysisDir + filename + ".csv" );
    const Vector<float> eqRadii = datasetNucleus.getValues<float>( "chromocenterRadius" );
    tempTriMeshSpatialModel.setDistanceToBorderRange( eqRadii.min() , eqRadii.max() +10  );
    tempTriMeshSpatialModel.initialize();
    Vertices<float> evaluationPositions = tempTriMeshSpatialModel.drawSample( 10000 );
    evaluationPositions.save( parentDir + "/spatial_models/" + filename + "_FModpattern.vx", true);
    SpatialDescriptorFunctionF<float>* spatialDescriptorFunctionF;
    spatialDescriptorFunctionF = new SpatialDescriptorFunctionF<float>();
    spatialDescriptorFunctionF->setEvaluationPositions( evaluationPositions );
    spatialDescriptor = spatialDescriptorFunctionF;
  }
  else //if ( function == "F" )
  {
    PRINT("F");
    TriMeshSpatialModel<float> tempTriMeshSpatialModel;
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

  float pValue = modelEvaluator.eval( vertices, &saveTest );
  EVAL( pValue );

  saveTest.save( analysisDir + iss.str() + "/" + function + "/" + filename + ".csv", true );
  const int row = dataSet.size()[0];
  dataSet.setValue( "nucleus", row, filename );
  dataSet.setValue( "class", row, classif );//classification: mutant, tissue, etc.
  dataSet.setValue( "pValues", row, pValue );


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
  TriMeshSpatialModel<float> triMeshSpatialModel;
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

  TriMeshSpatialModel<float> triMeshSpatialModel;
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
  //const DataSet datasetNucleus( analysisDir + filename + "_chromocenters.csv" );
  const DataSet datasetNucleus( analysisDir + filename + ".csv" );
  const Vector<float> distancesToBorder = datasetNucleus.getValues<float>( "distanceToTheBorder" );

  TriMeshSpatialModel<float> triMeshSpatialModel;
  triMeshSpatialModel.setRandomGenerator( randomGenerator );
  triMeshSpatialModel.setTriMesh( nucleusTriMesh );
  triMeshSpatialModel.setDistanceToBorder( distancesToBorder );
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

  TriMeshSpatialModel<float> triMeshSpatialModel;
  triMeshSpatialModel.setRandomGenerator( randomGenerator );
  triMeshSpatialModel.setTriMesh( nucleusTriMesh );
  triMeshSpatialModel.setDistanceToBorder( distancesToBorder );
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

  MaximalRepulsionTriMeshSpatialModel<float> triMeshSpatialModel;
  triMeshSpatialModel.setRandomGenerator( randomGenerator );
  triMeshSpatialModel.setTriMesh( nucleusTriMesh );
  triMeshSpatialModel.setHardcoreDistance( eqRadii );
  triMeshSpatialModel.initialize();

  evaluator(
    nucleusTriMesh, triMeshSpatialModel, filename, parentDir,
    function, constraints, dataSet, randomGenerator );
}

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
  }
}
