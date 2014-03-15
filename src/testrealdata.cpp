#include <spatialmodelevaluator.h>
#include <spatialdescriptorfunctionf.h>
#include <spatialdescriptorfunctiong.h>
#include <spatialdescriptorfunctionh.h>
#include "spatialdescriptorborder.h"
#include "trimeshspatialmodel.h"
#include <trimesh.h>
#include <voxelmatrix.h>
#include <fileinfo.h>
//#include <regionanalysis.h>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <string.h>

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
  else //if ( function == "F" )
  {
    PRINT("F");
    TriMeshSpatialModel<float> tempTriMeshSpatialModel;
    tempTriMeshSpatialModel.setRandomGenerator( randomGenerator );
    tempTriMeshSpatialModel.setTriMesh( nucleusTriMesh );
    tempTriMeshSpatialModel.initialize();
    Vertices<float> evaluationPositions = tempTriMeshSpatialModel.drawSample( 10000 );

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
  }

  ostringstream iss; //we suppose as much 99 labels
  iss << constraints;

  float pValue = modelEvaluator.eval( vertices, &saveTest );
  saveTest.save( analysisDir + iss.str() + "/" + function + "/" + iss.str() + "/" +filename + ".csv", true );
  dataSet.setValue( "pValues", dataSet.size()[0]+1, pValue );

  EVAL( pValue );

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

  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + "_nucleus.tm" );
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

  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + "_nucleus.tm" );
  const string analysisDir = parentDir + "/analysis/";
  DataSet datasetNucleus( analysisDir + filename + "_chromocenters.csv" );
  Vector<float> eqRadii = datasetNucleus.getValues<float>( "equivalentRadius_tm" );

  TriMeshSpatialModel<float> triMeshSpatialModel;
  triMeshSpatialModel.setRandomGenerator( randomGenerator );
  triMeshSpatialModel.setTriMesh( nucleusTriMesh );
  triMeshSpatialModel.setHardcoreDistance( eqRadii );
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

  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + "_nucleus.tm" );
  const string analysisDir = parentDir + "/analysis/";
  const DataSet datasetNucleus( analysisDir + filename + "_chromocenters.csv" );
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

  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + "_nucleus.tm" );
  const string analysisDir = parentDir + "/analysis/";
  const DataSet datasetNucleus( analysisDir + filename + "_chromocenters.csv" );
  const Vector<float> eqRadii = datasetNucleus.getValues<float>( "equivalentRadius_tm" );
  const Vector<float> distancesToBorder = datasetNucleus.getValues<float>( "distanceToTheBorder" );

  TriMeshSpatialModel<float> triMeshSpatialModel;
  triMeshSpatialModel.setRandomGenerator( randomGenerator );
  triMeshSpatialModel.setTriMesh( nucleusTriMesh );
  triMeshSpatialModel.setDistanceToBorder( distancesToBorder );
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
  }
}
