//#include <spatialmodelevaluator.h>
#include "spatialmodelevaluator2.h"
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


void spatialModelEvaluator(
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
  const int numSamples = 100;

  SpatialModelEvaluator<float,float> modelEvaluator;
  modelEvaluator.setModel( triMeshSpatialModel );
  modelEvaluator.setNumMonteCarloSamples( numPatterns );//to check uniformity
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

  //for (int i = 0; i < numSamples; ++i)
  for (int i = 0; i < 1; ++i)
  {
    EVAL( i );
    ostringstream oss; //we suppose as much 99 labels
    oss << setw(2) << setfill('0') << i;
    ostringstream iss; //we suppose as much 99 labels
    iss << constraints;
    Vertices<float> vertices = triMeshSpatialModel.drawSample( numPoints );
    vertices.save( parentDir + "/" + filename + ".vx", true );
    //float pValue = modelEvaluator.eval( vertices, &saveTest );
    //saveTest.save( analysisDir + iss.str() + "/" + function + "/" + filename + "_" + oss.str() + ".csv", true );
    //dataSet.setValue( "pValues", i, pValue );
    //EVAL( pValue );
  }
}


/*! Preparing constraints for the model
****************************************************************/
void spatialModelEvaluator_completeSpatialRandomness(
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

  spatialModelEvaluator(
    nucleusTriMesh,
    triMeshSpatialModel,
    filename, parentDir,
    function, constraints,
    dataSet, randomGenerator );
}

void spatialModelEvaluator_sizeConstrained(
  const string& filename, const string& parentDir,
  const string& function, const int& constraints,
  DataSet& dataSet,
  RandomGenerator& randomGenerator)
{
  PRINT("spatialModelEvaluator_sizeConstrained");

  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + ".tm" );
  const string analysisDir = parentDir + "/analysis/";
  DataSet datasetNucleus( analysisDir + filename + "_chromocenters.csv" );
  Vector<float> eqRadii = datasetNucleus.getValues<float>( "equivalentRadius_tm" );

  TriMeshSpatialModel<float> triMeshSpatialModel;
  triMeshSpatialModel.setRandomGenerator( randomGenerator );
  triMeshSpatialModel.setTriMesh( nucleusTriMesh );
  triMeshSpatialModel.setHardcoreDistances( eqRadii );
  triMeshSpatialModel.initialize();

  spatialModelEvaluator(
    nucleusTriMesh,
    triMeshSpatialModel,
    filename, parentDir,
    function, constraints,
    dataSet, randomGenerator );
}

void spatialModelEvaluator_distanceConstrained(
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

  spatialModelEvaluator(
    nucleusTriMesh,
    triMeshSpatialModel,
    filename, parentDir,
    function, constraints,
    dataSet, randomGenerator );
}

void spatialModelEvaluator_sizeAndDistanceConstrained(
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
  triMeshSpatialModel.setHardcoreDistances( eqRadii );
  triMeshSpatialModel.initialize();

  spatialModelEvaluator(
    nucleusTriMesh, triMeshSpatialModel, filename, parentDir,
    function, constraints, dataSet, randomGenerator );
}


/*! Chooses case depending on the constraints
****************************************************************/
void spatialModelEvaluator(
  const string& filename, const string& parentDir,
  const string& function, const int constraints,
  DataSet& dataSet,
  RandomGenerator& randomGenerator)
{
  switch( constraints )
  {
    case 0:
      spatialModelEvaluator_completeSpatialRandomness(
        filename, parentDir,
        function, constraints,
        dataSet, randomGenerator );
      break;
    case 1:
      spatialModelEvaluator_sizeConstrained(
        filename, parentDir,
        function, constraints,
        dataSet, randomGenerator );
      break;
    case 2:
      spatialModelEvaluator_distanceConstrained(
        filename, parentDir,
        function, constraints,
        dataSet, randomGenerator );
      break;
    case 3:
      spatialModelEvaluator_sizeAndDistanceConstrained(
        filename, parentDir,
        function, constraints,
        dataSet, randomGenerator );
      break;
  }
}
