#include <spatialmodelcompleterandomness3d.h>
#include <spatialmodelborderdistance3d.h>
#include <spatialmodelhardcoreborderdistance3d.h>
#include <spatialmodelhardcoredistance3d.h>
#include <spatialmodelmaximalrepulsion3d.h>
#include <spatialmodel.h>
#include <trimesh.h>
#include <voxelmatrix.h>
#include <fileinfo.h>
#include <dataset.h>

#define TRACE
#include <trace.h>
using namespace std;

/*! Preparing constraints for the model
****************************************************************/
void evaluator_completeSpatialRandomness(
    const string& filename, const string& parentDir,
    const int& numMCSimulations, RandomGenerator& randomGenerator)
{
  PRINT("generatePatterns_completeSpatialRandomness");

  //open data info
  const string analysisDir = parentDir + "/analysis/";
  const string patternsDir = parentDir + "/patterns/SpatialModelCompleteRandomness3D/";
  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + ".tm" );
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
    EVAL("Nucleus not found");
    return;
  }

  SpatialModelCompleteRandomness3D<float> triMeshSpatialModel;
  triMeshSpatialModel.setRandomGenerator( randomGenerator );
  triMeshSpatialModel.setTriMesh( nucleusTriMesh );
  triMeshSpatialModel.initialize();

  VertexStack<float> vertexStack;//( 3, numCCS, 2*numMCSimulations, 0, 0 );
  Vertices<float> vertices( 3, numCCS, 0, 0 );

  EVAL(vertexStack.getSize());

  for ( int jj = 0; jj < 2*numMCSimulations; ++jj )
  {
    vertices = triMeshSpatialModel.drawSample(numCCS);
    vertexStack.insert( jj, vertices );
  }

  vertexStack.save( patternsDir + filename + ".vs", true );
}

void evaluator_sizeConstrained(
    const string& filename, const string& parentDir,
    const int& numMCSimulations, RandomGenerator& randomGenerator)
{
  PRINT("generatePatterns_sizeConstrained");

  //open data info
  const string analysisDir = parentDir + "/analysis/";
  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + ".tm" );

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

  VertexStack<float> vertexStack;//( 3, numCCS, 2*numMCSimulations, 0, 0 );
  Vertices<float> vertices( 3, numCCS, 0, 0 );

  EVAL(vertexStack.getSize());

  for ( int jj = 0; jj < 2*numMCSimulations; ++jj )
  {
    vertices = triMeshSpatialModel.drawSample(numCCS);
    vertexStack.insert( jj, vertices );
  }

  const string patternsDir = parentDir + "/patterns/SpatialModelHardcoreDistance3D/";
  vertexStack.save( patternsDir + filename + ".vs", true );
}

void evaluator_distanceConstrained(
    const string& filename, const string& parentDir,
    const int& numMCSimulations, RandomGenerator& randomGenerator)
{
  PRINT("spatialModelEvaluator_distanceConstrained");

  //open data info
  const string analysisDir = parentDir + "/analysis/";
  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + ".tm" );
  const DataSet datasetNucleus( analysisDir + filename + "_chromocenters.csv" );
  const Vector<float> distancesToBorder = datasetNucleus.getValues<float>( "distanceToTheBorder" );
  EVAL(distancesToBorder);

  const int numCCS = distancesToBorder.getSize();

  SpatialModelBorderDistance3D<float> triMeshSpatialModel;
  triMeshSpatialModel.setRandomGenerator( randomGenerator );
  triMeshSpatialModel.setTriMesh( nucleusTriMesh );
  triMeshSpatialModel.setDistancesToBorder( distancesToBorder );
  triMeshSpatialModel.initialize();

  VertexStack<float> vertexStack;//( 3, numCCS, 2*numMCSimulations, 0, 0 );
  Vertices<float> vertices( 3, numCCS, 0, 0 );

  EVAL(vertexStack.getSize());

  for ( int jj = 0; jj < 2*numMCSimulations; ++jj )
  {
    vertices = triMeshSpatialModel.drawSample(numCCS);
    vertexStack.insert( jj, vertices );
  }

  const string patternsDir = parentDir + "/patterns/SpatialModelBorderDistance3D/";
  vertexStack.save( patternsDir + filename + ".vs", true );
}

void evaluator_sizeAndDistanceConstrained(
    const string& filename, const string& parentDir,
    const int& numMCSimulations, RandomGenerator& randomGenerator)
{
  PRINT("spatialModelEvaluator_sizeAndDistanceConstrained");

  //open data info
  const string analysisDir = parentDir + "/analysis/";
  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + ".tm" );

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

  VertexStack<float> vertexStack;//( 3, numCCS, 2*numMCSimulations, 0, 0 );
  Vertices<float> vertices( 3, numCCS, 0, 0 );

  EVAL(vertexStack.getSize());

  for ( int jj = 0; jj < 2*numMCSimulations; ++jj )
  {
    vertices = triMeshSpatialModel.drawSample(numCCS);
    vertexStack.insert( jj, vertices );
  }

  const string patternsDir = parentDir + "/patterns/SpatialModelHardcoreBorderDistance3D/";
  vertexStack.save( patternsDir + filename + ".vs", true );

}

void evaluator_MaximalRepulsionConstrained(
    const string& filename, const string& parentDir,
    const int& numMCSimulations, RandomGenerator& randomGenerator)
{
  PRINT("spatialModelEvaluator_MaximalRepulsionConstrained");

  //open data info
  const string analysisDir = parentDir + "/analysis/";
  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + ".tm" );

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

  VertexStack<float> vertexStack;//( 3, numCCS, 2*numMCSimulations, 0, 0 );
  Vertices<float> vertices( 3, numCCS, 0, 0 );

  EVAL(vertexStack.getSize());

  for ( int jj = 0; jj < 2*numMCSimulations; ++jj )
  {
    vertices = triMeshSpatialModel.drawSample(numCCS);
    vertexStack.insert( jj, vertices );
  }

  const string patternsDir = parentDir + "/patterns/SpatialModelMaximalRepulsion3D/";
  vertexStack.save( patternsDir + filename + ".vs", true );

}



/*! Chooses case depending on the constraints
****************************************************************/
void generatePatterns(
  const string& filename, const string& parentDir,
  const int& constraints, const int& monteCarloSimulations, RandomGenerator& randomGenerator)
{
  switch( constraints )
  {
    case 0:
      evaluator_completeSpatialRandomness(
        filename, parentDir,
        monteCarloSimulations, randomGenerator );
      break;
    case 1:
      evaluator_sizeConstrained(
        filename, parentDir,
        monteCarloSimulations, randomGenerator );
      break;
    case 2:
      evaluator_distanceConstrained(
        filename, parentDir,
        monteCarloSimulations, randomGenerator );
      break;
    case 3:
      evaluator_sizeAndDistanceConstrained(
        filename, parentDir,
        monteCarloSimulations, randomGenerator );
      break;
    case 4:
      evaluator_MaximalRepulsionConstrained(
        filename, parentDir,
        monteCarloSimulations, randomGenerator );
      break;
  }
}
