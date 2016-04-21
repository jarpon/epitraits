#include <spatialmodelcompleterandomness3d.h>
#include <spatialmodelborderdistance3d.h>
#include <spatialmodelhardcoreborderdistance3d.h>
#include <spatialmodelhardcoredistance3d.h>
//#include <spatialmodelmaximalrepulsion3d.h>
#include <spatialmodelmaximalrepulsion3d.h>
#include <spatialmodel.h>
#include <trimesh.h>
#include <voxelmatrix.h>
#include <fileinfo.h>
#include <dataset.h>
#include <sstream>

#define TRACE
#include <trace.h>
using namespace std;

/*! Preparing constraints for the model
****************************************************************/
void evaluatorLessCCS_completeSpatialRandomness(
    const string& filename, const string& parentDir,
    const int& numMCSimulations, RandomGenerator& randomGenerator)
{
  PRINT("generatePatterns_completeSpatialRandomness");
  PRINT("using 1 chromocenters less");

  //open data info
  const string analysisDir = parentDir + "/analysis/";
  //const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + ".tm" );
  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/nuclei/" + filename + ".tm" );
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

  SpatialModelCompleteRandomness3D<float> triMeshSpatialModel;
  triMeshSpatialModel.setRandomGenerator( randomGenerator );
  triMeshSpatialModel.setTriMesh( nucleusTriMesh );
  triMeshSpatialModel.initialize();

  const string patternsDir = parentDir + "/patterns/SpatialModelCompleteRandomness3D/";

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

void evaluatorLessCCS_sizeConstrained(
    const string& filename, const string& parentDir,
    const int& numMCSimulations, RandomGenerator& randomGenerator)
{
  PRINT("generatePatterns_sizeConstrained");

  //open data info
  const string analysisDir = parentDir + "/analysis/";
  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/nuclei/" + filename + ".tm" );
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

  //choose how many CCS are not taken into account
//  Vector<int> ccsToAvoid;
//  ccsToAvoid.setSize( 2 );
////  int ccsToAvoid;
////  ccsToAvoid = randomGenerator.uniformLF( 0, numCCS-1 );
//  ccsToAvoid[0] = randomGenerator.uniformLF( 0, numCCS-1 );
//  do
//  {
//    ccsToAvoid[1] = randomGenerator.uniformLF( 0, numCCS-1 );
//  } while ( ccsToAvoid[0] == ccsToAvoid[1] );


  for ( int i = 0; i < numCCS; ++i )
  {
//    EVAL( ccsToAvoid );

    Vertices<float> originalVertices( 3, 0, 0, 0 );
    Vector<float> eqRadiiTemp( numCCS );
    Vector<float> vertex( 3 );
    int k = 0;
    for ( int j = lastPos - numCCS + 1 ; j < lastPos + 1; ++j, ++k )
    {
        eqRadiiTemp[k] = ccsInfo.getValue<float>( "equivalentRadius_ZprojCorrection", j );
        vertex[X] = ccsInfo.getValue<float>( "centroidCoordX", j );
        vertex[Y] = ccsInfo.getValue<float>( "centroidCoordY", j );
        vertex[Z] = ccsInfo.getValue<float>( "centroidCoordZ", j );
        originalVertices.append( vertex );
    }

    EVAL(eqRadiiTemp);
    Vector<float> eqRadii;
    //eqRadii.setSize( numCCS - 3 );
    eqRadii.setSize( numCCS - 2 );
//      eqRadii.setSize( numCCS - 1 );
    int n = 0;

      Vector<int> ccsToAvoid;
      //ccsToAvoid.setSize( 3 );
      ccsToAvoid.setSize( 2 );
//        ccsToAvoid.setSize( 1 );
    //  int ccsToAvoid;
    //  ccsToAvoid = randomGenerator.uniformLF( 0, numCCS-1 );
      ccsToAvoid[0] = randomGenerator.uniformLF( 0, numCCS-1 );
      do
      {
        ccsToAvoid[1] = randomGenerator.uniformLF( 0, numCCS-1 );
//        ccsToAvoid[2] = randomGenerator.uniformLF( 0, numCCS-1 );
      //} while ( ccsToAvoid[0] == ccsToAvoid[1] == ccsToAvoid[2] );
        } while ( ccsToAvoid[0] == ccsToAvoid[1] );

    for ( int l = 0; l < numCCS; ++l )
      //if ( ( l != ccsToAvoid[0] ) && ( l != ccsToAvoid[1] ) && ( l != ccsToAvoid[2] ) )
      if ( ( l != ccsToAvoid[0] ) && ( l != ccsToAvoid[1] ) )
//      if ( ( l != ccsToAvoid[0] ) )
  //    if ( ( l != ccsToAvoid ) )
   //   if ( l != i )
      {
        eqRadii[n] = eqRadiiTemp[l];
        ++n;
      }
    originalVertices.remove( ccsToAvoid.max() );
//    originalVertices.remove( ccsToAvoid.max() );
    //originalVertices.remove( ccsToAvoid.min() );


    ostringstream oss;
    oss << i;
    originalVertices.save( parentDir + "/analysis/" + filename + "-2-" + oss.str() + ".vx", true );
//    originalVertices.save( parentDir + "/analysis/" + filename + "-1-" + oss.str() + ".vx", true );
    EVAL(eqRadii);

    SpatialModelHardcoreDistance3D<float> triMeshSpatialModel;
    triMeshSpatialModel.setRandomGenerator( randomGenerator );
    triMeshSpatialModel.setTriMesh( nucleusTriMesh );
    triMeshSpatialModel.setHardcoreDistances( eqRadii );
    triMeshSpatialModel.initialize();

    VertexStack<float> vertexStack;//( 3, numCCS, 2*numMCSimulations, 0, 0 );
//    Vertices<float> vertices( 3, numCCS-3, 0, 0 );
    Vertices<float> vertices( 3, numCCS-2, 0, 0 );
//    Vertices<float> vertices( 3, numCCS-1, 0, 0 );

    EVAL(vertexStack.getSize());

    for ( int jj = 0; jj < 2*numMCSimulations; ++jj )
    {
      //vertices = triMeshSpatialModel.drawSample(numCCS-3);
      vertices = triMeshSpatialModel.drawSample(numCCS-2);
//      vertices = triMeshSpatialModel.drawSample(numCCS-1);
      vertexStack.insert( jj, vertices );
    }

    const string patternsDir = parentDir + "/patterns/SpatialModelHardcoreDistance3DUsingLessCCS/";
    vertexStack.save( patternsDir + filename + "-2-" + oss.str() + ".vs", true );
//    vertexStack.save( patternsDir + filename + "-1-" + oss.str() + ".vs", true );
  }

}

void evaluatorLessCCS_distanceConstrained(
    const string& filename, const string& parentDir,
    const int& numMCSimulations, RandomGenerator& randomGenerator)
{
  PRINT("spatialModelevaluatorLessCCS_distanceConstrained");

  //open data info
  const string analysisDir = parentDir + "/analysis/";
  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/nuclei/" + filename + ".tm" );

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

  Vector<float> distancesToBorder( numCCS );
  int k = 0;
  for ( int j = lastPos - numCCS + 1 ; j < lastPos + 1; ++j, ++k )
    distancesToBorder[k] = ccsInfo.getValue<float>( "distanceToTheBorder", j );

  EVAL(distancesToBorder);

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

void evaluatorLessCCS_sizeAndDistanceConstrained(
    const string& filename, const string& parentDir,
    const int& numMCSimulations, RandomGenerator& randomGenerator)
{
  PRINT("spatialModelevaluatorLessCCS_sizeAndDistanceConstrained");

  //open data info
  const string analysisDir = parentDir + "/analysis/";
  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/nuclei/" + filename + ".tm" );

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

void evaluatorLessCCS_MaximalRepulsionConstrained(
    const string& filename, const string& parentDir,
    const int& numMCSimulations, RandomGenerator& randomGenerator)
{
  PRINT("spatialModelevaluatorLessCCS_MaximalRepulsionConstrained");

  //open data info
  const string analysisDir = parentDir + "/analysis/";
  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/nuclei/" + filename + ".tm" );

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
  triMeshSpatialModel.setNumMonteCarloCycles( 4000 );
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
void generatePatternsUsingLessObjects(
  const string& filename, const string& parentDir,
  const int& constraints, const int& monteCarloSimulations, RandomGenerator& randomGenerator)
{
  switch( constraints )
  {
    case 0:
      evaluatorLessCCS_completeSpatialRandomness(
        filename, parentDir,
        monteCarloSimulations, randomGenerator );
      break;
    case 1:
      evaluatorLessCCS_sizeConstrained(
        filename, parentDir,
        monteCarloSimulations, randomGenerator );
      break;
    case 2:
      evaluatorLessCCS_distanceConstrained(
        filename, parentDir,
        monteCarloSimulations, randomGenerator );
      break;
    case 3:
      evaluatorLessCCS_sizeAndDistanceConstrained(
        filename, parentDir,
        monteCarloSimulations, randomGenerator );
      break;
    case 4:
      evaluatorLessCCS_MaximalRepulsionConstrained(
        filename, parentDir,
        monteCarloSimulations, randomGenerator );
      break;
  }
}

