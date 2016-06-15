#include <spatialmodelcompleterandomness3d.h>
#include <spatialmodelborderdistance3d.h>
#include <spatialmodelhardcoreborderdistance3d.h>
#include <spatialmodelhardcoredistance3d.h>
#include <spatialmodelmaximalrepulsion3d.h>
//#include "spatialmodelmaximalrepulsion3d2.h"
#include <spatialmodel.h>
#include <trimesh.h>
#include <voxelmatrix.h>
#include <fileinfo.h>
#include <dataset.h>
#include <randompartitiongenerator.h>
#include <orbitalhardcoreterritorialspatialmodel.h>
#include <sstream>
#include <iomanip>

#define TRACE
#include <trace.h>
using namespace std;

/*! Preparing constraints for the model
****************************************************************/
void completeSpatialRandomness(
    const string& filename, const string& parentDir,
    const int& numMCSimulations, RandomGenerator& randomGenerator)
{
  PRINT("generatePatterns_completeSpatialRandomness");

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

void sizeConstrained(
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

  Vector<float> eqRadiiTemp( numCCS );
  int k = 0;
  for ( int j = lastPos - numCCS + 1 ; j < lastPos + 1; ++j, ++k )
    eqRadiiTemp[k] = ccsInfo.getValue<float>( "equivalentRadius_ZprojCorrection", j );

  EVAL(eqRadiiTemp);

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

void distanceConstrained(
    const string& filename, const string& parentDir,
    const int& numMCSimulations, RandomGenerator& randomGenerator)
{
  PRINT("spatialModeldistanceConstrained");

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

void sizeAndDistanceConstrained(
    const string& filename, const string& parentDir,
    const int& numMCSimulations, RandomGenerator& randomGenerator)
{
  PRINT("spatialModelsizeAndDistanceConstrained");

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
//    eqRadiiTemp[k] = ccsInfo.getValue<int>( "equivalentRadius_ZprojCorrection", j );
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
  //triMeshSpatialModel.setDistancesToBorder( eqRadii );
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
  //const string patternsDir = parentDir + "/patterns/SpatialModelHardcoreBorderDistance3D-Radii/";
  vertexStack.save( patternsDir + filename + ".vs", true );

}

void maximalRepulsionConstrained(
    const string& filename, const string& parentDir,
    const int& numMCSimulations, RandomGenerator& randomGenerator)
{
  PRINT("spatialModelMaximalRepulsionConstrained");

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
//    eqRadiiTemp[k] = ccsInfo.getValue<int>( "equivalentRadius_ZprojCorrection", j );
//    eqRadiiTemp[k] = ccsInfo.getValue<float>( "equivalentRadius_PSFVolCorrection", j );
    distancesToBorder[k] = ccsInfo.getValue<int>( "distanceToTheBorder", j );
  }


  for ( int i = 0; i < eqRadiiTemp.getSize(); ++i )
    if ( eqRadiiTemp[i] > distancesToBorder[i] )
      eqRadiiTemp[i] = distancesToBorder[i];

  Vector<float> eqRadii = eqRadiiTemp;
//  eqRadii.setZeros();
  EVAL(eqRadii);
//  EVAL(distancesToBorder);


//  SpatialModelMaximalRepulsion3D2 <float> triMeshSpatialModel;
  SpatialModelMaximalRepulsion3D <float> triMeshSpatialModel;
  randomGenerator.init( 101 );
  triMeshSpatialModel.setRandomGenerator( randomGenerator );
  triMeshSpatialModel.setTriMesh( nucleusTriMesh );
  triMeshSpatialModel.setNumMonteCarloCycles( 100000 );
  triMeshSpatialModel.setHardcoreDistances( eqRadii );
  triMeshSpatialModel.initialize();
  triMeshSpatialModel.initializeBeta( numCCS );
  DataSet energyProfile;
//  energyProfile.setValues<float> ( "energyProfile", triMeshSpatialModel.getEnergyProfile() );
//  energyProfile.save( filename + ".data", true );

  VertexStack<float> vertexStack;//( 3, numCCS, 2*numMCSimulations, 0, 0 );
  Vertices<float> vertices( 3, numCCS, 0, 0 );

  EVAL(vertexStack.getSize());

  for ( int jj = 0; jj < numMCSimulations; ++jj )
  //for ( int jj = 0; jj < 2 * numMCSimulations; ++jj )
  {
    vertices = triMeshSpatialModel.drawSample(numCCS);
    vertexStack.insert( jj, vertices );
  }

  energyProfile.setValues<float> ( "energyProfile", triMeshSpatialModel.getEnergyProfile() );
  energyProfile.save( filename + ".data", true );
  const string patternsDir = parentDir + "/patterns/SpatialModelMaximalRepulsion3D/";
  vertexStack.save( patternsDir + filename + ".vs", true );

}


void divideIntoTerritories(
    const string& filename, const string& parentDir,
    const int& numMCSimulations, RandomGenerator& randomGenerator)
{
  PRINT("divideIntoTerritories");

  //open data info
  const string analysisDir = parentDir + "/analysis/";
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

  const string nucleiDir = parentDir + "/segmented_nuclei/";
  VoxelMatrix<float> originalVoxelMatrixDomain;
  originalVoxelMatrixDomain.load( nucleiDir + filename + ".vm" );
  VoxelMatrix<float> voxelMatrixCT;

  RandomPartitionGenerator<float> randomPartitionGenerator;
  //randomPartitionGenerator.setRandomGenerator( *(new RandomGenerator()) );
  randomPartitionGenerator.setRandomGenerator( randomGenerator );
  randomPartitionGenerator.setNumRegions( numCCS );

  for ( int jj = 0; jj < numMCSimulations; ++jj )
  {
    ostringstream oss;
    oss << setw(2) << setfill('0') << jj +1;
    voxelMatrixCT = originalVoxelMatrixDomain;
    randomPartitionGenerator.run( voxelMatrixCT );
    voxelMatrixCT.save( parentDir + filename + "-territories-" + oss.str() + ".vm", true );
  }

}


  void sizeConstrainedIntoTerritories(
      const string& filename, const string& parentDir,
      const int& numMCSimulations, RandomGenerator& randomGenerator)
  {
    PRINT("divideIntoTerritories");

    //open data info
    const string analysisDir = parentDir + "/analysis/";
    const TriMesh<float> territoriesTriMesh ( parentDir + "/shapes/territories/" + filename + ".tm" );

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

    TriMesh<float> currentTerritory;
    Vector< TriMesh<float> >* allTerritoriesTrimeshes;
    allTerritoriesTrimeshes->setSize( numCCS );

    for ( int j = lastPos - numCCS + 1 ; j < lastPos + 1; ++j, ++k )
    {
      eqRadiiTemp[k] = ccsInfo.getValue<float>( "equivalentRadius_ZprojCorrection", j );
  //    eqRadiiTemp[k] = ccsInfo.getValue<int>( "equivalentRadius_ZprojCorrection", j );
  //    eqRadiiTemp[k] = ccsInfo.getValue<float>( "equivalentRadius_PSFVolCorrection", j );
      distancesToBorder[k] = ccsInfo.getValue<int>( "distanceToTheBorder", j );
      ostringstream oss;
      oss << setw(2) << setfill('0') << j +1;
      currentTerritory.load( parentDir + "/shapes/territories/" + filename + "-territories-" + oss.str() + ".tm" );
      allTerritoriesTrimeshes[k] = currentTerritory;
    }


    for ( int i = 0; i < eqRadiiTemp.getSize(); ++i )
      if ( eqRadiiTemp[i] > distancesToBorder[i] )
        eqRadiiTemp[i] = distancesToBorder[i];

    Vector<float> eqRadii = eqRadiiTemp;
  //  eqRadii.setZeros();
    EVAL(eqRadii);
  //  EVAL(distancesToBorder);


    OrbitalHardcoreTerritorialSpatialModel<float> orbitalHardcoreTerritorialSpatialModel;
    orbitalHardcoreTerritorialSpatialModel.setRandomGenerator( randomGenerator );
    orbitalHardcoreTerritorialSpatialModel.setHardcoreDistances( eqRadii );
    //orbitalHardcoreTerritorialSpatialModel.setDistancesToBorder( hardcoreDistances );
    orbitalHardcoreTerritorialSpatialModel.setTerritories( *allTerritoriesTrimeshes );
    orbitalHardcoreTerritorialSpatialModel.setTriMesh( territoriesTriMesh );
    orbitalHardcoreTerritorialSpatialModel.initialize();

    VertexStack<float> vertexStack;
    Vertices<float> vertices( 3, numCCS, 0, 0 );

    for ( int jj = 0; jj < numMCSimulations; ++jj )
    //for ( int jj = 0; jj < 2 * numMCSimulations; ++jj )
    {
      vertices = orbitalHardcoreTerritorialSpatialModel.drawSample( numCCS );
      vertexStack.insert( jj, vertices );
    }

    //const string patternsDir = parentDir + "/patterns/SpatialModelTerritories/";
    vertexStack.save( parentDir + filename + ".vs", true );

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
      completeSpatialRandomness(
        filename, parentDir,
        monteCarloSimulations, randomGenerator );
      break;
    case 1:
      sizeConstrained(
        filename, parentDir,
        monteCarloSimulations, randomGenerator );
      break;
    case 2:
      distanceConstrained(
        filename, parentDir,
        monteCarloSimulations, randomGenerator );
      break;
    case 3:
      sizeAndDistanceConstrained(
        filename, parentDir,
        monteCarloSimulations, randomGenerator );
      break;
    case 4:
      maximalRepulsionConstrained(
        filename, parentDir,
        monteCarloSimulations, randomGenerator );
      break;
    case 5:
      divideIntoTerritories(
        filename, parentDir,
        monteCarloSimulations, randomGenerator );
      break;
    case 6:
      sizeConstrainedIntoTerritories(
        filename, parentDir,
        monteCarloSimulations, randomGenerator );
      break;
  }
}
