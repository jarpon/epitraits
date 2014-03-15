#include <randomgenerator.h>
#include "trimeshspatialmodel.h"
#include <trimesh.h>
#include <voxelmatrix.h>
#include <fileinfo.h>
//#include <regionanalysis.h>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <dataset.h>

#define TRACE
#include <trace.h>
using namespace std;


void spatialModelAnalysis(TriMesh<float>& nucleusTriMesh,
                          const string& filename, const string& parentDir, const int& numPatterns)
{
  const string spatialModelsDir = parentDir + "/spatial_models/";
  const string analysisDir = parentDir + "/analysis/";

  DataSet datasetNucleus;
  datasetNucleus.load( analysisDir + filename + "_chromocenters.csv" );
  int labels = datasetNucleus.size()[0];

  Vector<float> vertexTriMesh(3);
  Vertices<float> centroidVertices ( 3, labels, 0, 0 );
  Vector<float> distancesToBorderVector ( labels );

  Vector<float> eqRadiusVector = datasetNucleus.getValues<float>( "equivalentRadius" );
  centroidVertices[X] = datasetNucleus.getValues<float>( "centroidCoordX" );
  centroidVertices[Y] = datasetNucleus.getValues<float>( "centroidCoordY" );
  centroidVertices[Z] = datasetNucleus.getValues<float>( "centroidCoordZ" );

  for (int numCC = 0; numCC < labels; ++numCC )
  {
    nucleusTriMesh.closestPoint( centroidVertices[numCC], vertexTriMesh );
    distancesToBorderVector[numCC-1] = centroidVertices[numCC].distance( vertexTriMesh );
  }


  TriMeshSpatialModel<float> triMeshSpatialModel;
  RandomGenerator randomGenerator;
  triMeshSpatialModel.setRandomGenerator( randomGenerator );
  triMeshSpatialModel.setTriMesh( nucleusTriMesh );
  triMeshSpatialModel.initialize();

  triMeshSpatialModel.setDistanceToBorder( distancesToBorderVector );
  triMeshSpatialModel.setHardcoreDistance( eqRadiusVector );

  Vertices<float> generatedResults;

  for ( int i = 1; i <= numPatterns; ++i )
  {
    ostringstream iss;
    generatedResults = triMeshSpatialModel.drawSample( labels );
    iss << setw(2) << setfill('0') << i;
    generatedResults.save( spatialModelsDir + filename + "_chromocenters-" + iss.str() + ".vx", true  );
  }

}


