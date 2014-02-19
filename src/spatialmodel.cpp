#include <randomgenerator.h>
#include "trimeshspatialmodel.h"
#include <trimesh.h>
#include <voxelmatrix.h>
#include <fileinfo.h>
//#include "regionanalysis2.h"
#include <regionanalysis.h>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <string.h>


#define TRACE
#include <trace.h>
using namespace std;


void spatialModelAnalysis(VoxelMatrix<float>& chromocentersVoxelMatrix, TriMesh<float>& nucleusTriMesh,
                          const string& filename,const string& parentDir, const int& numPatterns)
{
  const string shapesDir = parentDir + "/shapes/";
  const string spatialModelsDir = parentDir + "/spatial_models/";

  RegionAnalysis<float> regionAnalysisCCs;
  regionAnalysisCCs.setRegionMatrix( chromocentersVoxelMatrix );
  regionAnalysisCCs.run();

  Vector<float> vertex(3);
  Vector<float> vertexTriMesh(3);
  string currentFilename;
  int labels = regionAnalysisCCs.numRegions();
  Vertices<float> centroidVertices ( 3, labels, 0, 0 );
  Vector<float> distancesToBorderVector ( labels );
  Vector<float> eqRadiusVector ( labels );

  for (int numCC = 1; numCC <= labels; ++numCC )
  {
    ostringstream oss; //we suppose as much 99 labels
    oss << setw(2) << setfill('0') << numCC;
    currentFilename = shapesDir + filename + "_chromocenters-" + oss.str() + ".tm";
    TriMesh<float> currentChromocenterTriMesh ( currentFilename );

    eqRadiusVector[numCC-1] = currentChromocenterTriMesh.equivalentRadius();

    vertex = currentChromocenterTriMesh.cog();
    centroidVertices[numCC-1] = vertex;
    nucleusTriMesh.closestPoint(vertex,vertexTriMesh);
    distancesToBorderVector[numCC-1] = vertex.distance(vertexTriMesh);
  }


  TrimeshSpatialModel<float> trimeshSpatialModel;
  RandomGenerator randomGenerator;
  trimeshSpatialModel.setRandomGenerator( randomGenerator );
  trimeshSpatialModel.setTriMesh( nucleusTriMesh );
  trimeshSpatialModel.initialize();

  trimeshSpatialModel.setDistanceToBorder( distancesToBorderVector );
  trimeshSpatialModel.setHardcoreDistance( eqRadiusVector );

  Vertices<float> generatedResults;

  for ( int i = 1; i <= numPatterns; ++i )
  {
    ostringstream iss;
    generatedResults = trimeshSpatialModel.drawSample( labels );
    iss << setw(2) << setfill('0') << i;
    generatedResults.save( spatialModelsDir + filename + "_chromocenters-" + iss.str() + ".vx", true  );
  }

}


