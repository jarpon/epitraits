#include <spatialmodelevaluator.h>
#include <spatialdescriptorfunctionf.h>
#include <spatialdescriptorfunctiong.h>
#include <spatialdescriptorfunctionh.h>
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


void spatialModelEvaluator(VoxelMatrix<float>& chromocentersVoxelMatrix, TriMesh<float>& nucleusTriMesh,
                          const string& filename, const string& parentDir,
                          const int& numPatterns, const string& function,
                          DataSet& dataSet)
{
  const string shapesDir = parentDir + "/shapes/";

  RegionAnalysis<float> regionAnalysisCCs;
  regionAnalysisCCs.setRegionMatrix( chromocentersVoxelMatrix );
  regionAnalysisCCs.run();

  Vector<float> vertex(3);
  Vector<float> vertexTriMesh(3);
  string currentFilename;
  int labels = regionAnalysisCCs.numRegions();
  Vertices<float> centroidChromocenters ( 3, labels, 0, 0 );
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
    centroidChromocenters[numCC-1] = vertex;
    nucleusTriMesh.closestPoint(vertex,vertexTriMesh);
    distancesToBorderVector[numCC-1] = vertex.distance(vertexTriMesh);
  }

  EVAL(labels);
  //centroidChromocenters.save( shapesDir + filename + "_chromocenters.vx", true );


  TrimeshSpatialModel<float> trimeshSpatialModel;
  RandomGenerator randomGenerator;
  trimeshSpatialModel.setRandomGenerator( randomGenerator );
  trimeshSpatialModel.setTriMesh( nucleusTriMesh );
  //trimeshSpatialModel.setDistanceToBorder( distancesToBorderVector );
  //trimeshSpatialModel.setHardcoreDistance( eqRadiusVector );
  trimeshSpatialModel.initialize();

  SpatialModelEvaluator<float,float> modelEvaluator;
  modelEvaluator.setModel( trimeshSpatialModel );
  modelEvaluator.setNumRandomSamples( numPatterns );
  modelEvaluator.setPrecision( 0.05 );

  SpatialDescriptor<float>* spatialDescriptor;

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
  else //if ( function == "F" )
  {
    PRINT("F");
    TrimeshSpatialModel<float> tempSpatialModel;
    tempSpatialModel.setRandomGenerator( randomGenerator );
    tempSpatialModel.setTriMesh( nucleusTriMesh );
    tempSpatialModel.initialize();

    SpatialDescriptorFunctionF<float>* spatialDescriptorFunctionF;
    spatialDescriptorFunctionF = new SpatialDescriptorFunctionF<float>();
    spatialDescriptorFunctionF->setEvaluationPositions( tempSpatialModel.drawSample(10) );
    spatialDescriptor = spatialDescriptorFunctionF;
  }
  modelEvaluator.setDescriptor( *spatialDescriptor );

  PRINT( "let's go");
  const float pValue = modelEvaluator.eval( centroidChromocenters, &dataSet );
  EVAL(pValue);

}
