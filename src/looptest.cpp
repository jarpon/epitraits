#include <spatialmodelevaluator.h>
#include <spatialdescriptorfunctionf.h>
#include <spatialdescriptorfunctiong.h>
#include <spatialdescriptorfunctionh.h>
#include "trimeshspatialmodel.h"

#define TRACE
#include <trace.h>

void uniformTest(const string&filename, const string& parentDir, DataSet& dataSet)
{
  const string shapesDir = parentDir + "/shapes/";
  //const string analysisDir = parentDir + "/analysis/";
  //const string spatialModelsDir = parentDir + "/spatial_models/";

  const int numVertices = 10;
  const int numPatterns = 99;
  const int numSamples = 100;
  TriMesh<float> nucleusTriMesh( shapesDir + filename + "_nucleus.tm" );

  TrimeshSpatialModel<float> trimeshSpatialModel;
  RandomGenerator randomGenerator;
  trimeshSpatialModel.setRandomGenerator( randomGenerator );
  trimeshSpatialModel.setTriMesh( nucleusTriMesh );
  trimeshSpatialModel.initialize();

  SpatialDescriptorFunctionG<float> spatialDescriptor;

  SpatialModelEvaluator<float,float> modelEvaluator;
  modelEvaluator.setModel( trimeshSpatialModel );
  modelEvaluator.setNumRandomSamples( numPatterns );
  modelEvaluator.setPrecision( 0.05 );
  modelEvaluator.setDescriptor( spatialDescriptor );

  Vector<float> pValues( numSamples );

  for (int i = 0; i < numSamples; ++i)
  {
    EVAL( i );
    Vertices<float> vertices( 3, numVertices );
    vertices = trimeshSpatialModel.drawSample( numVertices );
    pValues[i] = modelEvaluator.eval( vertices );
    dataSet.setValue( "pValues", i, pValues[i] );
    EVAL( pValues[i] );
  }
}
