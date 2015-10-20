//#include <randomgenerator.h>
//#include "trimeshspatialmodel.h"
//#include <trimesh.h>
//#include <voxelmatrix.h>
//#include <fileinfo.h>
////#include <regionanalysis.h>
//#include <sstream>
//#include <iomanip>
//#include <stdio.h>
//#include <string.h>
//#include <dataset.h>

//#define TRACE
//#include <trace.h>
//using namespace std;


//void spatialModelAnalysis(TriMesh<float>& nucleusTriMesh,
//                          const string& filename, const string& parentDir, const int& numPatterns)
//{
//  PRINT("spatialModels generator");

//  const string spatialModelsDir = parentDir + "/spatial_models/";
//  const string analysisDir = parentDir + "/analysis/";

//  const DataSet datasetNucleus( analysisDir + filename + "_chromocenters.csv" );
////  const DataSet datasetNucleus( analysisDir + filename + ".csv" );

//  const int numPoints = datasetNucleus.size()[0];
//  const Vector<float> eqRadii = datasetNucleus.getValues<float>( "equivalentRadius_tm" );
//  //const Vector<float> eqRadii = datasetNucleus.getValues<float>( "chromocenterRadius" );
//  //const Vector<float> distancesToBorder = datasetNucleus.getValues<float>( "distanceToTheBorder" );

//  Vertices<float> vertices ( 3, numPoints, 0, 0 );
//  for ( int i = 0; i < numPoints; ++i )
//  {
//    vertices[i][0] = datasetNucleus.getValue<float>( "centroidCoordX", i );
//    vertices[i][1] = datasetNucleus.getValue<float>( "centroidCoordY", i );
//    vertices[i][2] = datasetNucleus.getValue<float>( "centroidCoordZ", i );
//  }

//  TriMeshSpatialModel<float> triMeshSpatialModel;
//  RandomGenerator randomGenerator;
//  triMeshSpatialModel.setRandomGenerator( randomGenerator );
//  triMeshSpatialModel.setTriMesh( nucleusTriMesh );
//  triMeshSpatialModel.initialize();

////  triMeshSpatialModel.setDistanceToBorder( distancesToBorder );
//  triMeshSpatialModel.setHardcoreDistances( eqRadii );

//  Vertices<float> generatedResults;

//  //for ( int i = 1; i <= numPatterns; ++i )
//  for ( int i = 1; i <= 1; ++i )
//  {
//    ostringstream iss;
//    generatedResults = triMeshSpatialModel.drawSample( numPoints );
//    iss << setw(2) << setfill('0') << i;
//    //generatedResults.save( spatialModelsDir + filename + "_chromocenters-" + iss.str() + ".vx", true  );
//    generatedResults.save( parentDir + "/" + filename + ".vx", true  );
//  }

//}


