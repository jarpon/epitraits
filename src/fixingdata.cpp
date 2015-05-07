#include <componentlabelling.h>
#include <dataset.h>
//#include "regionanalysis2.h"
#include <regionanalysis.h>
#include <cmath>
#include <dataset.h>
#include <thresholding.h>
#include <trimesh.h>
#include "trimeshspatialmodel.h"
#include "maximalrepulsion.h"
#include <curvestack.h>
#include <voxelmatrixdilatation.h>
#include <voxelmatrixerosion.h>

#define TRACE
#include <trace.h>


//void doIt2(const string& filename, const string& parentDir)
//{
//  TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + "_nucleus.tm" );
//  //TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + "-nucleus.tm" );

//  const string analysisDir = parentDir + "/analysis/";
//  EVAL(analysisDir);

//  Vector<float> centroid(3);
//  Vector<float> vertexTriMesh(3);
//  EVAL(analysisDir + filename + "_chromocenters.csv");
//  DataSet datasetNucleus( analysisDir + filename + "_chromocenters.csv" );
//  EVAL("ici");
//  const Vector<float> ids = datasetNucleus.getValues<float>( "idCC" );


///**///chromocenters individual information
//  for (int numCC = 0; numCC < ids.getSize(); numCC++ )
//  {
//    centroid[X] = datasetNucleus.getValue<float>( "centroidCoordX", numCC);
//    centroid[Y] = datasetNucleus.getValue<float>( "centroidCoordY", numCC);
//    centroid[Z] = datasetNucleus.getValue<float>( "centroidCoordZ", numCC);

//    float distanceToBorder = nucleusTriMesh.closestPoint( centroid, vertexTriMesh );

//    datasetNucleus.setValue("distanceToTheBorder", numCC, distanceToBorder);
//  }

//  datasetNucleus.save( analysisDir + filename + ".csv" );
//  EVAL("done!");
//}

void doIt2(const string& filename, const string& parentDir)
{
  VoxelMatrix<float> nucleusMask( parentDir + "/segmented_nuclei/" + filename + ".vm" );



  VoxelMatrix<float> structElement;
  structElement.setSize(3,3,3);
  structElement.setOnes();

//  VoxelMatrixDilatation<float> voxelDilatation;
//  voxelDilatation.setStructElt( structElement );
//  voxelDilatation.apply( nucleusMask );

  VoxelMatrixErosion<float> voxelErosion;
  voxelErosion.setStructElt( structElement );
  voxelErosion.apply( nucleusMask );


  Thresholding<float> thresholding;
  thresholding.setBackground( 0.0 );
  thresholding.setForeground( 1.0 );
  thresholding.setThreshold( -0.1 );
  thresholding.apply( nucleusMask );


  nucleusMask.save ( parentDir + "/" + filename + ".vm", true );





//  for (int k = 0; k < sizeZ; ++k)  holesFilling.apply( ccsMask[k] );

  //labeling the image
}


void doIt(const string& filename, const string& parentDir, RandomGenerator& randomGenerator)
{
  PRINT("maximaRepulsion");

  //const DataSet datasetNucleus( parentDir + "/analysis/" + filename + "_chromocenters.csv" );
  const DataSet datasetNucleus( parentDir + "/analysis/" + filename + ".csv" );
  const int numPoints = datasetNucleus.size()[0];

  //const Vector<float> eqRadii = datasetNucleus.getValues<float>( "equivalentRadius_tm" );
  const Vector<float> eqRadii = datasetNucleus.getValues<float>( "chromocenterRadius" );

  //const int numPoints = 5;

  Vertices<float> vertices ( 3, numPoints, 0, 0 );
//  for ( int i = 0; i < numPoints; ++i )
//  {
//    vertices[i][0] = datasetNucleus.getValue<float>( "centroidCoordX", i );
//    vertices[i][1] = datasetNucleus.getValue<float>( "centroidCoordY", i );
//    vertices[i][2] = datasetNucleus.getValue<float>( "centroidCoordZ", i );
//    EVAL(vertices[i]);
//  }

//  Vertices<float> vertices ( 3, 4, 0, 0 );
//  vertices[1][0] = 75;
//  vertices[1][1] = 75;
//  vertices[1][2] = 75;
//  vertices[2][0] = 77;
//  vertices[2][1] = 75;
//  vertices[2][2] = 75;
//  vertices[3][0] = 75;
//  vertices[3][1] = 75;
//  vertices[3][2] = 77;
//  vertices[4][0] = 77;
//  vertices[4][1] = 75;
//  vertices[4][2] = 77;


  Vertices<float> repulsiveVertices ( 3, numPoints, 0, 0 );

  //const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + "_nucleus.tm" );
  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + "-nucleus.tm" );
  //const TriMesh<float> nucleusTriMesh ( "/home/jarpon/Desktop/sphere.tm" );
//  TriMeshSpatialModel<float> tempTriMeshSpatialModel;
//  tempTriMeshSpatialModel.setRandomGenerator( randomGenerator );
//  tempTriMeshSpatialModel.setHardcoreDistances( eqRadii );
//  tempTriMeshSpatialModel.setTriMesh( nucleusTriMesh );
//  tempTriMeshSpatialModel.initialize();
//  vertices = tempTriMeshSpatialModel.drawSample( numPoints );
//  vertices.save( "/home/jarpon/Desktop/maxRepulsion_begin.vx", true );
  for ( int j = 0; j < numPoints; ++j )
  {
    EVAL(vertices[j]);
  }

  MaximalRepulsionTriMeshSpatialModel<float> triMeshSpatialModel;
  triMeshSpatialModel.setRandomGenerator( randomGenerator );
  triMeshSpatialModel.setTriMesh( nucleusTriMesh );
  triMeshSpatialModel.setHardcoreDistance( eqRadii );
  triMeshSpatialModel.initialize();


  //maximalRepulsion.setTriMeshSpatialModel( tempTriMeshSpatialModel );
  //maximalRepulsion.setTriMesh( nucleusTriMesh );
  //maximalRepulsion.setMaximaRepulsion();
  repulsiveVertices = triMeshSpatialModel.drawSample( numPoints );

}
