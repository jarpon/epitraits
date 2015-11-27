#include "spatialmodelhardcoredistance3ddifferentcompartments.h"
#include "spatialmodelborderdistance3ddifferentcompartments.h"
#include "spatialmodelhardcoreborderdistance3ddifferentcompartments.h"

#include <dataset.h>
#include <fileinfo.h>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <cmath>
#include <iostream>
#include <string>

#include <stringtools.h>
#include <fstream>
#include <sstream>
#include <unistd.h>
 #include <cstdlib>


#define TRACE
#include <trace.h>
using namespace std;


void evalDiffCompartments(
  const TriMesh<float>& nucleusTriMesh,
  TriMeshSpatialModelDifferentCompartments<float>& triMeshSpatialModelDifferentCompartments,
  const string& filename, const string& parentDir,
  const string& function, const int& constraints,
  DataSet& dataSet, RandomGenerator& randomGenerator)
{
  string classif = parentDir;
  classif = classif.substr( classif.find_last_of("/\\")+1, classif.length() );
  //open data info
  const string analysisDir = parentDir + "/analysis/";
//  const DataSet datasetNucleus( analysisDir + filename + "_chromocenters.csv" );
//  const DataSet datasetNucleus( analysisDir + filename + "_nucleoli.csv" );
  DataSet globalAnalysis( analysisDir + "nuclei.csv" );
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


  Vertices<float> vertices( 3, numCCS, 0, 0 );
  int k = 0;
  for ( int j = lastPos - numCCS + 1 ; j < lastPos + 1; ++j, ++k )
  {
    vertices[k][0] = ccsInfo.getValue<float>( "centroidCoordX", j );
    vertices[k][1] = ccsInfo.getValue<float>( "centroidCoordY", j );
    vertices[k][2] = ccsInfo.getValue<float>( "centroidCoordZ", j );
    EVAL(vertices[k]);
  }

  Vector<string> nucleiNames ;
  nucleiNames = globalAnalysis.getValues<string>( globalAnalysis.variableNames()[0] );

  //unifying datasets
  int numCurrentNucleus;

  for ( int j = 0; j < nucleiNames.getSize(); ++j )
    if ( nucleiNames[j] == filename )
      numCurrentNucleus = j;

  EVAL (numCurrentNucleus);

  //const int numPoints = vertices.getNumVertices();
  const int numPatterns = 99;

//  SpatialModelEvaluator<float,float> modelEvaluator;
//  modelEvaluator.setModel( triMeshSpatialModelDifferentCompartments );
//  modelEvaluator.setNumMonteCarloSamples( numPatterns ); //to check uniformity
//  modelEvaluator.setPrecision( 0.05 );

//  SpatialDescriptor<float>* spatialDescriptor;


//  //setting function parameters
//  if ( function == "all" )
//  {
//    PRINT("all functions");

//    SpatialModelCompleteRandomness3D<float> tempTriMeshSpatialModelDifferentCompartments;
//    tempTriMeshSpatialModelDifferentCompartments.setRandomGenerator( randomGenerator );
//    tempTriMeshSpatialModelDifferentCompartments.setTriMesh( nucleusTriMesh );
//    tempTriMeshSpatialModelDifferentCompartments.initialize();
//    Vertices<float> evaluationPositions = tempTriMeshSpatialModelDifferentCompartments.drawSample( 10000 );
//    //evaluationPositions.save( parentDir + "/" + filename + "_Fpattern-" + ".vx", true); //to check uniformity of the F patterns

//    SpatialDescriptorFunctionF<float>* spatialDescriptorFunctionF;
//    spatialDescriptorFunctionF = new SpatialDescriptorFunctionF<float>();
//    spatialDescriptorFunctionF->setEvaluationPositions( evaluationPositions );
//    spatialDescriptor = spatialDescriptorFunctionF;
//    modelEvaluator.addDescriptor( *spatialDescriptor );

//    spatialDescriptor = new SpatialDescriptorFunctionG<float>();
//    modelEvaluator.addDescriptor( *spatialDescriptor );

//    spatialDescriptor = new SpatialDescriptorFunctionH<float>();
//    modelEvaluator.addDescriptor( *spatialDescriptor );

//    SpatialDescriptorFunctionB<float>* spatialDescriptorFunctionB;
//    spatialDescriptorFunctionB = new SpatialDescriptorFunctionB<float>();
//    spatialDescriptorFunctionB->setTriMesh( nucleusTriMesh );
//    spatialDescriptor = spatialDescriptorFunctionB;
//    modelEvaluator.addDescriptor( *spatialDescriptor );

//    SpatialDescriptorFunctionC<float>* spatialDescriptorFunctionC;
//    spatialDescriptorFunctionC = new SpatialDescriptorFunctionC<float>();
//    spatialDescriptorFunctionC->setCenter( nucleusTriMesh.cog() );
//    spatialDescriptor = spatialDescriptorFunctionC;
//    modelEvaluator.addDescriptor( *spatialDescriptor );
//  }
//  else if ( function == "G" )
//  {
//    PRINT("G");
//    spatialDescriptor = new SpatialDescriptorFunctionG<float>();
//  }
//  else if ( function == "H" )
//  {
//    PRINT("H");
//    spatialDescriptor = new SpatialDescriptorFunctionH<float>();
//  }
//  else if ( function == "B" )
//  {
//    PRINT("B");
//    SpatialDescriptorFunctionB<float>* spatialDescriptorFunctionB;
//    spatialDescriptorFunctionB = new SpatialDescriptorFunctionB<float>();
//    spatialDescriptorFunctionB->setTriMesh( nucleusTriMesh );
//    spatialDescriptor = spatialDescriptorFunctionB;
//  }
//  else if ( function == "C" )
//  {
//    PRINT("C");
//    SpatialDescriptorFunctionC<float>* spatialDescriptorFunctionC;
//    spatialDescriptorFunctionC = new SpatialDescriptorFunctionC<float>();
//    spatialDescriptorFunctionC->setCenter( nucleusTriMesh.cog() );
//    spatialDescriptor = spatialDescriptorFunctionC;
//  }
//  else //if ( function == "F" )
//  {
//    PRINT("F");
//    SpatialModelCompleteRandomness3D<float> tempTriMeshSpatialModelDifferentCompartments;
//    tempTriMeshSpatialModelDifferentCompartments.setRandomGenerator( randomGenerator );
//    tempTriMeshSpatialModelDifferentCompartments.setTriMesh( nucleusTriMesh );
//    tempTriMeshSpatialModelDifferentCompartments.initialize();
//    Vertices<float> evaluationPositions = tempTriMeshSpatialModelDifferentCompartments.drawSample( 10000 );
//    //evaluationPositions.save( parentDir + "/spatial_models/" + filename + "_Fpattern-" + ".vx", true);
//    SpatialDescriptorFunctionF<float>* spatialDescriptorFunctionF;
//    spatialDescriptorFunctionF = new SpatialDescriptorFunctionF<float>();
//    spatialDescriptorFunctionF->setEvaluationPositions( evaluationPositions );
//    spatialDescriptor = spatialDescriptorFunctionF;
//  }


//  //processing data
//  if ( function != "all" )
//  {
//    string newVariable;
//    switch ( constraints )
//    {
//      case 0:
//        newVariable = ( "SpatialModelCompleteRandomness3D_" + function + "-SDI" );
//        break;
//      case 1:
//        newVariable = ( "SpatialModelHardcoreDistance3D_" + function + "-SDI" );
//        break;
//      case 2:
//        newVariable = ( "spatialModelBorderDistance3D_" + function + "-SDI" );
//        break;
//      case 3:
//        newVariable = ( "SpatialModelBorderHardcoreDistance3D_" + function + "-SDI" );
//        break;
//      case 4:
//        newVariable = ( "SpatialModelMaximalRepulsion3D_" + function + "-SDI" );
//        break;
//    }

//    //globalAnalysis.setValues<float>( newVariable, -1 );
//    float sdi;
//    int row;

//    modelEvaluator.setDescriptor( *spatialDescriptor );

//    ostringstream iss; //we suppose as much 99 labels
//    iss << constraints;
//    DataSet saveTest;
//    try {
//      sdi = modelEvaluator.eval( vertices, &saveTest );

//  //    Vector<float> output = modelevalDiffCompartments.evalSDIandMaxDiff( vertices, &saveTest );
//  //    EVAL( output[0] );
//  //    EVAL( output[1] );


//      saveTest.save( analysisDir + iss.str() + "/" + function + "/" + filename + ".csv", true );
//    //  saveTest.save( analysisDir + iss.str() + "/" + function + "/" + filename + "_random.csv", true );
//      row = dataSet.size()[0];

//      dataSet.setValue( "nucleus", row, filename );
//      dataSet.setValue( "class", row, classif );//classification: mutant, tissue, etc.
//      dataSet.setValue( "descriptor", row, function );//spatial descriptor
//      dataSet.setValue( "index", row, sdi );
//      //dataSet.setValue( "index", row, output[0] );
//      //dataSet.setValue( "signedMaxDiff", row, output[1] );

//      globalAnalysis.setValue( newVariable, numCurrentNucleus, sdi );
//    }
//    catch( Exception exception ) {
//      EVAL( exception.getWhat() );
//      EVAL( exception.getWhere() );

//      const int row = dataSet.size()[0];

//      dataSet.setValue( "nucleus", row, filename );
//      dataSet.setValue( "class", row, classif );//classification: mutant, tissue, etc.
//      dataSet.setValue( "descriptor", row, function );//spatial descriptor
//      dataSet.setValue( "index", row, sqrt(-1) );


//    }

//  }
//  else
//  {
////    Vertices<float> vertices ( 3, numPoints, 0, 0 );
////    for ( int i = 0; i < numPoints; ++i )
////    {
////      vertices[i][0] = datasetNucleus.getValue<float>( "centroidCoordX", i );
////      vertices[i][1] = datasetNucleus.getValue<float>( "centroidCoordY", i );
////      vertices[i][2] = datasetNucleus.getValue<float>( "centroidCoordZ", i );
////      EVAL(vertices[i]);
////    }

//    const int row = dataSet.numRows();
//    vector<float> sdis;
//    vector<int> ranks;
////    vector<float> maxDiff;

//    try {
////    modelevalDiffCompartments.evalSDIandMaxDiff( vertices, pValues, ranks, maxDiff);
//    modelEvaluator.eval( vertices, sdis, ranks);

//    EVAL( sdis[0] );
//    EVAL( sdis[1] );
//    EVAL( sdis[2] );
//    EVAL( sdis[3] );
//    EVAL( sdis[4] );

//    string newVariable;
//    switch ( constraints )
//    {
//      case 0: newVariable = "SpatialModelCompleteRandomness3D_";
//      case 1: newVariable = "SpatialModelHardcoreDistance3D_";
//      case 2: newVariable = "spatialModelBorderDistance3D_";
//      case 3: newVariable = "SpatialModelBorderHardcoreDistance3D_";
//      case 4: newVariable = "SpatialModelMaximalRepulsion3D_";
//    }

//    //unifying datasets
//    for ( int jj = 0; jj < sdis.size(); ++jj )
//    {
//      if ( ( sdis[jj] < 0 || sdis[jj] > 1 )  && ( jj = 0 ) )
//        globalAnalysis.setValue( newVariable + "F-SDI", numCurrentNucleus, sqrt(-1) );
//      else
//        globalAnalysis.setValue( newVariable + "F-SDI", numCurrentNucleus, sdis[jj] );
//      if ( ( sdis[jj] < 0 || sdis[jj] > 1 )  && ( jj = 1 ) )
//        globalAnalysis.setValue( newVariable + "G-SDI", numCurrentNucleus, sqrt(-1) );
//      else
//        globalAnalysis.setValue( newVariable + "G-SDI", numCurrentNucleus, sdis[jj] );
//      if ( ( sdis[jj] < 0 || sdis[jj] > 1 )  && ( jj = 2 ) )
//        globalAnalysis.setValue( newVariable + "H-SDI", numCurrentNucleus, sqrt(-1) );
//      else
//        globalAnalysis.setValue( newVariable + "H-SDI", numCurrentNucleus, sdis[jj] );
//      if ( ( sdis[jj] < 0 || sdis[jj] > 1 )  && ( jj = 3 ) )
//        globalAnalysis.setValue( newVariable + "B-SDI", numCurrentNucleus, sqrt(-1) );
//      else
//        globalAnalysis.setValue( newVariable + "B-SDI", numCurrentNucleus, sdis[jj] );
//      if ( ( sdis[jj] < 0 || sdis[jj] > 1 )  && ( jj = 4 ) )
//        globalAnalysis.setValue( newVariable + "C-SDI", numCurrentNucleus, sqrt(-1) );
//      else
//        globalAnalysis.setValue( newVariable + "C-SDI", numCurrentNucleus, sdis[jj] );
//    }

//    dataSet.setValue( "nucleus", row, filename );
//    dataSet.setValue( "class", row, classif );//classification: mutant, tissue, etc.
//    dataSet.setValue( "F-SDI", row, sdis[0] );
////    dataSet.setValue( "F-maxDiff", row, maxDiff[0] );
//    dataSet.setValue( "G-SDI", row, sdis[1] );
////    dataSet.setValue( "G-maxDiff", row, maxDiff[1] );
//    dataSet.setValue( "H-SDI", row, sdis[2] );
////    dataSet.setValue( "H-maxDiff", row, maxDiff[2] );
//    dataSet.setValue( "B-SDI", row, sdis[3] );
////    dataSet.setValue( "B-maxDiff", row, maxDiff[3] );
//    dataSet.setValue( "C-SDI", row, sdis[4] );
////    dataSet.setValue( "C-maxDiff", row, maxDiff[4] );

//    }
//    catch( Exception exception ) {
//      EVAL( exception.getWhat() );
//      EVAL( exception.getWhere() );
//      const int row = dataSet.size()[0];

//      dataSet.setValue( "nucleus", row, filename );
//      dataSet.setValue( "class", row, classif );//classification: mutant, tissue, etc.
//      dataSet.setValue( "descriptor", row, function );//spatial descriptor
//      dataSet.setValue( "index", row, sqrt(-1) );
//      //dataSet.setValue( "index", row, output[0] );
//      //dataSet.setValue( "signedMaxDiff", row, output[1] );

//    }

//  }

//  globalAnalysis.save( analysisDir + "/" + "nuclei_complete2.csv", true );


}


/*! Preparing constraints for the model
****************************************************************/
void evalDiffCompartments_completeSpatialRandomness(
  const string& filename, const string& parentDir,
  const string& function, const int& constraints,
  DataSet& dataSet,
  RandomGenerator& randomGenerator)
{
//  PRINT("spatialModelevalDiffCompartments_completeSpatialRandomness");

//  //open data info
//  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + ".tm" );

//  SpatialModelCompleteRandomness3D<float> triMeshSpatialModelDifferentCompartments;
//  triMeshSpatialModelDifferentCompartments.setRandomGenerator( randomGenerator );
//  triMeshSpatialModelDifferentCompartments.setTriMesh( nucleusTriMesh );
//  triMeshSpatialModelDifferentCompartments.initialize();

//  evalDiffCompartments(
//    nucleusTriMesh,
//    triMeshSpatialModelDifferentCompartments,
//    filename, parentDir,
//    function, constraints,
//    dataSet, randomGenerator );
}

void evalDiffCompartments_sizeConstrained(
  const string& filename, const string& parentDir,
  const string& function, const int& constraints,
  DataSet& dataSet,
  RandomGenerator& randomGenerator)
{
  PRINT("evalDiffCompartments_sizeConstrained");

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

  SpatialModelHardcoreDistance3DDifferentCompartments<float> triMeshSpatialModelDifferentCompartments;
  triMeshSpatialModelDifferentCompartments.setRandomGenerator( randomGenerator );
  triMeshSpatialModelDifferentCompartments.setTriMesh( nucleusTriMesh );
  triMeshSpatialModelDifferentCompartments.setHardcoreDistances( eqRadii, eqRadii );
  triMeshSpatialModelDifferentCompartments.initialize();
  triMeshSpatialModelDifferentCompartments.drawSample( numCCS, numCCS );

  Vertices<float> dist1, dist2;
  dist1 = triMeshSpatialModelDifferentCompartments.getVerticesDistribution1();
  dist2 = triMeshSpatialModelDifferentCompartments.getVerticesDistribution2();

  dist1.save( "/home/jarpon/Desktop/dist1.vx", true );
  dist2.save( "/home/jarpon/Desktop/dist2.vx", true );

//  evalDiffCompartments(
//    nucleusTriMesh,
//    triMeshSpatialModelDifferentCompartments,
//    filename, parentDir,
//    function, constraints,
//    dataSet, randomGenerator );
}

void evalDiffCompartments_distanceConstrained(
  const string& filename, const string& parentDir,
  const string& function, const int& constraints,
  DataSet& dataSet,
  RandomGenerator& randomGenerator)
{
  PRINT("spatialModelevalDiffCompartments_distanceConstrained");

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

  Vector<float> distancesToBorder( numCCS );
  int k = 0;
  for ( int j = lastPos - numCCS + 1 ; j < lastPos + 1; ++j, ++k )
    distancesToBorder[k] = ccsInfo.getValue<float>( "distanceToTheBorder", j );
  EVAL(distancesToBorder);


  SpatialModelBorderDistance3DDifferentCompartments<float> triMeshSpatialModelDifferentCompartments;
  triMeshSpatialModelDifferentCompartments.setRandomGenerator( randomGenerator );
  triMeshSpatialModelDifferentCompartments.setTriMesh( nucleusTriMesh );
  triMeshSpatialModelDifferentCompartments.setDistancesToBorder( distancesToBorder, distancesToBorder );
  triMeshSpatialModelDifferentCompartments.initialize();
  triMeshSpatialModelDifferentCompartments.drawSample( numCCS, numCCS );

  Vertices<float> dist1, dist2;
  dist1 = triMeshSpatialModelDifferentCompartments.getVerticesDistribution1();
  dist2 = triMeshSpatialModelDifferentCompartments.getVerticesDistribution2();

  dist1.save( "/home/jarpon/Desktop/dist1b.vx", true );
  dist2.save( "/home/jarpon/Desktop/dist2b.vx", true );

//  evalDiffCompartments(
//    nucleusTriMesh,
//    triMeshSpatialModelDifferentCompartments,
//    filename, parentDir,
//    function, constraints,
//    dataSet, randomGenerator );
}

void evalDiffCompartments_sizeAndDistanceConstrained(
  const string& filename, const string& parentDir,
  const string& function, const int& constraints,
  DataSet& dataSet,
  RandomGenerator& randomGenerator)
{
  PRINT("spatialModelevalDiffCompartments_sizeAndDistanceConstrained");

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

//  SpatialModelHardcoreBorderDistance3DDifferentCompartments<float>* triMeshSpatialModelDifferentCompartments;
//  triMeshSpatialModelDifferentCompartments->setRandomGenerator( randomGenerator );
//  triMeshSpatialModelDifferentCompartments->setTriMesh( nucleusTriMesh );
//  triMeshSpatialModelDifferentCompartments->setDistancesToBorder( distancesToBorder, distancesToBorder );
//  triMeshSpatialModelDifferentCompartments->setHardcoreDistances( eqRadii, eqRadii );
//  triMeshSpatialModelDifferentCompartments->initialize();
//  triMeshSpatialModelDifferentCompartments->drawSample( distancesToBorder.getSize(), distancesToBorder.getSize() );

//  evalDiffCompartments(
//    nucleusTriMesh, *triMeshSpatialModelDifferentCompartments, filename, parentDir,
//    function, constraints, dataSet, randomGenerator );
}

void evalDiffCompartments_maximalRepulsionConstrained(
  const string& filename, const string& parentDir,
  const string& function, const int& constraints,
  DataSet& dataSet,
  RandomGenerator& randomGenerator)
{
//  PRINT("spatialModelevalDiffCompartments_MaximalRepulsionConstrained");

//  //open data info
//  const string analysisDir = parentDir + "/analysis/";
//  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + ".tm" );
// // const DataSet datasetNucleus( analysisDir + filename + "_chromocenters.csv" );
////  const Vector<float> eqRadii = datasetNucleus.getValues<float>( "equivalentRadius_tm" );
// // const Vector<float> eqRadii = datasetNucleus.getValues<float>( "equivalentRadius_vm" );
////  EVAL(eqRadii);

//  const DataSet ccsInfo( analysisDir + "ccs.csv" );

//  Vector<string> tempFileNames;
//  tempFileNames = ccsInfo.getValues<string>( ccsInfo.variableNames()[0] );

//  int lastPos, numCCS = 0;

//  for ( int j = 0; j < tempFileNames.getSize(); ++j )
//    if ( tempFileNames[j] == filename )
//    {
//      lastPos = j;
//      ++ numCCS;
//    }

//  if ( numCCS == 0 )
//  {
//    EVAL("Nucleus not found");
//    return;
//  }

//  Vector<float> eqRadiiTemp( numCCS );
//  Vector<float> distancesToBorder( numCCS );
//  int k = 0;
//  for ( int j = lastPos - numCCS + 1 ; j < lastPos + 1; ++j, ++k )
//  {
//    eqRadiiTemp[k] = ccsInfo.getValue<float>( "equivalentRadius_ZprojCorrection", j );
////    eqRadiiTemp[k] = ccsInfo.getValue<float>( "equivalentRadius_PSFVolCorrection", j );
//    distancesToBorder[k] = ccsInfo.getValue<float>( "distanceToTheBorder", j );
//  }


//  for ( int i = 0; i < eqRadiiTemp.getSize(); ++i )
//    if ( eqRadiiTemp[i] > distancesToBorder[i] )
//      eqRadiiTemp[i] = distancesToBorder[i];

//  const Vector<float> eqRadii = eqRadiiTemp;
//  EVAL(eqRadii);
//  EVAL(distancesToBorder);


//  SpatialModelMaximalRepulsion3D<float> triMeshSpatialModelDifferentCompartments;
//  triMeshSpatialModelDifferentCompartments.setRandomGenerator( randomGenerator );
//  triMeshSpatialModelDifferentCompartments.setTriMesh( nucleusTriMesh );
//  triMeshSpatialModelDifferentCompartments.setNumMonteCarloCycles( 2000 );
//  triMeshSpatialModelDifferentCompartments.setHardcoreDistances( eqRadii );
//  triMeshSpatialModelDifferentCompartments.initialize();
//  triMeshSpatialModelDifferentCompartments.initializeBeta( eqRadii.getSize() );

//  evalDiffCompartments(
//    nucleusTriMesh, triMeshSpatialModelDifferentCompartments, filename, parentDir,
//    function, constraints, dataSet, randomGenerator );
}


//void evalDiffCompartments_MaxRepulsionWithDistanceToTheBorderConstrained(
//  const string& filename, const string& parentDir,
//  const string& function, const int& constraints,
//  DataSet& dataSet,
//  RandomGenerator& randomGenerator)
//{
//  PRINT("spatialModelevalDiffCompartments_MaxRepulsionWithDistanceToTheBorderConstrained");

//  const string analysisDir = parentDir + "/analysis/";

//  //new data
//  //const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + "_nucleus.tm" );
//  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + ".tm" );
//  //old data
//  //const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + "-nucleus.tm" );

//  //new data
//  const DataSet datasetNucleus( analysisDir + filename + "_chromocenters.csv" );
//  //old data
//  //const DataSet datasetNucleus( analysisDir + filename + ".csv" );
//  const Vector<float> distancesToBorder = datasetNucleus.getValues<float>( "distanceToTheBorder" );

//  //new data
//  const Vector<float> eqRadii = datasetNucleus.getValues<float>( "equivalentRadius_tm" );
//  //old data
//  //const Vector<float> eqRadii = datasetNucleus.getValues<float>( "chromocenterRadius" );
//  EVAL(eqRadii);

//  MaxRepulsionWithDistancesTriMeshSpatialModelDifferentCompartments<float> TriMeshSpatialModelDifferentCompartments;
//  TriMeshSpatialModelDifferentCompartments.setRandomGenerator( randomGenerator );
//  TriMeshSpatialModelDifferentCompartments.setTriMesh( nucleusTriMesh );
//  TriMeshSpatialModelDifferentCompartments.setDistanceToBorder( distancesToBorder );
//  TriMeshSpatialModelDifferentCompartments.setHardcoreDistance( eqRadii );
//  TriMeshSpatialModelDifferentCompartments.initialize();

//  evalDiffCompartments(
//    nucleusTriMesh, TriMeshSpatialModelDifferentCompartments, filename, parentDir,
//    function, constraints, dataSet, randomGenerator );
//}


/*! Chooses case depending on the constraints
****************************************************************/
void twoCompartmentsEvaluator(
  const string& filename, const string& parentDir,
  const string& function, const int& constraints,
  DataSet& dataSet,
  RandomGenerator& randomGenerator)
{
  ENTER("twoCompartmentsEvaluator(...)")
  switch( constraints )
  {
    case 0:
      evalDiffCompartments_completeSpatialRandomness(
        filename, parentDir,
        function, constraints,
        dataSet, randomGenerator );
      break;
    case 1:
      evalDiffCompartments_sizeConstrained(
        filename, parentDir,
        function, constraints,
        dataSet, randomGenerator );
      break;
    case 2:
      evalDiffCompartments_distanceConstrained(
        filename, parentDir,
        function, constraints,
        dataSet, randomGenerator );
      break;
    case 3:
      evalDiffCompartments_sizeAndDistanceConstrained(
        filename, parentDir,
        function, constraints,
        dataSet, randomGenerator );
      break;
    case 4:
      evalDiffCompartments_maximalRepulsionConstrained(
        filename, parentDir,
        function, constraints,
        dataSet, randomGenerator );
      break;
//    case 5:
//    evalDiffCompartments_MaxRepulsionWithDistanceToTheBorderConstrained(
//      filename, parentDir,
//      function, constraints,
//      dataSet, randomGenerator );
//    break;
  }
  LEAVE();
}
