#include <spatialmodelevaluator.h>
//#include "spatialmodelevaluator2.h"
#include <spatialdescriptorfunctionf.h>
#include <spatialdescriptorfunctiong.h>
#include <spatialdescriptorfunctionh.h>
#include <spatialdescriptorfunctionb.h>
#include <spatialdescriptorfunctionc.h>
#include "spatialdescriptorfunctionz.h"
#include <spatialmodelmaximalrepulsion3d.h>
#include "maxrepulsionwithdistances.h"
#include <spatialmodelcompleterandomness3d.h>
#include <spatialmodelborderdistance3d.h>
#include <spatialmodelhardcoreborderdistance3d.h>
#include <spatialmodelhardcoredistance3d.h>
#include <spatialmodel.h>
#include <trimesh.h>
#include <voxelmatrix.h>
#include <fileinfo.h>
//#include <regionanalysis.h>
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


void nucleoliEvaluator(
  const TriMesh<float>& nucleusTriMesh,
  TriMeshSpatialModel<float>& triMeshSpatialModel,
  const string& filename, const string& parentDir,
  const string& function, const int& constraints,
  DataSet& dataSet, RandomGenerator& randomGenerator)
{
  string classif = parentDir;
  const string analysisDir = parentDir + "/analysis/";

  classif = classif.substr(classif.find_last_of("/\\")+1,classif.length());

  //new data
//  const DataSet datasetNucleus( analysisDir + filename + "_chromocenters.csv" );
  //const DataSet datasetNucleus( analysisDir + filename + "_nucleoli.csv" );
  DataSet globalAnalysis( analysisDir + "nuclei_extended.csv" );
  const DataSet nucleoliInfo( analysisDir + "nucleoli.csv" );

  Vector<string> tempFileNames;
  tempFileNames = nucleoliInfo.getValues<string>( nucleoliInfo.variableNames()[0] );

  int lastPos, numNucleoli = 0;

  for ( int j = 0; j < tempFileNames.getSize(); ++j )
    if ( tempFileNames[j] == filename )
    {
      lastPos = j;
      ++ numNucleoli;
    }

  if ( numNucleoli == 0 )
  {
    EVAL("Nucleus not found");
    return;
  }


  Vertices<float> vertices( 3, numNucleoli, 0, 0 );
  int k = 0;
  for ( int j = lastPos - numNucleoli + 1 ; j < lastPos + 1; ++j, ++k )
  {
    vertices[k][0] = nucleoliInfo.getValue<float>( "centroidCoordX", j );
    vertices[k][1] = nucleoliInfo.getValue<float>( "centroidCoordY", j );
    vertices[k][2] = nucleoliInfo.getValue<float>( "centroidCoordZ", j );
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

  SpatialModelEvaluator<float,float> modelEvaluator;
  modelEvaluator.setModel( triMeshSpatialModel );
  modelEvaluator.setNumMonteCarloSamples( numPatterns ); //to check uniformity
  modelEvaluator.setPrecision( 0.05 );

  SpatialDescriptor<float>* spatialDescriptor;

  DataSet saveTest;

  //setting function parameters
  if ( function == "all" )
  {
    PRINT("all functions");

    SpatialModelCompleteRandomness3D<float> tempTriMeshSpatialModel;
    tempTriMeshSpatialModel.setRandomGenerator( randomGenerator );
    tempTriMeshSpatialModel.setTriMesh( nucleusTriMesh );
    tempTriMeshSpatialModel.initialize();
    Vertices<float> evaluationPositions = tempTriMeshSpatialModel.drawSample( 10000 );
    //evaluationPositions.save( parentDir + "/" + filename + "_Fpattern-" + ".vx", true); //to check uniformity of the F patterns

    SpatialDescriptorFunctionF<float>* spatialDescriptorFunctionF;
    spatialDescriptorFunctionF = new SpatialDescriptorFunctionF<float>();
    spatialDescriptorFunctionF->setEvaluationPositions( evaluationPositions );
    spatialDescriptor = spatialDescriptorFunctionF;
    modelEvaluator.addDescriptor( *spatialDescriptor );

    spatialDescriptor = new SpatialDescriptorFunctionG<float>();
    modelEvaluator.addDescriptor( *spatialDescriptor );

    spatialDescriptor = new SpatialDescriptorFunctionH<float>();
    modelEvaluator.addDescriptor( *spatialDescriptor );

    SpatialDescriptorFunctionB<float>* spatialDescriptorFunctionB;
    spatialDescriptorFunctionB = new SpatialDescriptorFunctionB<float>();
    spatialDescriptorFunctionB->setTriMesh( nucleusTriMesh );
    spatialDescriptor = spatialDescriptorFunctionB;
    modelEvaluator.addDescriptor( *spatialDescriptor );

    SpatialDescriptorFunctionC<float>* spatialDescriptorFunctionC;
    spatialDescriptorFunctionC = new SpatialDescriptorFunctionC<float>();
    spatialDescriptorFunctionC->setCenter( nucleusTriMesh.cog() );
    spatialDescriptor = spatialDescriptorFunctionC;
    modelEvaluator.addDescriptor( *spatialDescriptor );

    spatialDescriptor = new SpatialDescriptorFunctionZ<float>();
    modelEvaluator.addDescriptor( *spatialDescriptor );
  }
  else if ( function == "G" )
  {
    PRINT("G");
    spatialDescriptor = new SpatialDescriptorFunctionG<float>();
  }
  else if ( function == "H" )
  {
    PRINT("H");
    spatialDescriptor = new SpatialDescriptorFunctionH<float>();
  }
  else if ( function == "B" )
  {
    PRINT("B");
    SpatialDescriptorFunctionB<float>* spatialDescriptorFunctionB;
    spatialDescriptorFunctionB = new SpatialDescriptorFunctionB<float>();
    spatialDescriptorFunctionB->setTriMesh( nucleusTriMesh );
    spatialDescriptor = spatialDescriptorFunctionB;
  }
  else if ( function == "C" )
  {
    PRINT("C");
    SpatialDescriptorFunctionC<float>* spatialDescriptorFunctionC;
    spatialDescriptorFunctionC = new SpatialDescriptorFunctionC<float>();
    spatialDescriptorFunctionC->setCenter( nucleusTriMesh.cog() );
    spatialDescriptor = spatialDescriptorFunctionC;
  }
  else if ( function == "Z" )
  {
    PRINT("Z");
    spatialDescriptor = new SpatialDescriptorFunctionZ<float>();
  }
  else //if ( function == "F" )
  {
    PRINT("F");
    SpatialModelCompleteRandomness3D<float> tempTriMeshSpatialModel;
    tempTriMeshSpatialModel.setRandomGenerator( randomGenerator );
    tempTriMeshSpatialModel.setTriMesh( nucleusTriMesh );
    tempTriMeshSpatialModel.initialize();
    Vertices<float> evaluationPositions = tempTriMeshSpatialModel.drawSample( 10000 );
    //evaluationPositions.save( parentDir + "/spatial_models/" + filename + "_Fpattern-" + ".vx", true);
    SpatialDescriptorFunctionF<float>* spatialDescriptorFunctionF;
    spatialDescriptorFunctionF = new SpatialDescriptorFunctionF<float>();
    spatialDescriptorFunctionF->setEvaluationPositions( evaluationPositions );
    spatialDescriptor = spatialDescriptorFunctionF;
  }


  //processing data
  if ( function != "all" )
  {
    string newVariable;
    switch ( constraints )
    {
      case 0:
        newVariable = ( "SpatialModelCompleteRandomness3D_" + function + "-SDI" );
        break;
      case 1:
        newVariable = ( "SpatialModelHardcoreDistance3D_" + function + "-SDI" );
        break;
      case 2:
        newVariable = ( "spatialModelBorderDistance3D_" + function + "-SDI" );
        break;
      case 3:
        newVariable = ( "SpatialModelBorderHardcoreDistance3D_" + function + "-SDI" );
        break;
      case 4:
        newVariable = ( "SpatialModelMaximalRepulsion3D_" + function + "-SDI" );
        break;
    }

    //globalAnalysis.setValues<float>( newVariable, -1 );
    float sdi;
    int row;

    modelEvaluator.setDescriptor( *spatialDescriptor );

    ostringstream iss; //we suppose as much 99 labels
    iss << constraints;
    DataSet saveTest;
    try {
      sdi = modelEvaluator.eval( vertices, &saveTest );
      EVAL(sdi);

  //    Vector<float> output = modelEvaluator.evalSDIandMaxDiff( vertices, &saveTest );
  //    EVAL( output[0] );
  //    EVAL( output[1] );


      saveTest.save( analysisDir + iss.str() + "/" + function + "/" + filename + ".csv", true );
    //  saveTest.save( analysisDir + iss.str() + "/" + function + "/" + filename + "_random.csv", true );
      row = dataSet.size()[0];

      dataSet.setValue( "nucleus", row, filename );
      dataSet.setValue( "class", row, classif );//classification: mutant, tissue, etc.
      dataSet.setValue( "descriptor", row, function );//spatial descriptor
      dataSet.setValue( "index", row, sdi );
      //dataSet.setValue( "index", row, output[0] );
      //dataSet.setValue( "signedMaxDiff", row, output[1] );

      globalAnalysis.setValue( newVariable, numCurrentNucleus, sdi );
    }
    catch( Exception exception ) {
      EVAL( exception.getWhat() );
      EVAL( exception.getWhere() );

      const int row = dataSet.size()[0];

      dataSet.setValue( "nucleus", row, filename );
      dataSet.setValue( "class", row, classif );//classification: mutant, tissue, etc.
      dataSet.setValue( "descriptor", row, function );//spatial descriptor
      dataSet.setValue( "index", row, sqrt(-1) );


    }

  }
  else
  {
//    Vertices<float> vertices ( 3, numPoints, 0, 0 );
//    for ( int i = 0; i < numPoints; ++i )
//    {
//      vertices[i][0] = datasetNucleus.getValue<float>( "centroidCoordX", i );
//      vertices[i][1] = datasetNucleus.getValue<float>( "centroidCoordY", i );
//      vertices[i][2] = datasetNucleus.getValue<float>( "centroidCoordZ", i );
//      EVAL(vertices[i]);
//    }

    const int row = dataSet.numRows();
    vector<float> sdis;
    vector<int> ranks;
    DataSet saveTest;

//    vector<float> maxDiff;

    try {
//    modelEvaluator.evalSDIandMaxDiff( vertices, pValues, ranks, maxDiff);

    modelEvaluator.eval( vertices, sdis, ranks, &saveTest );

    EVAL( sdis[0] );
    EVAL( sdis[1] );
    EVAL( sdis[2] );
    EVAL( sdis[3] );
    EVAL( sdis[4] );
    EVAL( sdis[5] );

    string newVariable;
    switch ( constraints )
    {
      case 0: newVariable = "SpatialModelCompleteRandomness3D_";
      case 1: newVariable = "SpatialModelHardcoreDistance3D_";
      case 2: newVariable = "spatialModelBorderDistance3D_";
      case 3: newVariable = "SpatialModelBorderHardcoreDistance3D_";
      case 4: newVariable = "SpatialModelMaximalRepulsion3D_";
    }

    //unifying datasets
    for ( int jj = 0; jj < sdis.size(); ++jj )
    {
      if ( ( sdis[jj] < 0 || sdis[jj] > 1 )  && ( jj = 0 ) )
        globalAnalysis.setValue( newVariable + "F-SDI", numCurrentNucleus, sqrt(-1) );
      else
        globalAnalysis.setValue( newVariable + "F-SDI", numCurrentNucleus, sdis[jj] );
      if ( ( sdis[jj] < 0 || sdis[jj] > 1 )  && ( jj = 1 ) )
        globalAnalysis.setValue( newVariable + "G-SDI", numCurrentNucleus, sqrt(-1) );
      else
        globalAnalysis.setValue( newVariable + "G-SDI", numCurrentNucleus, sdis[jj] );
      if ( ( sdis[jj] < 0 || sdis[jj] > 1 )  && ( jj = 2 ) )
        globalAnalysis.setValue( newVariable + "H-SDI", numCurrentNucleus, sqrt(-1) );
      else
        globalAnalysis.setValue( newVariable + "H-SDI", numCurrentNucleus, sdis[jj] );
      if ( ( sdis[jj] < 0 || sdis[jj] > 1 )  && ( jj = 3 ) )
        globalAnalysis.setValue( newVariable + "B-SDI", numCurrentNucleus, sqrt(-1) );
      else
        globalAnalysis.setValue( newVariable + "B-SDI", numCurrentNucleus, sdis[jj] );
      if ( ( sdis[jj] < 0 || sdis[jj] > 1 )  && ( jj = 4 ) )
        globalAnalysis.setValue( newVariable + "C-SDI", numCurrentNucleus, sqrt(-1) );
      else
        globalAnalysis.setValue( newVariable + "C-SDI", numCurrentNucleus, sdis[jj] );
      if ( ( sdis[jj] < 0 || sdis[jj] > 1 )  && ( jj = 5 ) )
        globalAnalysis.setValue( newVariable + "Z-SDI", numCurrentNucleus, sqrt(-1) );
      else
        globalAnalysis.setValue( newVariable + "Z-SDI", numCurrentNucleus, sdis[jj] );
    }

    dataSet.setValue( "nucleus", row, filename );
    dataSet.setValue( "class", row, classif );//classification: mutant, tissue, etc.
    dataSet.setValue( "F-SDI", row, sdis[0] );
//    dataSet.setValue( "F-maxDiff", row, maxDiff[0] );
    dataSet.setValue( "G-SDI", row, sdis[1] );
//    dataSet.setValue( "G-maxDiff", row, maxDiff[1] );
    dataSet.setValue( "H-SDI", row, sdis[2] );
//    dataSet.setValue( "H-maxDiff", row, maxDiff[2] );
    dataSet.setValue( "B-SDI", row, sdis[3] );
//    dataSet.setValue( "B-maxDiff", row, maxDiff[3] );
    dataSet.setValue( "C-SDI", row, sdis[4] );
//    dataSet.setValue( "C-maxDiff", row, maxDiff[4] );
    dataSet.setValue( "Z-SDI", row, sdis[4] );
//    dataSet.setValue( "Z-maxDiff", row, maxDiff[5] );

    ostringstream iss; //we have 4 constraints
    iss << constraints;
    saveTest.save( analysisDir + iss.str() + "/" + filename + ".csv", true );

    }
    catch( Exception exception ) {
      EVAL( exception.getWhat() );
      EVAL( exception.getWhere() );
      const int row = dataSet.size()[0];

      dataSet.setValue( "nucleus", row, filename );
      dataSet.setValue( "class", row, classif );//classification: mutant, tissue, etc.
      dataSet.setValue( "descriptor", row, function );//spatial descriptor
      dataSet.setValue( "index", row, sqrt(-1) );
      //dataSet.setValue( "index", row, output[0] );
      //dataSet.setValue( "signedMaxDiff", row, output[1] );

    }

  }

  globalAnalysis.save( analysisDir + "/" + "nuclei_complete2.csv", true );


}


/*! Preparing constraints for the model
****************************************************************/
//void nucleoliEvaluator_completeSpatialRandomness(
//  const string& filename, const string& parentDir,
//  const string& function, const int& constraints,
//  DataSet& dataSet,
//  RandomGenerator& randomGenerator)
//{
//  PRINT("SpatialModelEvaluator_completeSpatialRandomness");

//  //const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + "_nucleus.tm" );
//  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + ".tm" );
//  TriMeshSpatialModel<float> triMeshSpatialModel;
//  triMeshSpatialModel.setRandomGenerator( randomGenerator );
//  triMeshSpatialModel.setTriMesh( nucleusTriMesh );
//  triMeshSpatialModel.initialize();

//  nucleoliEvaluator(
//    nucleusTriMesh,
//    triMeshSpatialModel,
//    filename, parentDir,
//    function, constraints,
//    dataSet, randomGenerator );
//}

void nucleoliEvaluator_sizeConstrained(
  const string& filename, const string& parentDir,
  const string& function, const int& constraints,
  DataSet& dataSet,
  RandomGenerator& randomGenerator)
{
  PRINT("SpatialModelEvaluator_sizeConstrained");

  //open data info
  const string analysisDir = parentDir + "/analysis/";
  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + ".tm" );

  const DataSet nucleoliInfo( analysisDir + "nucleoli.csv" );

  Vector<string> tempFileNames;
  tempFileNames = nucleoliInfo.getValues<string>( nucleoliInfo.variableNames()[0] );

  int lastPos, numNucleoli = 0;

  for ( int j = 0; j < tempFileNames.getSize(); ++j )
    if ( tempFileNames[j] == filename )
    {
      lastPos = j;
      ++ numNucleoli;
    }

  if ( numNucleoli == 0 )
  {
    EVAL("Nucleus not found");
    return;
  }

  Vector<float> eqRadiiTemp( numNucleoli );
  Vector<float> distancesToBorder( numNucleoli );
  int k = 0;
  for ( int j = lastPos - numNucleoli + 1 ; j < lastPos + 1; ++j, ++k )
  {
    eqRadiiTemp[k] = nucleoliInfo.getValue<float>( "equivalentRadius_vm", j );
    //eqRadiiTemp[k] = nucleoliInfo.getValue<float>( "equivalentRadius_ZprojCorrection", j );
//    eqRadiiTemp[k] = nucleoliInfo.getValue<float>( "equivalentRadius_PSFVolCorrection", j );
//    distancesToBorder[k] = nucleoliInfo.getValue<float>( "distanceToTheBorder", j );
  }


  for ( int i = 0; i < eqRadiiTemp.getSize(); ++i )
    if ( eqRadiiTemp[i] > distancesToBorder[i] )
      eqRadiiTemp[i] = distancesToBorder[i];

  const Vector<float> eqRadii = eqRadiiTemp;
  EVAL(eqRadii);
//  EVAL(distancesToBorder);

  SpatialModelBorderDistance3D<float> triMeshSpatialModel;
  triMeshSpatialModel.setRandomGenerator( randomGenerator );
  triMeshSpatialModel.setTriMesh( nucleusTriMesh );
  triMeshSpatialModel.setDistancesToBorder( distancesToBorder );
  triMeshSpatialModel.initialize();

  nucleoliEvaluator(
        nucleusTriMesh, triMeshSpatialModel, filename, parentDir,
        function, constraints, dataSet, randomGenerator );
}

void nucleoliEvaluator_distanceConstrained(
  const string& filename, const string& parentDir,
  const string& function, const int& constraints,
  DataSet& dataSet,
  RandomGenerator& randomGenerator)
{
  PRINT("SpatialModelEvaluator_distanceConstrained");

  //open data info
  const string analysisDir = parentDir + "/analysis/";
  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + ".tm" );

  const DataSet nucleoliInfo( analysisDir + "nucleoli.csv" );

  Vector<string> tempFileNames;
  tempFileNames = nucleoliInfo.getValues<string>( nucleoliInfo.variableNames()[0] );

  int lastPos, numNucleoli = 0;

  for ( int j = 0; j < tempFileNames.getSize(); ++j )
    if ( tempFileNames[j] == filename )
    {
      lastPos = j;
      ++ numNucleoli;
    }

  if ( numNucleoli == 0 )
  {
    EVAL("Nucleus not found");
    return;
  }

//  Vector<float> eqRadiiTemp( numNucleoli );
  Vector<float> distancesToBorder( numNucleoli );
  int k = 0;
  for ( int j = lastPos - numNucleoli + 1 ; j < lastPos + 1; ++j, ++k )
  {
//    eqRadiiTemp[k] = nucleoliInfo.getValue<float>( "equivalentRadius_vm", j );
    //eqRadiiTemp[k] = nucleoliInfo.getValue<float>( "equivalentRadius_ZprojCorrection", j );
//    eqRadiiTemp[k] = nucleoliInfo.getValue<float>( "equivalentRadius_PSFVolCorrection", j );
    distancesToBorder[k] = nucleoliInfo.getValue<float>( "distanceToTheBorder", j );
  }

//  const Vector<float> eqRadii = eqRadiiTemp;
//  EVAL(eqRadii);
  EVAL(distancesToBorder);

  SpatialModelBorderDistance3D<float> triMeshSpatialModel;
  triMeshSpatialModel.setRandomGenerator( randomGenerator );
  triMeshSpatialModel.setTriMesh( nucleusTriMesh );
  triMeshSpatialModel.setDistancesToBorder( distancesToBorder );
  triMeshSpatialModel.initialize();

  nucleoliEvaluator(
        nucleusTriMesh, triMeshSpatialModel, filename, parentDir,
        function, constraints, dataSet, randomGenerator );
}

void nucleoliEvaluator_sizeAndDistanceConstrained(
  const string& filename, const string& parentDir,
  const string& function, const int& constraints,
  DataSet& dataSet,
  RandomGenerator& randomGenerator)
{
  PRINT("SpatialModelEvaluator_sizeAndDistanceConstrained");

  //open data info
  const string analysisDir = parentDir + "/analysis/";
  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + ".tm" );

  const DataSet nucleoliInfo( analysisDir + "nucleoli.csv" );

  Vector<string> tempFileNames;
  tempFileNames = nucleoliInfo.getValues<string>( nucleoliInfo.variableNames()[0] );

  int lastPos, numNucleoli = 0;

  for ( int j = 0; j < tempFileNames.getSize(); ++j )
    if ( tempFileNames[j] == filename )
    {
      lastPos = j;
      ++ numNucleoli;
    }

  if ( numNucleoli == 0 )
  {
    EVAL("Nucleus not found");
    return;
  }

  Vector<float> eqRadiiTemp( numNucleoli );
  Vector<float> distancesToBorder( numNucleoli );
  int k = 0;
  for ( int j = lastPos - numNucleoli + 1 ; j < lastPos + 1; ++j, ++k )
  {
    eqRadiiTemp[k] = nucleoliInfo.getValue<float>( "equivalentRadius_vm", j );
    //eqRadiiTemp[k] = nucleoliInfo.getValue<float>( "equivalentRadius_ZprojCorrection", j );
//    eqRadiiTemp[k] = nucleoliInfo.getValue<float>( "equivalentRadius_PSFVolCorrection", j );
    distancesToBorder[k] = nucleoliInfo.getValue<float>( "distanceToTheBorder", j );
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

  nucleoliEvaluator(
        nucleusTriMesh, triMeshSpatialModel, filename, parentDir,
        function, constraints, dataSet, randomGenerator );
}

void nucleoliEvaluator_MaximalRepulsionConstrained(
  const string& filename, const string& parentDir,
  const string& function, const int& constraints,
  DataSet& dataSet,
  RandomGenerator& randomGenerator)
{
  PRINT("SpatialModelEvaluator_MaximalRepulsionConstrained");

  //open data info
  const string analysisDir = parentDir + "/analysis/";
  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + ".tm" );

  const DataSet nucleoliInfo( analysisDir + "nucleoli.csv" );

  Vector<string> tempFileNames;
  tempFileNames = nucleoliInfo.getValues<string>( nucleoliInfo.variableNames()[0] );

  int lastPos, numNucleoli = 0;

  for ( int j = 0; j < tempFileNames.getSize(); ++j )
    if ( tempFileNames[j] == filename )
    {
      lastPos = j;
      ++ numNucleoli;
    }

  if ( numNucleoli == 0 )
  {
    EVAL("Nucleus not found");
    return;
  }

  Vector<float> eqRadiiTemp( numNucleoli );
  Vector<float> distancesToBorder( numNucleoli );
  int k = 0;
  for ( int j = lastPos - numNucleoli + 1 ; j < lastPos + 1; ++j, ++k )
  {
    eqRadiiTemp[k] = nucleoliInfo.getValue<float>( "equivalentRadius_vm", j );
    //eqRadiiTemp[k] = nucleoliInfo.getValue<float>( "equivalentRadius_ZprojCorrection", j );
//    eqRadiiTemp[k] = nucleoliInfo.getValue<float>( "equivalentRadius_PSFVolCorrection", j );
    distancesToBorder[k] = nucleoliInfo.getValue<float>( "distanceToTheBorder", j );
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
  triMeshSpatialModel.initializeBeta( numNucleoli );

  nucleoliEvaluator(
        nucleusTriMesh, triMeshSpatialModel, filename, parentDir,
        function, constraints, dataSet, randomGenerator );
}


//void nucleoliEvaluator_MaxRepulsionWithDistanceToTheBorderConstrained(
//  const string& filename, const string& parentDir,
//  const string& function, const int& constraints,
//  DataSet& dataSet,
//  RandomGenerator& randomGenerator)
//{
//  PRINT("SpatialModelEvaluator_MaxRepulsionWithDistanceToTheBorderConstrained");

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

//  MaxRepulsionWithDistancesTriMeshSpatialModel<float> triMeshSpatialModel;
//  triMeshSpatialModel.setRandomGenerator( randomGenerator );
//  triMeshSpatialModel.setTriMesh( nucleusTriMesh );
//  triMeshSpatialModel.setDistanceToBorder( distancesToBorder );
//  triMeshSpatialModel.setHardcoreDistance( eqRadii );
//  triMeshSpatialModel.initialize();

//  nucleoliEvaluator(
//    nucleusTriMesh, triMeshSpatialModel, filename, parentDir,
//    function, constraints, dataSet, randomGenerator );
//}


/*! Chooses case depending on the constraints
****************************************************************/
void nucleoliEvaluator(
  const string& filename, const string& parentDir,
  const string& function, const int& constraints,
  DataSet& dataSet,
  RandomGenerator& randomGenerator)
{
  switch( constraints )
  {
//    case 0:
//      nucleoliEvaluator_completeSpatialRandomness(
//        filename, parentDir,
//        function, constraints,
//        dataSet, randomGenerator );
//      break;
    case 1:
      nucleoliEvaluator_sizeConstrained(
        filename, parentDir,
        function, constraints,
        dataSet, randomGenerator );
      break;
    case 2:
      nucleoliEvaluator_distanceConstrained(
        filename, parentDir,
        function, constraints,
        dataSet, randomGenerator );
      break;
    case 3:
      nucleoliEvaluator_sizeAndDistanceConstrained(
        filename, parentDir,
        function, constraints,
        dataSet, randomGenerator );
      break;
    case 4:
      nucleoliEvaluator_MaximalRepulsionConstrained(
        filename, parentDir,
        function, constraints,
        dataSet, randomGenerator );
      break;
//    case 5:
//    nucleoliEvaluator_MaxRepulsionWithDistanceToTheBorderConstrained(
//      filename, parentDir,
//      function, constraints,
//      dataSet, randomGenerator );
//    break;
  }
}

