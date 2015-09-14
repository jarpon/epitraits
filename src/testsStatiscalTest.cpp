#include <ask4.h>
#include <dataset.h>
#include <pixelmatrix.h>
#include <spatialdescriptorfunctionf.h>
#include <spatialdescriptorfunctiong.h>
#include <spatialdescriptorfunctionh.h>
#include "spatialdescriptorborder.h"
#include "spatialdescriptorcentroid.h"
#include "spatialmodelboundaryinteraction.h"
#include <spatialmodelcompleterandomness.h>
#include "spatialmodelevaluator2.h"
#include "spatialdescriptorborder2D.h"


#include <cmath>

#define CoordType float
#define PixelType float

#define TRACE
#include <trace.h>

extern void testsStatisticalTests();

// 100 iterations
//
// 1 sample of 100 patterns with 50 points using spatial model 1
// 1 sample of 100 patterns with p {10,20,30,40,50,60,70,80,90,100} points using spatial model 1

// compare samples against CSR model

// test sdi indexes using KS test



//vector<int> initNumMonteCarloSamples()
//{
//  vector<int> numMonteCarloSamples;
//  numMonteCarloSamples.push_back( 9 );
//  numMonteCarloSamples.push_back( 19 );
//  numMonteCarloSamples.push_back( 39 );
//  return numMonteCarloSamples;
//}

//vector<int> initSampleSizes()
//{
//  vector<int> sampleSizes;
//  sampleSizes.push_back( 10 );
//  for (int s = 50; s <= 1000; s += 50)
//    sampleSizes.push_back( s );
//  return sampleSizes;
//}

vector<int> initNumObjectsS2()
{
  vector<int> numObjectsS2;
  numObjectsS2.push_back( 10 );
  for (int s = 10; s <= 100; s += 10)
    numObjectsS2.push_back( s );
  return numObjectsS2;
}

Curve<CoordType> initSquareBoundary()
{
  Curve<CoordType> boundary( 2, 4 );
  boundary[0][X] = 0.0;
  boundary[0][Y] = 0.0;
  boundary[1][X] = 100.0;
  boundary[1][Y] = 0.0;
  boundary[2][X] = 100.0;
  boundary[2][Y] = 100.0;
  boundary[3][X] = 0.0;
  boundary[3][Y] = 100.0;
  return boundary;
}

PixelMatrix<PixelType> initSquareMask()
{
  PixelMatrix<PixelType> pixelMatrix( 200, 200 );
  pixelMatrix.fill( 0 );
  pixelMatrix.fillBorder( 50, 255 );
  return pixelMatrix;
}

Curve<CoordType> initSinusoidalBoundary()
{
  const int n = 100;
  Curve<CoordType> boundary( 2, n );
  CoordType r;

  for (int v = 0; v < n; ++v)
  {
    r = 8.0 + 2.0*cos( static_cast<float>(16.0*M_PI*v)/n );
    boundary[v][X] = r * cos( static_cast<float>(2.0*M_PI*v)/n );
    boundary[v][Y] = r * sin( static_cast<float>(2.0*M_PI*v)/n );
  }

  // on redimensionne le contour afin d'avoir la même aire
  // que le carré de taille 100 par 100
  boundary.scale( 100.0*sqrt(2.0/boundary.twiceArea()) );

  return boundary;
}

//void testPatternGeneration(
//  SpatialModel<CoordType,PixelType>& spatialModel,
//  const int numObjects,
//  const string& id)
//{
//  Vertices<CoordType> pattern;
//  pattern = spatialModel.drawSample( numObjects );
//  pattern.save( "pattern-"+id, true );

//  EVAL( pattern.squareNearestNeighborDistances() );
//  EVAL( pattern.squareNearestNeighborDistances().min() );

//  SpatialDescriptorFunctionF<CoordType> functionF;
//  SpatialDescriptorFunctionG<CoordType> functionG;
//  SpatialDescriptorFunctionH<CoordType> functionH;
//  SpatialDescriptorDistanceToBorder<CoordType> functionB;
//  SpatialDescriptorDistanceToCentroid<CoordType> functionC;
//  SpatialModelEvaluator<CoordType,PixelType> spatialModelEvaluator;
//  spatialModelEvaluator.setModel( spatialModel );
//  spatialModelEvaluator.setNumMonteCarloSamples( 99 );
//  spatialModelEvaluator.addDescriptor( functionF );
//  spatialModelEvaluator.addDescriptor( functionG );
//  spatialModelEvaluator.addDescriptor( functionH );
//  spatialModelEvaluator.addDescriptor( functionB );
//  spatialModelEvaluator.addDescriptor( functionC );
//  spatialModelEvaluator.setPrecision( 1.0 );
//  functionF.setEvaluationPositions( spatialModel.drawSample(1000) );

//  vector<float> pValues;
//  vector<int> ranks;
//  spatialModelEvaluator.eval( pattern, pValues, ranks );
//  EVAL( pValues[0] );
//  EVAL( pValues[1] );
//  EVAL( pValues[2] );
//  EVAL( pValues[3] );
//  EVAL( pValues[4] );

//}


//void testPatternGeneration()
//{
//  SpatialModelCompleteRandomness<CoordType,PixelType> spatialModel;
//  const Curve<CoordType> boundary = initSinusoidalBoundary();
//  RandomGenerator randomGenerator;
//  //const int numObjects = ask4int( "Number of objects ");
//  const int numObjectsS1 = 50;

//  spatialModel1.setRandomGenerator( randomGenerator );
//  spatialModel1.setBoundary( boundary );
//  spatialModel1.initialize();
//  testPatternGeneration( spatialModel1, numObjectsS1, "sinusoidal" );
//  boundary.save( "sinusoidalBoundary", true );
//}

void testPattern(
  const Vertices<CoordType>& pattern1,
  const Vertices<CoordType>& pattern2,
//  SpatialModel<CoordType,PixelType>& spatialModel1,
//  SpatialModel<CoordType,PixelType>& spatialModel2,
  SpatialModel<CoordType,PixelType>& csrModel,
  const Curve<CoordType> boundary,
  DataSet& gDataSet)
{
  ENTER( "void testPattern(...)" );
//  SpatialDescriptorFunctionF<CoordType> functionF;
//  SpatialDescriptorFunctionG<CoordType> functionG;
//  SpatialDescriptorFunctionH<CoordType> functionH;
//  SpatialDescriptorDistanceToBorder2D<CoordType>* functionB;
//  functionB = new SpatialDescriptorDistanceToBorder2D<CoordType>();
  SpatialDescriptorFunctionF<CoordType> functionF;
  SpatialDescriptorFunctionG<CoordType> functionG;
  SpatialDescriptorFunctionH<CoordType> functionH;
//  SpatialDescriptorDistanceToBorder2D<float>* functionB;
//  functionB->setBoundary( boundary );
//  SpatialModelEvaluator<CoordType,PixelType> spatialModelEvaluator;
  SpatialModelEvaluator<CoordType,PixelType> spatialModelEvaluator;
  spatialModelEvaluator.setModel( csrModel );
  spatialModelEvaluator.setNumMonteCarloSamples( 99 );
  spatialModelEvaluator.addDescriptor( functionF );
  spatialModelEvaluator.addDescriptor( functionG );
  spatialModelEvaluator.addDescriptor( functionH );
  //spatialModelEvaluator.addDescriptor( functionB );
  spatialModelEvaluator.setPrecision( 1.0 );
  functionF.setEvaluationPositions( csrModel.drawSample(1000) );

	EVAL("spatialmodel 1");
  const int r = gDataSet.numRows();
  vector<float> pValues;
  vector<int> ranks;
  vector<float> maxDiff;

  //spatialModelEvaluator.eval( pattern1, pValues, ranks );
  spatialModelEvaluator.evalSDIandMaxDiff( pattern1, pValues, ranks, maxDiff);
  gDataSet.setValue( "F1-SDI", r, pValues[0] );
  gDataSet.setValue( "F1-maxDiff", r, maxDiff[0] );
  gDataSet.setValue( "G1-SDI", r, pValues[1] );
  gDataSet.setValue( "G1-maxDiff", r, maxDiff[1] );
  gDataSet.setValue( "H1-SDI", r, pValues[2] );
  gDataSet.setValue( "H1-maxDiff", r, maxDiff[2] );
  //gDataSet.setValue( "B1-SDI", r, pValues[3] );



  	EVAL("ok");
 // SpatialDescriptorFunctionF<CoordType> functionF2;
 // SpatialDescriptorFunctionG<CoordType> functionG2;
 // SpatialDescriptorFunctionH<CoordType> functionH2;
 // SpatialDescriptorDistanceToBorder<CoordType> functionB2;
 // SpatialDescriptorDistanceToCentroid<CoordType> functionC2;
  SpatialModelEvaluator<CoordType,PixelType> spatialModelEvaluator2;
  spatialModelEvaluator2.setModel( csrModel );
  spatialModelEvaluator2.setNumMonteCarloSamples( 99 );
  spatialModelEvaluator2.addDescriptor( functionF );
  spatialModelEvaluator2.addDescriptor( functionG );
  spatialModelEvaluator2.addDescriptor( functionH );
  //spatialModelEvaluator2.addDescriptor( functionB );
  spatialModelEvaluator2.setPrecision( 1.0 );
  functionF.setEvaluationPositions( csrModel.drawSample(1000) );

	EVAL("spatialmodel 2");

  vector<float> pValues2;
  vector<int> ranks2;
  vector<float> maxDiff2;
  //spatialModelEvaluator2.eval( pattern2, pValues2, ranks2 );
  spatialModelEvaluator2.evalSDIandMaxDiff( pattern2, pValues2, ranks2, maxDiff2 );
  gDataSet.setValue( "F2-SDI", r, pValues2[0] );
  gDataSet.setValue( "F2-maxDiff", r, maxDiff2[0] );
  gDataSet.setValue( "G2-SDI", r, pValues2[1] );
  gDataSet.setValue( "G2-maxDiff", r, maxDiff2[1] );
  gDataSet.setValue( "H2-SDI", r, pValues2[2] );
  gDataSet.setValue( "H2-maxDiff", r, maxDiff2[2] );


  	EVAL("ok");

  LEAVE();
}

void testModel(
  SpatialModel<CoordType,PixelType>& spatialModel1,
  SpatialModel<CoordType,PixelType>& spatialModel2,
  SpatialModel<CoordType,PixelType>& csrModel,
  const int numObjectsS1,
  const int numObjectsS2,
  const Curve<CoordType> boundary,
  DataSet& gDataSet)
{
  ENTER( "void testModel(...)" );
//  EVAL( sampleSize );
  const int sampleSize = 100;
  for (int i = 0; i < sampleSize; ++i)
  {
    const Vertices<CoordType> pattern1 = spatialModel1.drawSample( numObjectsS1 );
	//EVAL(pattern1.size());
	EVAL("1 ok");
    const Vertices<CoordType> pattern2 = spatialModel2.drawSample( numObjectsS2 );
	//EVAL(pattern2.size());
	EVAL("2 ok");
	EVAL(numObjectsS2);
    pattern1.save( "current-pattern1", true );
    pattern2.save( "current-pattern2", true );
    //testPattern( pattern1, pattern2, spatialModel1, spatialModel2, csrModel, boundary, gDataSet );
    testPattern( pattern1, pattern2, csrModel, boundary, gDataSet );
  }
  LEAVE();
}

void testModel(
  SpatialModel<CoordType,PixelType>& spatialModel1,
  SpatialModel<CoordType,PixelType>& spatialModel2,
  SpatialModel<CoordType,PixelType>& csrModel,
  const int numObjectsS1,
  const int numObjectsS2,
  const Curve<CoordType> boundary,
  const int numBatches,
  DataSet& gDataSet)
{
  ENTER( "void testModel(...)" );
  for (int b = 1; b <= numBatches; ++b)
  {
    EVAL( b );
    DataSet dataSet;
    testModel( spatialModel1, spatialModel2, csrModel, numObjectsS1, numObjectsS2, boundary, dataSet );
    dataSet.setValues( "batch", b );
    gDataSet.append( dataSet );
  }
  LEAVE();
}

void testModel(
  SpatialModel<CoordType,PixelType>& spatialModel1,
  SpatialModel<CoordType,PixelType>& spatialModel2,
  SpatialModel<CoordType,PixelType>& csrModel,
  const int numObjectsS1,
  const vector<int>& numObjectsS2,
  const Curve<CoordType> boundary,
  const int numBatches,
  DataSet& gDataSet)
{
  for (size_t s = 0; s < numObjectsS2.size(); ++s)
    {
      DataSet dataSet;
      testModel( spatialModel1, spatialModel2, csrModel, numObjectsS1, numObjectsS2[s], boundary, numBatches, dataSet );
      dataSet.setValues( "numObjectsS2", numObjectsS2[s] );
      gDataSet.append( dataSet );
      gDataSet.save( "current-dataset.data", true );
    }
}

//void testCSROverBoundary(const Curve<CoordType>& boundary, const string& id)
//{
//  ENTER( "void testCSROverBoundary(...)" );
//  RandomGenerator randomGenerator;
//  SpatialModelCompleteRandomness<CoordType,PixelType> spatialModel;
//  spatialModel.setBoundary( boundary );
//  spatialModel.setRandomGenerator( randomGenerator );
//  spatialModel.initialize();

//  DataSet dataSet;
//  const int numObjectsS1 = 50;
//  const int numBatches = 100;
//  testModel( spatialModel, spatialModel, numObjectsS1, initSampleSizes(), initNumMonteCarloSamples(), numBatches, dataSet );
//  dataSet.save( id+".data", true );
//  LEAVE();
//}

//void testCSROverSquareBoundary()
//{
//  ENTER( "void testCSROverSquareBoundary()" );
//  const Curve<CoordType> boundary = initSquareBoundary();
//  boundary.save( "squareBoundary", true );
//  testCSROverBoundary( boundary, "testCSROverSquareBoundary" );
//  LEAVE();
//}

//void testCSROverSinusoidalBoundary()
//{
//  ENTER( "void testCSROverSinusoidalBoundary()" );
//  const Curve<CoordType> boundary = initSinusoidalBoundary();
//  boundary.save( "sinusoidalBoundary", true );
//  testCSROverBoundary( boundary, "testCSROverSinusoidalBoundary" );
//  LEAVE();
//}

//void testCSROverSquareMask()
//{
//  ENTER( "void testCSROverSquare()" );
//  const PixelMatrix<PixelType> pixelMatrix = initSquareMask();
//  RandomGenerator randomGenerator;

//  SpatialModelCompleteRandomness<CoordType,PixelType> spatialModel;
//  spatialModel.setMask( pixelMatrix );
//  spatialModel.setRandomGenerator( randomGenerator );
//  spatialModel.initialize();

//  DataSet dataSet;
//  const int numObjects = 50;
//  const int numBatches = 100;
//  testModel( spatialModel, spatialModel, numObjectsS1, initSampleSizes(), initNumMonteCarloSamples(), numBatches, dataSet );
//  dataSet.save( "dataSet.data", true );
//  LEAVE();
//}

//void testCSR()
//{
//  ENTER( "void testCSR()" );
//  //testCSROverSquareMask();
//  //testCSROverSquareBoundary();
//  testCSROverSinusoidalBoundary();
//  LEAVE();
//}

void testBoundaryInteractionModel(const float marginProb, DataSet& gDataSet)
{
  ENTER( "void createSpatialModels(const float,DataSet&)" );
  const Curve<CoordType> boundary = initSquareBoundary();
  SpatialModelBoundaryInteraction<float> spatialModel1;
  RandomGenerator randomGenerator;
  spatialModel1.setBoundary( boundary );
  spatialModel1.setRandomGenerator( randomGenerator );
//  spatialModel1.setMargin( 0.1 );
  spatialModel1.setMargin( 0.4 );
  spatialModel1.setMarginProb( marginProb );
  spatialModel1.initialize();

  SpatialModelBoundaryInteraction<float> spatialModel2;
  spatialModel2.setBoundary( boundary );
  spatialModel2.setRandomGenerator( randomGenerator );
  //  spatialModel2.setMargin( 0.1 );
  spatialModel2.setMargin( 0.4 );
  spatialModel2.setMarginProb( marginProb );
  spatialModel2.initialize();


  SpatialModelCompleteRandomness<CoordType,PixelType> csrSpatialModel;
  csrSpatialModel.setBoundary( boundary );
  csrSpatialModel.setRandomGenerator( randomGenerator );
  csrSpatialModel.initialize();

  boundary.save( "squareBoundary", true );
  spatialModel1.getInnerBoundary().save( "innerSquareBoundary1", true );
  spatialModel1.drawSample( 50 ).save( "pattern-squareBoundaryInteractionModel1", true );

  //const int numObjectsS2 = initNumObjectsS2();
  //spatialModel2.getInnerBoundary().save( "innerSquareBoundary2", true );
  //spatialModel2.drawSample( numObjectsS2 ).save( "pattern-squareBoundaryInteractionModel2", true );

  DataSet dataSet;
  const int numObjectsS1 = 50;
//  const vector<int> sampleSizes(1,100);
  const int numBatches = 1;
  testModel( spatialModel1, spatialModel2, csrSpatialModel, numObjectsS1, initNumObjectsS2(), boundary, numBatches, dataSet );
  dataSet.setValues( "marginProb", marginProb );
  gDataSet.append( dataSet );
  LEAVE();
}

void powerUnderBoundaryInteractionModel()
{
  ENTER( "void powerUnderBoundaryInteractionModel()" );

  DataSet dataSet;
  const string id = "power-peripheral-square-objects100";
  for (int i = 0; i < 24; ++i)
  {
    const float marginProb = 0.04 + i * 0.04;
    EVAL( marginProb );
    testBoundaryInteractionModel( marginProb, dataSet );
  }
  dataSet.save( id+".data", true );

  LEAVE();
}

void testsStatisticalTests()
{
  ENTER( "void testSpatialAnalysis()" );
  //testPatternGeneration();
  //testCSR();
  powerUnderBoundaryInteractionModel();
  LEAVE();
}