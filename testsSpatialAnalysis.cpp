#include <ask4.h>
#include <dataset.h>
#include <pixelmatrix.h>
#include <spatialdescriptorfunctionf.h>
#include <spatialdescriptorfunctiong.h>
#include <spatialdescriptorfunctionh.h>
#include "spatialmodelboundaryinteraction.h"
#include <spatialmodelcompleterandomness.h>
#include <spatialmodelevaluator.h>

#include <cmath>

#define CoordType float
#define PixelType float

#define TRACE
#include <trace.h>

vector<int> initNumMonteCarloSamples()
{
  vector<int> numMonteCarloSamples;
  numMonteCarloSamples.push_back( 9 );
  numMonteCarloSamples.push_back( 19 );
  numMonteCarloSamples.push_back( 39 );
  return numMonteCarloSamples;
}

vector<int> initSampleSizes()
{
  vector<int> sampleSizes;
  sampleSizes.push_back( 10 );
  for (int s = 50; s <= 1000; s += 50)
    sampleSizes.push_back( s );
  return sampleSizes;
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

void testPatternGeneration(
  SpatialModel<CoordType,PixelType>& spatialModel,
  const int numObjects,
  const string& id)
{
  Vertices<CoordType> pattern;
  pattern = spatialModel.drawSample( numObjects );
  pattern.save( "pattern-"+id, true );

  EVAL( pattern.squareNearestNeighborDistances() );
  EVAL( pattern.squareNearestNeighborDistances().min() );

  SpatialDescriptorFunctionF<CoordType> functionF;
  SpatialDescriptorFunctionG<CoordType> functionG;
  SpatialDescriptorFunctionH<CoordType> functionH;
  SpatialModelEvaluator<CoordType,PixelType> spatialModelEvaluator;
  spatialModelEvaluator.setModel( spatialModel );
  spatialModelEvaluator.setNumMonteCarloSamples( 39 );
  spatialModelEvaluator.addDescriptor( functionF );
  spatialModelEvaluator.addDescriptor( functionG );
  spatialModelEvaluator.addDescriptor( functionH );
  spatialModelEvaluator.setPrecision( 1.0 );
  functionF.setEvaluationPositions( spatialModel.drawSample(1000) );

  vector<float> pValues;
  vector<int> ranks;
  spatialModelEvaluator.eval( pattern, pValues, ranks );
  EVAL( pValues[0] );
  EVAL( pValues[1] );
  EVAL( pValues[2] );
}


void testPatternGeneration()
{
  SpatialModelCompleteRandomness<CoordType,PixelType> spatialModel;
  const Curve<CoordType> boundary = initSinusoidalBoundary();
  RandomGenerator randomGenerator;
  const int numObjects = ask4int( "Number of objects ");

  spatialModel.setRandomGenerator( randomGenerator );
  spatialModel.setBoundary( boundary );
  spatialModel.initialize();
  testPatternGeneration( spatialModel, numObjects, "sinusoidal" );
  boundary.save( "sinusoidalBoundary", true );
}

void testPattern(
  const Vertices<CoordType>& pattern,
  SpatialModel<CoordType,PixelType>& spatialModel,
  SpatialModel<CoordType,PixelType>& csrModel,
  const int numMonteCarloSamples,
  DataSet& gDataSet)
{
  SpatialDescriptorFunctionF<CoordType> functionF;
  SpatialDescriptorFunctionG<CoordType> functionG;
  SpatialDescriptorFunctionH<CoordType> functionH;
  SpatialModelEvaluator<CoordType,PixelType> spatialModelEvaluator;
  spatialModelEvaluator.setModel( spatialModel );
  spatialModelEvaluator.setNumMonteCarloSamples( numMonteCarloSamples );
  spatialModelEvaluator.addDescriptor( functionF );
  spatialModelEvaluator.addDescriptor( functionG );
  spatialModelEvaluator.addDescriptor( functionH );
  spatialModelEvaluator.setPrecision( 1.0 );
  functionF.setEvaluationPositions( csrModel.drawSample(1000) );

  const int r = gDataSet.numRows();
  vector<float> pValues;
  vector<int> ranks;
  spatialModelEvaluator.eval( pattern, pValues, ranks );
  gDataSet.setValue( "F-SDI", r, pValues[0] );
  gDataSet.setValue( "G-SDI", r, pValues[1] );
  gDataSet.setValue( "H-SDI", r, pValues[2] );
}

void testModel(
  SpatialModel<CoordType,PixelType>& spatialModel,
  SpatialModel<CoordType,PixelType>& csrModel,
  const int numObjects,
  const int sampleSize,
  const int numMonteCarloSamples,
  DataSet& gDataSet)
{
  ENTER( "void testModel(...)" );
  EVAL( numMonteCarloSamples );
  EVAL( sampleSize );
  for (int i = 0; i < sampleSize; ++i)
  {
    const Vertices<CoordType> pattern = spatialModel.drawSample( numObjects );
    pattern.save( "current-pattern", true );
    testPattern( pattern, csrModel, csrModel, numMonteCarloSamples, gDataSet );
  }
  LEAVE();
}

void testModel(
  SpatialModel<CoordType,PixelType>& spatialModel,
  SpatialModel<CoordType,PixelType>& csrModel,
  const int numObjects,
  const int sampleSize,
  const int numMonteCarloSamples,
  const int numBatches,
  DataSet& gDataSet)
{
  ENTER( "void testModel(...)" );  
  for (int b = 1; b <= numBatches; ++b)
  {
    EVAL( b );
    DataSet dataSet;
    testModel( spatialModel, csrModel, numObjects, sampleSize, numMonteCarloSamples, dataSet );
    dataSet.setValues( "batch", b );
    gDataSet.append( dataSet );
  }
  LEAVE();
}

void testModel(
  SpatialModel<CoordType,PixelType>& spatialModel,
  SpatialModel<CoordType,PixelType>& csrModel,
  const int numObjects,
  const vector<int>& sampleSizes,
  const vector<int>& numMonteCarloSamples,
  const int numBatches,
  DataSet& gDataSet)
{
  for (size_t s = 0; s < sampleSizes.size(); ++s)
    for (size_t m = 0; m < numMonteCarloSamples.size(); ++m)
    {
      DataSet dataSet;
      testModel( spatialModel, csrModel, numObjects, sampleSizes[s], numMonteCarloSamples[m], numBatches, dataSet );
      dataSet.setValues( "sampleSize", sampleSizes[s] );
      dataSet.setValues( "numMonteCarloSamples", numMonteCarloSamples[m] );
      gDataSet.append( dataSet );
      gDataSet.save( "current-dataset.data", true );
    }
}

void testCSROverBoundary(const Curve<CoordType>& boundary, const string& id)
{
  ENTER( "void testCSROverBoundary(...)" );
  RandomGenerator randomGenerator;
  SpatialModelCompleteRandomness<CoordType,PixelType> spatialModel;
  spatialModel.setBoundary( boundary );
  spatialModel.setRandomGenerator( randomGenerator );
  spatialModel.initialize();

  DataSet dataSet;
  const int numObjects = 10;
  const int numBatches = 30;
  testModel( spatialModel, spatialModel, numObjects, initSampleSizes(), initNumMonteCarloSamples(), numBatches, dataSet );
  dataSet.save( id+".data", true );
  LEAVE();
}

void testCSROverSquareBoundary()
{
  ENTER( "void testCSROverSquareBoundary()" );
  const Curve<CoordType> boundary = initSquareBoundary();
  boundary.save( "squareBoundary", true );
  testCSROverBoundary( boundary, "testCSROverSquareBoundary" );
  LEAVE();
}

void testCSROverSinusoidalBoundary()
{
  ENTER( "void testCSROverSinusoidalBoundary()" );
  const Curve<CoordType> boundary = initSinusoidalBoundary();
  boundary.save( "sinusoidalBoundary", true );
  testCSROverBoundary( boundary, "testCSROverSinusoidalBoundary" );
  LEAVE();
}

void testCSROverSquareMask()
{
  ENTER( "void testCSROverSquare()" );
  const PixelMatrix<PixelType> pixelMatrix = initSquareMask();
  RandomGenerator randomGenerator;

  SpatialModelCompleteRandomness<CoordType,PixelType> spatialModel;
  spatialModel.setMask( pixelMatrix );
  spatialModel.setRandomGenerator( randomGenerator );
  spatialModel.initialize();

  DataSet dataSet;
  const int numObjects = 10;
  const int numBatches = 10;
  testModel( spatialModel, spatialModel, numObjects, initSampleSizes(), initNumMonteCarloSamples(), numBatches, dataSet );
  dataSet.save( "dataSet.data", true );
  LEAVE();
}

void testCSR()
{
  ENTER( "void testCSR()" );
  //testCSROverSquareMask();
  //testCSROverSquareBoundary();
  testCSROverSinusoidalBoundary();
  LEAVE();
}

void testBoundaryInteractionModel(const float marginProb, DataSet& gDataSet)
{
  ENTER( "void testBoundaryInteractionModel(const float,DataSet&)" );
  const Curve<CoordType> boundary = initSquareBoundary();
  SpatialModelBoundaryInteraction<float> spatialModel;
  RandomGenerator randomGenerator;
  spatialModel.setBoundary( boundary );
  spatialModel.setRandomGenerator( randomGenerator );
  spatialModel.setMargin( 0.1 );
  spatialModel.setMarginProb( marginProb );
  spatialModel.initialize();

  SpatialModelCompleteRandomness<CoordType,PixelType> csrSpatialModel;
  csrSpatialModel.setBoundary( boundary );
  csrSpatialModel.setRandomGenerator( randomGenerator );
  csrSpatialModel.initialize();

  boundary.save( "squareBoundary", true );
  spatialModel.getInnerBoundary().save( "innerSquareBoundary", true );
  spatialModel.drawSample( 100 ).save( "pattern-squareBoundaryInteractionModel", true );

  DataSet dataSet;
  const int numObjects = 100;
  const vector<int> sampleSizes(1,50);
  const vector<int> numMonteCarloSamples(1,39);
  const int numBatches = 50;
  testModel( spatialModel, csrSpatialModel, numObjects, sampleSizes, numMonteCarloSamples, numBatches, dataSet );
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

void testsSpatialAnalysis()
{
  ENTER( "void testSpatialAnalysis()" );
  //testPatternGeneration();
  //testCSR();
  powerUnderBoundaryInteractionModel();
  LEAVE();
}
