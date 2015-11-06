/*!
 * \class  SpatialModelEvaluator
 * \author Philippe Andrey (pa), INRA
 * \date   2012.01.16 - creation (pa)
 * \brief  Statistical evaluation of spatial point process models
 * \details
 * Example:
 * \code
  SpatialModelEvaluator<CoordType,PixelType> spatialModelEvaluator;
  spatialModelEvaluator.setModel( spatialModel );
  spatialModelEvaluator.setDescriptor( spatialDescriptor );
  spatialModelEvaluator.setNumRandomSamples( numRandomSamples );
  const float pValue = spatialModelEvaluator.eval( vertices );
  * \endcode
****************************************************************/

#include "spatialmodelevaluator2.h"

//#include <cdftools.h>
#include "cdftools2.h"
//#define TRACE
#include <trace.h>

template<class CoordType,class PixelType>
SpatialModelEvaluator<CoordType,PixelType>::SpatialModelEvaluator()
{
  _model = 0;
  _numMonteCarloSamples = 0;
  _precision = 1.0;
}

template<class CoordType,class PixelType>
void SpatialModelEvaluator<CoordType,PixelType>::setModel(SpatialModel<CoordType, PixelType>& model)
{
  _model = &model;
}

template<class CoordType,class PixelType>
void SpatialModelEvaluator<CoordType,PixelType>::setDescriptor(SpatialDescriptor<CoordType>& descriptor)
{
  _descriptors.clear();
  addDescriptor( descriptor );
}


template<class CoordType,class PixelType>
void SpatialModelEvaluator<CoordType,PixelType>::addDescriptor(SpatialDescriptor<CoordType>& descriptor)
{
  _descriptors.push_back( &descriptor );
}

template<class CoordType,class PixelType>
void SpatialModelEvaluator<CoordType,PixelType>::setNumMonteCarloSamples(const int numMonteCarloSamples)
{
  _numMonteCarloSamples = numMonteCarloSamples;
}

template<class CoordType,class PixelType>
void SpatialModelEvaluator<CoordType,PixelType>::setPrecision(const CoordType& precision)
{
  _precision = precision;
}

///*! Evaluates and returns the p-value for the pattern \c vertices.
//****************************************************************/
//template<class CoordType,class PixelType>
//float SpatialModelEvaluator<CoordType,PixelType>::evalArea(
//  const Vertices<CoordType>& vertices,
//  DataSet* dataSet)
//{
////  vector<float> pValue;
////  vector<int> rank;
////  vector<float> maxDiff;
//  vector<float> pValue;

//  evalArea( vertices, pValue, dataSet );
//  //eval( vertices, pValue, rank, dataSet );

////  Vector<float> output( 2 );
////  output[0] = pValue[0];
////  output[1] = maxDiff[0];
////  return pValue[0];
//  return pValue[0];
//}

/*! Evaluates and returns the p-value for the pattern \c vertices.
****************************************************************/
template<class CoordType,class PixelType>
Vector<float> SpatialModelEvaluator<CoordType,PixelType>::evalSDIandMaxDiff(
  const Vertices<CoordType>& vertices,
  DataSet* dataSet)
{
  vector<float> pValue;
  vector<int> rank;
  vector<float> maxDiff;

  evalSDIandMaxDiff( vertices, pValue, rank, maxDiff, dataSet );
  //eval( vertices, pValue, rank, dataSet );

  Vector<float> output( 2 );
  output[0] = pValue[0];
  output[1] = maxDiff[0];

//  return pValue[0];
  return output;
}

/*! Evaluates and returns the p-value for the pattern \c vertices.
****************************************************************/
template<class CoordType,class PixelType>
float SpatialModelEvaluator<CoordType,PixelType>::eval(
  const Vertices<CoordType>& vertices,
  DataSet* dataSet)
{
  vector<float> pValue;
  vector<int> rank;

  eval( vertices, pValue, rank, dataSet );

  return pValue[0];
}

/*! Evaluates the p-value for the pattern \c vertices.
 *
 * The method also provides the rank from which the p-value was computed.
****************************************************************/
template<class CoordType,class PixelType>
void SpatialModelEvaluator<CoordType,PixelType>::eval(
  const Vertices<CoordType>& vertices,
  SpatialDescriptor<CoordType>& descriptor,
  const ShapeSet<CoordType>& monteCarloSamples1,
  const ShapeSet<CoordType>& monteCarloSamples2,
  float& pValue,
  int& rank,
  DataSet* dataSet)
{
  Vector<CoordType> x, y, yAverage;
  Matrix<CoordType> x1, x2;
  CDFTools<CoordType> cdfTools;

  descriptor.eval( monteCarloSamples1, x1 );
  descriptor.eval( monteCarloSamples2, x2 );
  descriptor.eval( vertices, x, y );

  const CoordType max1 = x1.max().max();
  const CoordType max2 = x2.max().max();
  const CoordType max = x.max();
  const CoordType max12 = max1 > max2? max1: max2;
  const CoordType allMax = max > max12? max: max12;
  Vector<CoordType> xEvals( static_cast<int>(allMax/_precision) + 1 );
  xEvals.setSequence( _precision );
  yAverage = cdfTools.average( x1, xEvals );

  rank = cdfTools.rank( x, x2, xEvals, yAverage );
  pValue = 1.0 - static_cast<float>(rank) / (1.0+_numMonteCarloSamples); // OK, checked (pa)

  if ( dataSet != 0 )
  {
    Vector<float> percents( 2 );
    percents[0] = 2.5;
    percents[1] = 97.5;
    Matrix<CoordType> percentiles = cdfTools.percentiles( x2, xEvals, percents );
    dataSet->setValues( "x", xEvals );
    dataSet->setValues( "observed", cdfTools.cdf(x,xEvals) );
    dataSet->setValues( "average", yAverage );
    dataSet->setValues( "lower", percentiles[0] );
    dataSet->setValues( "upper", percentiles[1] );
  }
}

/*! Evaluates the p-values for the pattern \c vertices.
 * A p-value is computed for each descriptor that has been added.
 * The method also provides the ranks from which the p-values were computed.
 *
 * Importantly, this method evaluates the p-values from the different descriptors
 * using the same set of Monte Carlo patterns generated by the model.
 *
 * The vectors of p-values and of ranks are resized to the number of descriptors.
 *
 * If a dataset is specified, it is filled with data from the first descriptor only.
 * This behavior may change in the future (eg vector of DataSet*).
****************************************************************/
template<class CoordType,class PixelType>
void SpatialModelEvaluator<CoordType,PixelType>::eval(
  const Vertices<CoordType>& vertices,
  vector<float>& pValues,
  vector<int>& ranks,
  DataSet* dataSet)
{
  const int numVertices = vertices.getSize();
  ShapeSet<CoordType> monteCarloSamples1;
  ShapeSet<CoordType> monteCarloSamples2;

  pValues.resize( _descriptors.size() );
  ranks.resize( _descriptors.size() );

  monteCarloSamples1 = _model->drawSamples( _numMonteCarloSamples, numVertices );
  monteCarloSamples2 = _model->drawSamples( _numMonteCarloSamples, numVertices );
  for (size_t i = 0; i < _descriptors.size(); ++i)
    eval( vertices, *_descriptors[i], monteCarloSamples1, monteCarloSamples2, pValues[i], ranks[i], i==0?dataSet:0 );
  monteCarloSamples1.clear();
  monteCarloSamples2.clear();
}

/*!Same getting also the maximum signed difference between the curves
 * Evaluates the p-value for the pattern \c vertices.
 *
 * The method also provides the rank from which the p-value was computed.
****************************************************************/
template<class CoordType,class PixelType>
void SpatialModelEvaluator<CoordType,PixelType>::eval(
  const Vertices<CoordType>& vertices,
  SpatialDescriptor<CoordType>& descriptor,
  const ShapeSet<CoordType>& monteCarloSamples1,
  const ShapeSet<CoordType>& monteCarloSamples2,
  float& pValue,
  int& rank,
  float& maxDiff,
  DataSet* dataSet)
{
  Vector<CoordType> x, y, yAverage;
  Matrix<CoordType> x1, x2;
  CDFTools<CoordType> cdfTools;

  EVAL("here2a1");
  descriptor.eval( monteCarloSamples1, x1 );
  EVAL("here2a2");
  descriptor.eval( monteCarloSamples2, x2 );
  EVAL("here2a3");
  descriptor.eval( vertices, x, y );

  const CoordType max1 = x1.max().max();
  const CoordType max2 = x2.max().max();
  const CoordType max = x.max();
  const CoordType max12 = max1 > max2? max1: max2;
  const CoordType allMax = max > max12? max: max12;
  Vector<CoordType> xEvals( static_cast<int>(allMax/_precision) + 1 );
  xEvals.setSequence( _precision );
  yAverage = cdfTools.average( x1, xEvals );

  //rank = cdfTools.rank( x, x2, xEvals, yAverage );
  Vector<float> output = cdfTools.rankAndMaxDiff( x, x2, xEvals, yAverage );
  rank = output[0];
  maxDiff = output[1];

  pValue = 1.0 - static_cast<float>(rank) / (1.0+_numMonteCarloSamples); // OK, checked (pa)

  if ( dataSet != 0 )
  {
    Vector<float> percents( 2 );
    percents[0] = 2.5;
    percents[1] = 97.5;
    Matrix<CoordType> percentiles = cdfTools.percentiles( x2, xEvals, percents );
    dataSet->setValues( "x", xEvals );
    dataSet->setValues( "observed", cdfTools.cdf(x,xEvals) );
    dataSet->setValues( "average", yAverage );
    dataSet->setValues( "lower", percentiles[0] );
    dataSet->setValues( "upper", percentiles[1] );
  }
}

/*!Same getting also the maximum signed difference between the curves
 * Evaluates the p-value for the pattern \c vertices.
 *
 * The method also provides the rank from which the p-value was computed.
****************************************************************/
template<class CoordType,class PixelType>
void SpatialModelEvaluator<CoordType,PixelType>::evalArea(
  const Vertices<CoordType>& vertices,
  SpatialDescriptor<CoordType>& descriptor,
  const ShapeSet<CoordType>& monteCarloSamples1,
  const ShapeSet<CoordType>& monteCarloSamples2,
  float& pValue,
  Vector<CoordType>& values,
  DataSet* dataSet)
{
  Vector<CoordType> x, y, yAverage;
  Matrix<CoordType> x1, x2;
  CDFTools<CoordType> cdfTools;

  descriptor.eval( monteCarloSamples1, x1 );
  descriptor.eval( monteCarloSamples2, x2 );
  descriptor.eval( vertices, x, y );

  const CoordType max1 = x1.max().max();
  const CoordType max2 = x2.max().max();
  const CoordType max = x.max();
  const CoordType max12 = max1 > max2? max1: max2;
  const CoordType allMax = max > max12? max: max12;
  Vector<CoordType> xEvals( static_cast<int>(allMax/_precision) + 1 );
  xEvals.setSequence( _precision );
  yAverage = cdfTools.average( x1, xEvals );

  // [0] area1
  // [1] area2
  // [2] is the difference between areas (could be more than 1)
  // [3] is the difference between areas^2
  // [4] is the coefficient between areas
  // [5] is the coefficient between areas^2
  Vector<CoordType> allValues;
  allValues = cdfTools.areasDifference( cdfTools.cdf(x), x, yAverage, xEvals );

  pValue = allValues[4];

  if ( values.getSize() != 0 )
  {
    values[0] = allValues[0];
    values[1] = allValues[1];
  }
  EVAL(values[0]);
  EVAL(values[1]);
  if ( dataSet != 0 )
  {
    Vector<float> percents( 2 );
    percents[0] = 2.5;
    percents[1] = 97.5;
    Matrix<CoordType> percentiles = cdfTools.percentiles( x2, xEvals, percents );
    dataSet->setValues( "x", xEvals );
    dataSet->setValues( "observed", cdfTools.cdf(x,xEvals) );
    dataSet->setValues( "average", yAverage );
    dataSet->setValues( "lower", percentiles[0] );
    dataSet->setValues( "upper", percentiles[1] );
  }
}

/*! Evaluates the p-values for the pattern \c vertices.
 * A p-value is computed for each descriptor that has been added.
 * The method also provides the ranks from which the p-values were computed.
 *
 * Importantly, this method evaluates the p-values from the different descriptors
 * using the same set of Monte Carlo patterns generated by the model.
 *
 * The vectors of p-values and of ranks are resized to the number of descriptors.
 *
 * If a dataset is specified, it is filled with data from the first descriptor only.
 * This behavior may change in the future (eg vector of DataSet*).
****************************************************************/
template<class CoordType,class PixelType>
void SpatialModelEvaluator<CoordType,PixelType>::evalArea(
  const Vertices<CoordType>& vertices,
  vector<float>& pValues,
  Matrix<CoordType>& values,
  DataSet* dataSet)
{
  const int numVertices = vertices.getSize();
  ShapeSet<CoordType> monteCarloSamples1;
  ShapeSet<CoordType> monteCarloSamples2;

  pValues.resize( _descriptors.size() );
  if ( ( values.getSize1() != _descriptors.size() ) || ( values.getSize1() != 2 ) )
    values.setSize( _descriptors.size(), 2 );

  monteCarloSamples1 = _model->drawSamples( _numMonteCarloSamples, numVertices );
  monteCarloSamples2 = _model->drawSamples( _numMonteCarloSamples, numVertices );
  EVAL(_descriptors.size());
  for (size_t i = 0; i < _descriptors.size(); ++i)
  {
    //evalArea( vertices, *_descriptors[i], monteCarloSamples1, monteCarloSamples2, pValues[i], tempValues[i], i==0?dataSet:0 );
    evalArea( vertices, *_descriptors[i], monteCarloSamples1, monteCarloSamples2, pValues[i], values[i], i==0?dataSet:0 );
    //values.setRow( i, tempValues[i] );
  }

  monteCarloSamples1.clear();
  monteCarloSamples2.clear();
}

/*! Evaluates the p-values for the pattern \c vertices.
 * A p-value is computed for each descriptor that has been added.
 * The method also provides the ranks from which the p-values were computed.
 *
 * Importantly, this method evaluates the p-values from the different descriptors
 * using the same set of Monte Carlo patterns generated by the model.
 *
 * The vectors of p-values and of ranks are resized to the number of descriptors.
 *
 * If a dataset is specified, it is filled with data from the first descriptor only.
 * This behavior may change in the future (eg vector of DataSet*).
****************************************************************/
template<class CoordType,class PixelType>
void SpatialModelEvaluator<CoordType,PixelType>::evalSDIandMaxDiff(
  const Vertices<CoordType>& vertices,
  vector<float>& pValues,
  vector<int>& ranks,
  vector<float>& maxDiff,
  DataSet* dataSet)
{
  const int numVertices = vertices.getSize();
  ShapeSet<CoordType> monteCarloSamples1;
  ShapeSet<CoordType> monteCarloSamples2;

  pValues.resize( _descriptors.size() );
  ranks.resize( _descriptors.size() );
  maxDiff.resize( _descriptors.size() );

  monteCarloSamples1 = _model->drawSamples( _numMonteCarloSamples, numVertices );
  monteCarloSamples2 = _model->drawSamples( _numMonteCarloSamples, numVertices );
  EVAL(_descriptors.size());
  for (size_t i = 0; i < _descriptors.size(); ++i)
    eval( vertices, *_descriptors[i], monteCarloSamples1, monteCarloSamples2, pValues[i], ranks[i], maxDiff[i], i==0?dataSet:0 );
  EVAL("here3");

  monteCarloSamples1.clear();
  monteCarloSamples2.clear();
}

template class SpatialModelEvaluator<float,float>;
template class SpatialModelEvaluator<double,double>;
template class SpatialModelEvaluator<long double,long double>;
