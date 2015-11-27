/*!
 * \class  SpatialModelHardcoreDistance3DDifferentCompartments
 * \author Javier Arpón (ja), INRA
 * \date   2015.11.26 - creation (ja)
 * \brief  3D hardcore point process of two kind
****************************************************************/

#include "spatialmodelhardcoredistance3ddifferentcompartments.h"

#include <programerror.h>
#include <cmath>

//#define TRACE
#include <trace.h>

/*! Constructor.
****************************************************************/
template<class CoordType>
SpatialModelHardcoreDistance3DDifferentCompartments<CoordType>::SpatialModelHardcoreDistance3DDifferentCompartments() : TriMeshSpatialModelDifferentCompartments<CoordType>()
{
}

/*! Sets the hardcore distances.
****************************************************************/
template<class CoordType>
void SpatialModelHardcoreDistance3DDifferentCompartments<CoordType>::setHardcoreDistances(const Vector<CoordType>& hardcoreDistancesDistribution1, const Vector<CoordType>& hardcoreDistancesDistribution2)
{
  _hardcoreDistancesDistribution1 = hardcoreDistancesDistribution1;
  _hardcoreDistancesDistribution2 = hardcoreDistancesDistribution2;
  _hardcoreDistances = hardcoreDistancesDistribution1;
  _hardcoreDistances.append( hardcoreDistancesDistribution2 );
  EVAL(_hardcoreDistances);
}

/*! Returns the hardcore distances of objects kind 1.
****************************************************************/
template<class CoordType>
const Vector<CoordType>& SpatialModelHardcoreDistance3DDifferentCompartments<CoordType>::getHardcoreDistancesDistribution1() const
{
  return _hardcoreDistancesDistribution1;
}

/*! Returns the hardcore distances of objects kind 2.
****************************************************************/
template<class CoordType>
const Vector<CoordType>& SpatialModelHardcoreDistance3DDifferentCompartments<CoordType>::getHardcoreDistancesDistribution2() const
{
  return _hardcoreDistancesDistribution2;
}

/*! Returns the vertices of objects kind 1.
****************************************************************/
template<class CoordType>
const Vertices<CoordType>& SpatialModelHardcoreDistance3DDifferentCompartments<CoordType>::getVerticesDistribution1() const
{
  return _verticesDist1;
}

/*! Returns the vertices of objects kind 2.
****************************************************************/
template<class CoordType>
const Vertices<CoordType>& SpatialModelHardcoreDistance3DDifferentCompartments<CoordType>::getVerticesDistribution2() const
{
  return _verticesDist2;
}

/*! Generates a sample of points respecting the distance constraints.
****************************************************************/
template<class CoordType>
Vertices<CoordType> SpatialModelHardcoreDistance3DDifferentCompartments<CoordType>::drawSample(const int numVerticesDist1, const int numVerticesDist2)
{
  ENTER( "Vertices<CoordType> SpatialModelHardcoreDistance3DDifferentCompartments<CoordType>::drawSample(const int)" );

  _numVerticesDist1 = numVerticesDist1;
  _numVerticesDist2 = numVerticesDist2;

  if ( ( _numVerticesDist1 != _hardcoreDistancesDistribution1.getSize() ) || ( _numVerticesDist2 != _hardcoreDistancesDistribution2.getSize() ) )
  {
    ProgramError programError;
    programError.setWhere( "void SpatialModelHardcoreDistance3DDifferentCompartments<CoordType>::drawSample(const int)" );
    programError.setWhat( "At least one of the number of objects differs from their number of hardcore distances" );
    throw programError;
  }

  //new part corresponding to 2 compartments
  const int numVertices = numVerticesDist1 + numVerticesDist2;

  //keep order of classes 1 and 2
  _classCoordinates.setSize( numVertices, 2 );

  for ( int i = 0; i < numVerticesDist1; ++i )
  {
    _classCoordinates[i][0] = 1;
    _classCoordinates[i][1] = i;
  }
  for ( int j = numVerticesDist1; j < numVertices; ++j )
  {
    _classCoordinates[j][0] = 2;
    _classCoordinates[j][1] = j;
  }


  const int numPermutations = numVertices; // this is arbitrary but works...
  const int maxAttempts = 200;
  Vertices<CoordType> vertices( 3, 0 );
  Vector<CoordType> vertex( 3 );
  int attempts;
  bool success = false;

  Vector<CoordType> hardcoreDistancesTemp;
  hardcoreDistancesTemp = _hardcoreDistances;

  for (int p = 0; p < numPermutations && !success; ++p)
  {
    shuffleObjectsOrder( hardcoreDistancesTemp );
    vertices.setSize( 0 );
    success = true;

    for (int v = 0; v < numVertices; ++v)
    {
      attempts = 0;
      do {
        this->drawPosition( vertex );
      } while (
               ( !validInterObjectDistances(vertex,vertices,hardcoreDistancesTemp)
                 || !validObjectToBorderDistance(vertex,hardcoreDistancesTemp[v]) )
               && ++attempts < maxAttempts );

      if ( attempts < maxAttempts )
        vertices.append( vertex );
      else
      {
        success = false;
        break;
      }
    }
  }

  if ( !success )
  {
    Exception exception;
    exception.setWhere( "void SpatialModelHardcoreDistance3DDifferentCompartments<CoordType>::drawSample(const int)" );
    exception.setWhat( "Too many unsuccessful attemps to generate a vertex" );
    throw exception;
  }

  LEAVE();

  return sortVertices( vertices, hardcoreDistancesTemp );
}

/*! Generates a sample of points respecting the distance constraints.
****************************************************************/
template<class CoordType>
Vertices<CoordType> SpatialModelHardcoreDistance3DDifferentCompartments<CoordType>::sortVertices( const Vertices<CoordType>& vertices, const Vector<CoordType>& randomOrder )
{
  Vector<CoordType> tempOrder = randomOrder;
  const int numVertices = vertices.getNumVertices();
  int j = 0;
  Vertices<CoordType> tempVerticesDist1 ( 3, 0 );
  Vertices<CoordType> tempVerticesDist2 ( 3, 0 );

  for ( int i = 0; i < numVertices; ++i )
  {
    EVAL(vertices[i]);
  }

  for ( int i = 0; i < numVertices; ++i )
  {
    j = tempOrder.find( _hardcoreDistances[i] );
    EVAL(_hardcoreDistances[i]);
    EVAL(j);
    EVAL(vertices[j]);
    EVAL(_classCoordinates[j][0]);
    if ( _classCoordinates[j][1] < _numVerticesDist1 )
      tempVerticesDist1.append( vertices[j] );
    else
      tempVerticesDist2.append( vertices[j] );

    tempOrder[j] = sqrt(-1);
  }

  _verticesDist1.setSize( _numVerticesDist1 );
  _verticesDist1 = tempVerticesDist1;
  _verticesDist2.setSize( _numVerticesDist2 );
  _verticesDist2 = tempVerticesDist2;

  Vertices<CoordType> sortedVertices;
  sortedVertices = _verticesDist1;
  sortedVertices.append( _verticesDist2 );

  for ( int i = 0; i < numVertices; ++i )
  {
    EVAL(sortedVertices[i]);
  }

  return sortedVertices;
}

template<class CoordType>
void SpatialModelHardcoreDistance3DDifferentCompartments<CoordType>::shuffleObjectsOrder(Vector<CoordType>& distances)
{
  const int n = distances.getSize();
  int i1, i2, i;
  CoordType tmp;
  Vector<int> tmpClass( 2 );

  for (i = 0; i < n; ++i)
  {
    i1 = this->getRandomGenerator().uniformL( n );
    i2 = this->getRandomGenerator().uniformL( n );
    tmp = distances[i1];
    distances[i1] = distances[i2];
    distances[i2] = tmp;

    //save correspondence to the distribution kind
    //at the end we want to split the proper objects in the two correct distributions
    tmpClass = _classCoordinates.getRow( i1 );
    _classCoordinates.setRow( i1, _classCoordinates.getRow( i2 ) );
    _classCoordinates.setRow( i2, tmpClass );

  }
  EVAL( distances );
}

/*! Returns \c true if \c vertex respects the distance constraints with the
 * already present \c vertices.
****************************************************************/
template<class CoordType>
bool SpatialModelHardcoreDistance3DDifferentCompartments<CoordType>::validInterObjectDistances(
  const Vector<CoordType>& vertex,
  const Vertices<CoordType>& vertices,
  const Vector<CoordType>& minimumDistances)
{
  const int k = vertices.getNumVertices();

  for (int i = 0; i < k; ++i)
    if ( vertices[i].distance(vertex) < minimumDistances[i]+minimumDistances[k] )
      return false;

  return true;
}

/*! Returns \c true if \c vertex is farther from the border than the specified distance.
****************************************************************/
template<class CoordType>
bool SpatialModelHardcoreDistance3DDifferentCompartments<CoordType>::validObjectToBorderDistance(
  const Vector<CoordType>& vertex,
  const CoordType& minimumDistance)
{
  Vector<CoordType> triMeshVertex; // dummy
  return this->getTriMeshQuery().closestPoint(vertex,triMeshVertex) >= minimumDistance;
}

template class SpatialModelHardcoreDistance3DDifferentCompartments<float>;
template class SpatialModelHardcoreDistance3DDifferentCompartments<double>;
template class SpatialModelHardcoreDistance3DDifferentCompartments<long double>;
