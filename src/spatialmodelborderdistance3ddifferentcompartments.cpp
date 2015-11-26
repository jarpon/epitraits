/*!
 * \class  SpatialModelHardcoreDistance3DDifferentCompartments
 * \author Javier Arpón (ja), INRA
 * \date   2015.11.26 - creation (ja)
 * \brief  3D point process of two classes with distance-to-border constraint
****************************************************************/

#include "spatialmodelborderdistance3ddifferentcompartments.h"

#include <programerror.h>

#include <cmath>

//#define TRACE
#include <trace.h>

/*! Constructor.
****************************************************************/
template<class CoordType>
SpatialModelBorderDistance3DDifferentCompartments<CoordType>::SpatialModelBorderDistance3DDifferentCompartments() : TriMeshSpatialModelDifferentCompartments<CoordType>()
{
}

/*! Sets the distances between objects and border.
****************************************************************/
template<class CoordType>
void SpatialModelBorderDistance3DDifferentCompartments<CoordType>::setDistancesToBorder(const Vector<CoordType>& distancesToBorderDistribution1, const Vector<CoordType>& distancesToBorderDistribution2)
{
  _distancesToBorderDistribution1 = distancesToBorderDistribution1;
  _distancesToBorderDistribution2 = distancesToBorderDistribution2;
  _distancesToBorder = distancesToBorderDistribution1;
  _distancesToBorder.append( distancesToBorderDistribution2 );
}

/*! Returns the distances between objects and border of objects class 1.
****************************************************************/
template<class CoordType>
const Vector<CoordType>& SpatialModelBorderDistance3DDifferentCompartments<CoordType>::getDistancesToBorderDistribution1() const
{
  return _distancesToBorderDistribution1;
}

/*! Returns the distances between objects and border of objects class 2.
****************************************************************/
template<class CoordType>
const Vector<CoordType>& SpatialModelBorderDistance3DDifferentCompartments<CoordType>::getDistancesToBorderDistribution2() const
{
  return _distancesToBorderDistribution2;
}

/*! Returns the vertices of objects kind 1.
****************************************************************/
template<class CoordType>
const Vertices<CoordType>& SpatialModelBorderDistance3DDifferentCompartments<CoordType>::getVerticesDistribution1() const
{
  return _verticesDist1;
}

/*! Returns the vertices of objects kind 2.
****************************************************************/
template<class CoordType>
const Vertices<CoordType>& SpatialModelBorderDistance3DDifferentCompartments<CoordType>::getVerticesDistribution2() const
{
  return _verticesDist2;
}

/*! Generates vertices into the trimesh with a fixed distance to the border.
****************************************************************/
template<class CoordType>
Vertices<CoordType> SpatialModelBorderDistance3DDifferentCompartments<CoordType>::drawSample(const int numVerticesDist1, const int numVerticesDist2)
{
  ENTER("Vertices<CoordType> SpatialModelBorderDistance3DDifferentCompartments<CoordType>::drawSample(...)");

  //new part corresponding to 2 compartments
  _numVerticesDist1 = numVerticesDist1;
  _numVerticesDist2 = numVerticesDist2;

  const int numVertices = numVerticesDist1 + numVerticesDist2;

  if ( ( _numVerticesDist1 != _distancesToBorderDistribution1.getSize() ) || ( _numVerticesDist2 != _distancesToBorderDistribution2.getSize() ) )
  {
    ProgramError programError;
    programError.setWhere( "void SpatialModelBorderDistance3DDifferentCompartments<CoordType>::drawSample(const int)" );
    programError.setWhat( "The number of objects differs from the number of distances to the border" );
    throw programError;
  }

  Vector<CoordType> borderDistancesTemp;
  borderDistancesTemp = _distancesToBorder;

  shuffleObjects( borderDistancesTemp );

  Vertices<CoordType> vertices( 3, numVertices );
  for (int v = 0; v < numVertices; ++v)
    drawPositionFromBorder( vertices[v] , borderDistancesTemp[v] );


  LEAVE();

  return sortVertices( vertices, borderDistancesTemp );
}

/*! Generates a sample of points respecting the distance constraints.
****************************************************************/
template<class CoordType>
Vertices<CoordType> SpatialModelBorderDistance3DDifferentCompartments<CoordType>::sortVertices( const Vertices<CoordType>& vertices, const Vector<CoordType>& randomOrder )
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
    j = tempOrder.find( _distancesToBorder[i] );
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

/*! Generates a random vertex into the trimesh taking into account a distance to the border.
****************************************************************/
template<class CoordType>
void SpatialModelBorderDistance3DDifferentCompartments<CoordType>::drawPositionFromBorder(Vector<CoordType>& position, const CoordType distanceToBorder)
{
  const int maxAttempts = 200;
  Vector<CoordType> triMeshVertex;
  Vector<CoordType> normalVector;
  int attempts = 0;
  bool found = false;

  do
  {
    this->drawPosition( position );

    // from the random point we find the closest point to it on the triMesh
    // the segment between them will be the normal vector
    // from this normal vector and using the distance to the border we obtain
    // the "position" vertex which will be our wanted vertex
    this->getTriMeshQuery().closestPoint( position, triMeshVertex );

    normalVector = position - triMeshVertex;
    normalVector /= position.distance( triMeshVertex );
    normalVector *= distanceToBorder;
    position = triMeshVertex + normalVector;

    // check if the obtained point is inside the triMesh and if it keeps
    // the condition of being at the specified distance to the border
    if ( this->getTriMeshQuery().contains(position) )
      found = fabs( distanceToBorder-this->getTriMeshQuery().closestPoint(position,triMeshVertex) ) < Vector<CoordType>::epsilon();

  } while ( found == false && ++attempts <= maxAttempts );

  if ( attempts > maxAttempts )
  {
    Exception exception;
    exception.setWhere( "void SpatialModelBorderDistance3DDifferentCompartments<CoordType>::drawPositionFromBorder(...)" );
    exception.setWhat( "Too many unsuccessful attemps to generate a vertex" );
    throw exception;
  }
}

template<class CoordType>
void SpatialModelBorderDistance3DDifferentCompartments<CoordType>::shuffleObjects(Vector<CoordType>& distances)
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
}

template class SpatialModelBorderDistance3DDifferentCompartments<float>;
template class SpatialModelBorderDistance3DDifferentCompartments<double>;
template class SpatialModelBorderDistance3DDifferentCompartments<long double>;
