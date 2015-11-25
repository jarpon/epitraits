/*!
 * \class  SpatialModelBorderDistance3DDifferentCompartments
 * \author Javier Arpón (ja), INRA
 * \author Philippe Andrey (pa), INRA
 * \date   XXXX.XX.XX - creation (ja)
 * \date   2015.10.12 - integration (pa)
 * \brief  3D point process with distance-to-border constraint
****************************************************************/

#include "spatialmodelborderdistance3ddifferentcompartments.h"

#include <programerror.h>

#include <cmath>

/*! Constructor.
****************************************************************/
template<class CoordType>
SpatialModelBorderDistance3DDifferentCompartments<CoordType>::SpatialModelBorderDistance3DDifferentCompartments() : TriMeshSpatialModelDifferentCompartments<CoordType>()
{
}

/*! Sets the distances between objects and border.
****************************************************************/
template<class CoordType>
void SpatialModelBorderDistance3DDifferentCompartments<CoordType>::setDistancesToBorder(const Vector<CoordType>& distancesToBorder)
{
  _distancesToBorder = distancesToBorder;
}

/*! Returns the distances between objects and border.
****************************************************************/
template<class CoordType>
const Vector<CoordType>& SpatialModelBorderDistance3DDifferentCompartments<CoordType>::getDistancesToBorder() const
{
  return _distancesToBorder;
}

/*! Generates vertices into the trimesh with a fixed distance to the border.
****************************************************************/
template<class CoordType>
Vertices<CoordType> SpatialModelBorderDistance3DDifferentCompartments<CoordType>::drawSample(const int numVertices)
{
  if ( numVertices != _distancesToBorder.getSize() )
  {
    ProgramError programError;
    programError.setWhere( "void SpatialModelBorderDistance3DDifferentCompartments<CoordType>::drawSample(const int)" );
    programError.setWhat( "The number of objects differs from the number of distances to the border" );
    throw programError;
  }

  Vertices<CoordType> vertices( 3, numVertices );
  for (int v = 0; v < numVertices; ++v)
    drawPositionFromBorder( vertices[v] , _distancesToBorder[v] );
  return vertices;
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

template class SpatialModelBorderDistance3DDifferentCompartments<float>;
template class SpatialModelBorderDistance3DDifferentCompartments<double>;
template class SpatialModelBorderDistance3DDifferentCompartments<long double>;
