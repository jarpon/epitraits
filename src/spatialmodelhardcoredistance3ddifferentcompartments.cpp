/*!
 * \class  SpatialModelHardcoreDistance3DDifferentCompartments
 * \author Javier Arpón (ja), INRA
 * \author Philippe Andrey (pa), INRA
 * \date   XXXX.XX.XX - creation (ja)
 * \date   2015.10.12 - integration (pa)
 * \brief  3D hardcore point process
****************************************************************/

#include "spatialmodelhardcoredistance3ddifferentcompartments.h"

#include <programerror.h>

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
void SpatialModelHardcoreDistance3DDifferentCompartments<CoordType>::setHardcoreDistances(const Vector<CoordType>& hardcoreDistances)
{
  _hardcoreDistances = hardcoreDistances;
}

/*! Returns the hardcore distances.
****************************************************************/
template<class CoordType>
const Vector<CoordType>& SpatialModelHardcoreDistance3DDifferentCompartments<CoordType>::getHardcoreDistances() const
{
  return _hardcoreDistances;
}

/*! Generates a sample of points respecting the distance constraints.
****************************************************************/
template<class CoordType>
Vertices<CoordType> SpatialModelHardcoreDistance3DDifferentCompartments<CoordType>::drawSample(const int numVertices)
{
  ENTER( "Vertices<CoordType> SpatialModelHardcoreDistance3DDifferentCompartments<CoordType>::drawSample(const int)" );

  if ( numVertices != _hardcoreDistances.getSize() )
  {
    ProgramError programError;
    programError.setWhere( "void SpatialModelHardcoreDistance3DDifferentCompartments<CoordType>::drawSample(const int)" );
    programError.setWhat( "The number of objects differs from the number of hardcore distances" );
    throw programError;
  }

  const int numPermutations = numVertices; // this is arbitrary but works...
  const int maxAttempts = 200;
  Vertices<CoordType> vertices( 3, 0 );
  Vector<CoordType> vertex( 3 );
  int attempts;
  bool success = false;

  for (int p = 0; p < numPermutations && !success; ++p)
  {
    shuffleDistances( _hardcoreDistances );
    vertices.setSize( 0 );
    success = true;

    for (int v = 0; v < numVertices; ++v)
    {
      attempts = 0;
      do {
        this->drawPosition( vertex );
      } while (
               ( !validInterObjectDistances(vertex,vertices,_hardcoreDistances)
                 || !validObjectToBorderDistance(vertex,_hardcoreDistances[v]) )
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

  return vertices;
}

template<class CoordType>
void SpatialModelHardcoreDistance3DDifferentCompartments<CoordType>::shuffleDistances(Vector<CoordType>& distances)
{
  const int n = distances.getSize();
  int i1, i2, i;
  CoordType tmp;

  for (i = 0; i < n; ++i)
  {
    i1 = this->getRandomGenerator().uniformL( n );
    i2 = this->getRandomGenerator().uniformL( n );
    tmp = distances[i1];
    distances[i1] = distances[i2];
    distances[i2] = tmp;
  }
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
