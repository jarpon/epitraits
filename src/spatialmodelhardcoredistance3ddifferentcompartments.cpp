/*!
 * \class  SMHardcoreDistance3DDifferentCompartments
 * \author Javier Arpón (ja), INRA
 * \author Philippe Andrey (pa), INRA
 * \date   XXXX.XX.XX - creation (ja)
 * \date   2015.10.12 - integration (pa)
 * \brief  3D hardcore point process
****************************************************************/

#include "spatialmodelhardcoredistance3ddifferentcompartments2.h"

#include <programerror.h>

#define TRACE
#include <trace.h>

/*! Constructor.
****************************************************************/
template<class CoordType>
SMHardcoreDistance3DDifferentCompartments<CoordType>::SMHardcoreDistance3DDifferentCompartments() : TriMeshSpatialModelDifferentCompartments<CoordType>()
{
}

/*! Sets the hardcore distances.
****************************************************************/
template<class CoordType>
void SMHardcoreDistance3DDifferentCompartments<CoordType>::setHardcoreDistances(const Vector<CoordType>& hardcoreDistancesDistribution1, const Vector<CoordType>& hardcoreDistancesDistribution2)
{
  _hardcoreDistancesDistribution1 = hardcoreDistancesDistribution1;
  _hardcoreDistancesDistribution2 = hardcoreDistancesDistribution2;
  _hardcoreDistances = hardcoreDistancesDistribution1;
  _hardcoreDistances.append( hardcoreDistancesDistribution2 );
}

/*! Returns the hardcore distances of objects kind 1.
****************************************************************/
template<class CoordType>
const Vector<CoordType>& SMHardcoreDistance3DDifferentCompartments<CoordType>::getHardcoreDistancesDistribution1() const
{
  return _hardcoreDistancesDistribution1;
}

/*! Returns the hardcore distances of objects kind 2.
****************************************************************/
template<class CoordType>
const Vector<CoordType>& SMHardcoreDistance3DDifferentCompartments<CoordType>::getHardcoreDistancesDistribution2() const
{
  return _hardcoreDistancesDistribution2;
}

/*! Generates a sample of points respecting the distance constraints.
****************************************************************/
template<class CoordType>
Vertices<CoordType> SMHardcoreDistance3DDifferentCompartments<CoordType>::drawSample(const int numVerticesDist1, const int numVerticesDist2)
{
  ENTER( "Vertices<CoordType> SMHardcoreDistance3DDifferentCompartments<CoordType>::drawSample(const int)" );

  if ( ( numVerticesDist1 != _hardcoreDistancesDistribution1.getSize() ) || ( numVerticesDist2 != _hardcoreDistancesDistribution2.getSize() ) )
  {
    ProgramError programError;
    programError.setWhere( "void SMHardcoreDistance3DDifferentCompartments<CoordType>::drawSample(const int)" );
    programError.setWhat( "At least one of the number of objects differs from their number of hardcore distances" );
    throw programError;
  }

  //new part corresponding to 2 compartments
  const int numVertices = numVerticesDist1 + numVerticesDist2;
  //class 1 and 2
//  _classBelongings.setSize( numVertices );
//  _classBelongings.setOnes();
//  for ( int i = numVerticesDist1; i < numVertices; ++i )
//    _classBelongings[i] = 2;

  //keep order of classes 1 and 2
  _classBelongings.setSize( numVertices );
  for ( int i = 0; i < numVerticesDist1; ++i )
    _classBelongings[i] = i;
  for ( int j = numVerticesDist1; j < numVertices; ++j )
    _classBelongings[j] = j;
  EVAL(_classBelongings);

//  Vertices<int> _classCoordinates( 2, numVertices, 0, 0 );
  Matrix<int> _classCoordinates( numVertices, 2 );

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
    exception.setWhere( "void SMHardcoreDistance3DDifferentCompartments<CoordType>::drawSample(const int)" );
    exception.setWhat( "Too many unsuccessful attemps to generate a vertex" );
    throw exception;
  }

  LEAVE();

  return vertices;
}

template<class CoordType>
void SMHardcoreDistance3DDifferentCompartments<CoordType>::shuffleDistances(Vector<CoordType>& distances)
{
  const int n = distances.getSize();
  int i1, i2, i;
  CoordType tmp;
  int o1, o2, tmp2;


  for (i = 0; i < n; ++i)
  {
    i1 = this->getRandomGenerator().uniformL( n );
    i2 = this->getRandomGenerator().uniformL( n );
    tmp = distances[i1];
    distances[i1] = distances[i2];
    distances[i2] = tmp;

    //save correspondence to the distribution kind
    //at the end we want to split the proper objects in the two correct distributions
    tmp2 = _classBelongings[i1];
    _classBelongings[i1] = _classBelongings[i2];
    _classBelongings[i2] = tmp2;


  }
}

/*! Returns \c true if \c vertex respects the distance constraints with the
 * already present \c vertices.
****************************************************************/
template<class CoordType>
bool SMHardcoreDistance3DDifferentCompartments<CoordType>::validInterObjectDistances(
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
bool SMHardcoreDistance3DDifferentCompartments<CoordType>::validObjectToBorderDistance(
  const Vector<CoordType>& vertex,
  const CoordType& minimumDistance)
{
  Vector<CoordType> triMeshVertex; // dummy
  return this->getTriMeshQuery().closestPoint(vertex,triMeshVertex) >= minimumDistance;
}

template class SMHardcoreDistance3DDifferentCompartments<float>;
template class SMHardcoreDistance3DDifferentCompartments<double>;
template class SMHardcoreDistance3DDifferentCompartments<long double>;
