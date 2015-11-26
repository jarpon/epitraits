/*!
 * \class  TriMeshSpatialModel
 * \author Javier Arpón (ja), INRA
 * \author Philippe Andrey (pa), INRA
 * \date   XXXX.XX.XX - creation (ja)
 * \date   2015.10.12 - integration (pa)
 * \brief  Base class for 3D point processes
****************************************************************/

#include "trimeshspatialmodeldifferentcompartments.h"

#include <programerror.h>

#include <cmath>

//#define TRACE
#include <trace.h>
#include <sstream>

using namespace std;



/*! Initializes all the parameters.
****************************************************************/
template<class CoordType>
TriMeshSpatialModelDifferentCompartments<CoordType>::TriMeshSpatialModelDifferentCompartments()
{
}

/*! Destroys it.
****************************************************************/
template<class CoordType>
TriMeshSpatialModelDifferentCompartments<CoordType>::~TriMeshSpatialModelDifferentCompartments()
{
}

/*! Sets the triMesh to work with.
****************************************************************/
template<class CoordType>
void TriMeshSpatialModelDifferentCompartments<CoordType>::setTriMesh(const TriMesh<CoordType>& triMesh)
{
  _triMeshQuery.setTriMesh( triMesh );
}

template<class CoordType>
const TriMesh<CoordType>& TriMeshSpatialModelDifferentCompartments<CoordType>::getTriMesh() const
{
  return _triMeshQuery.getTriMesh();
}

template<class CoordType>
const TriMeshQuery<CoordType>& TriMeshSpatialModelDifferentCompartments<CoordType>::getTriMeshQuery() const
{
  return _triMeshQuery;
}


/*! Initializes the class.
****************************************************************/
template<class CoordType>
void TriMeshSpatialModelDifferentCompartments<CoordType>::initialize()
{
  _boundingBox = _triMeshQuery.getTriMesh().boundingBox();

#if 0
  //shows the maximum distance to the border that can be used
  float maxDistance = ( boundingBox.max(0) - boundingBox.min(0) )/2;
  for ( int i = 1; i < 3; ++i)
  {
    if ( ( boundingBox.max(i) - boundingBox.min(i) )/2 < maxDistance )
      maxDistance = ( boundingBox.max(i) - boundingBox.min(i) )/2;
  }
  EVAL(maxDistance);
#endif
}

/*! Generates a random vertex into the trimesh.
****************************************************************/
template<class CoordType>
void TriMeshSpatialModelDifferentCompartments<CoordType>::drawPosition(Vector<CoordType>& position)
{
  RandomGenerator& randomGenerator = this->getRandomGenerator();
  const Vector<CoordType> v1 = _boundingBox.getVertex1();
  const Vector<CoordType> v2 = _boundingBox.getVertex2();

  do
  {
    position[X] = randomGenerator.uniformLF( v1[X], v2[X] );
    position[Y] = randomGenerator.uniformLF( v1[Y], v2[Y] );
    position[Z] = randomGenerator.uniformLF( v1[Z], v2[Z] );
  } while ( _triMeshQuery.contains(position) == false );
}


template class TriMeshSpatialModelDifferentCompartments<double>;
template class TriMeshSpatialModelDifferentCompartments<float>;
template class TriMeshSpatialModelDifferentCompartments<long double>;

