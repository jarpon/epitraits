/*!
 * \class  SpatialDescriptorMaximaRepulsion
 * \author Javier Arpón (ja), INRA
 * \date   2014.02.15 - creation (ja)
 * \brief  Study the repulsive position among studied points into a shape.
****************************************************************/

#include "spatialdescriptormaxima.h"

#include <cdftools.h>
#include <trimesh.h>

//#define TRACE
#include <trace.h>

template<class CoordType>
SpatialDescriptorMaximaRepulsion<CoordType>::SpatialDescriptorMaximaRepulsion() : SpatialDescriptor<CoordType>()
{
  _triMesh = 0;
}

template<class CoordType>
void SpatialDescriptorMaximaRepulsion<CoordType>::setTriMesh(const TriMesh<CoordType>& triMesh)
{
  _triMesh = &triMesh;
}

template<class CoordType>
void SpatialDescriptorMaximaRepulsion<CoordType>::eval(
  const Vertices<CoordType>& vertices,
  Vector<CoordType>& xvalues,
  Vector<CoordType>& yvalues)
{
  ENTER( "void SpatialDescriptorMaximaRepulsion<CoordType>::eval(...)" );
    PRINT("here");
  const int numVertices = vertices.getSize();

  Vector<CoordType> testVertex(3);
  Vector<CoordType> triMeshVertex(3);
  xvalues.setSize( numVertices );
  PRINT("here");
  for (int i = 0; i < numVertices; ++i)
  {
    testVertex = vertices[i];
    _triMesh->closestPoint( testVertex, triMeshVertex );
    PRINT("here!!");

    xvalues[i] = testVertex.distance( triMeshVertex );
  }
  xvalues.sort();
  PRINT("here!!!!!!!");

  CDFTools<CoordType> cdftools;
  yvalues = cdftools.cdf( xvalues );
  PRINT("here!!!!!!!!!!!!");

  LEAVE();
}

template class SpatialDescriptorMaximaRepulsion<float>;
template class SpatialDescriptorMaximaRepulsion<double>;
template class SpatialDescriptorMaximaRepulsion<long double>;

