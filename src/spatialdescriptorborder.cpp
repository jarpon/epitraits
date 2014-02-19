/*!
 * \class  SpatialDescriptorFunctionF
 * \author Philippe Andrey (pa), INRA
 * \date   2012.01.16 - creation (pa)
 * \brief  F-function statistics for spatial point processes
****************************************************************/

#include "spatialdescriptorborder.h"

#include <cdftools.h>
#include <trimesh.h>

//#define TRACE
#include <trace.h>

template<class CoordType>
SpatialDescriptorDistanceToBorder<CoordType>::SpatialDescriptorDistanceToBorder() : SpatialDescriptor<CoordType>()
{
  _triMesh = 0;
}

template<class CoordType>
void SpatialDescriptorDistanceToBorder<CoordType>::setTriMesh(const TriMesh<CoordType>& triMesh)
{
  _triMesh = &triMesh;
}

template<class CoordType>
void SpatialDescriptorDistanceToBorder<CoordType>::eval(
  const Vertices<CoordType>& vertices,
  Vector<CoordType>& xvalues,
  Vector<CoordType>& yvalues)
{
  ENTER( "void SpatialDescriptorDistanceToBorder<CoordType>::eval(...)" );
  const int numVertices = vertices.getSize();

  Vector<CoordType> triMeshVertex;
  xvalues.setSize( numVertices );

  for (int i = 0; i < numVertices; ++i)
  {
    _triMesh->closestPoint( vertices[i], triMeshVertex );
    xvalues[i] = vertices[i].distance( triMeshVertex );
  }
  xvalues.sort();

  CDFTools<CoordType> cdftools;
  yvalues = cdftools.cdf( xvalues );
  LEAVE();
}

template class SpatialDescriptorDistanceToBorder<float>;
template class SpatialDescriptorDistanceToBorder<double>;
template class SpatialDescriptorDistanceToBorder<long double>;
