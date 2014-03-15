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

template class SpatialDescriptorDistanceToBorder<float>;
template class SpatialDescriptorDistanceToBorder<double>;
template class SpatialDescriptorDistanceToBorder<long double>;
