/*!
 * \class  SpatialDescriptorDistanceToCentroid
 * \author Javier Arpón (ja), INRA
 * \date   2015.04.13 - creation (ja)
 * \brief  Study of the central positioning of objects within a volume.
****************************************************************/

#include "spatialdescriptorcentroid.h"

#include <cdftools.h>
#include <trimesh.h>

//#define TRACE
#include <trace.h>

template<class CoordType>
SpatialDescriptorDistanceToCentroid<CoordType>::SpatialDescriptorDistanceToCentroid() : SpatialDescriptor<CoordType>()
{
  _triMesh = 0;
}

template<class CoordType>
void SpatialDescriptorDistanceToCentroid<CoordType>::setTriMesh(const TriMesh<CoordType>& triMesh)
{
  _triMesh = &triMesh;
}

template<class CoordType>
void SpatialDescriptorDistanceToCentroid<CoordType>::eval(
  const Vertices<CoordType>& vertices,
  Vector<CoordType>& xvalues,
  Vector<CoordType>& yvalues)
{
  ENTER( "void SpatialDescriptorDistanceToBorder<CoordType>::eval(...)" );
    PRINT("here");
  const int numVertices = vertices.getSize();
  Vector<CoordType> triMeshCentroid(3);
  triMeshCentroid = _triMesh->cog();

  xvalues.setSize( numVertices );
  PRINT("here");
  for (int i = 0; i < numVertices; ++i)
  {
    //testVertex = vertices[i];
    xvalues[i] = vertices[i].distance( triMeshCentroid );
  }
  xvalues.sort();

  CDFTools<CoordType> cdftools;
  yvalues = cdftools.cdf( xvalues );


  LEAVE();
}

template class SpatialDescriptorDistanceToCentroid<float>;
template class SpatialDescriptorDistanceToCentroid<double>;
template class SpatialDescriptorDistanceToCentroid<long double>;
