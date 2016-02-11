/*!
 * \class  SpatialDescriptorFunctionNotNeighbors
 * \author Javier Arpon (ja), INRA
 * \date   2016.02.10 - creation (ja)
 * \brief  Not-Neighbors-(NN)-function statistics for spatial point processes
 * It measures all the inter-distances among objects without distances between neighbors
****************************************************************/

#include "spatialdescriptorfunctionnn.h"

#include <cdftools.h>

#include <cmath>

#define TRACE
#include <trace.h>

template<class CoordType>
SpatialDescriptorFunctionNN<CoordType>::SpatialDescriptorFunctionNN() : SpatialDescriptor<CoordType>()
{
}

template<class CoordType>
void SpatialDescriptorFunctionNN<CoordType>::eval(
  const Vertices<CoordType>& vertices,
  Vector<CoordType>& x,
  Vector<CoordType>& y)
{
  const int numVertices = vertices.getSize();
  int i, j = 0;
  Vector<CoordType> temp, gFunction;

  x.setSize( 0 );
  temp.setSize( 1 );

  gFunction = vertices.squareNearestNeighborDistances( );
  gFunction.apply( sqrt );

  for (i = 0; i < numVertices; ++i)
    for (j = i+1; j < numVertices; ++j)
    {
      temp[0] = vertices[i].distance( vertices[j] );
      if ( temp[0] != gFunction[i] )
        x.append( temp );
    }

  x.sort();
  y.setSize( x.getSize() );

  CDFTools<CoordType> cdfTools;
  y = cdfTools.cdf( x );
}

template class SpatialDescriptorFunctionNN<float>;
template class SpatialDescriptorFunctionNN<double>;
template class SpatialDescriptorFunctionNN<long double>;

