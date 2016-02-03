/*!
 * \class  SpatialDescriptorFunctionLongRangeDistance
 * \author Javier Arpon (ja), INRA
 * \date   2016.01.18 - creation (ja)
 * \brief  LRD-function statistics for spatial point processes
 * It studies the interacion between objects at long range distances
****************************************************************/

#include "spatialdescriptorfunctionlrd.h"

#include <cdftools.h>

#include <cmath>

#define TRACE
#include <trace.h>

template<class CoordType>
SpatialDescriptorFunctionLRD<CoordType>::SpatialDescriptorFunctionLRD() : SpatialDescriptor<CoordType>()
{
}

template<class CoordType>
void SpatialDescriptorFunctionLRD<CoordType>::eval(
  const Vertices<CoordType>& vertices,
  Vector<CoordType>& x,
  Vector<CoordType>& y)
{
  const int numVertices = vertices.getSize();
  int i, j, k = 0;
  Vector<CoordType> gFunction, temp;
  CoordType threshold;
  gFunction = vertices.squareNearestNeighborDistances();
  gFunction.apply( sqrt );
  threshold = gFunction.max();

  x.setSize( 0 );
  temp.setSize( 1 );

  for (i = 0; i < numVertices; ++i)
    for (j = i+1; j < numVertices; ++j)
    {
      temp[0] = vertices[i].distance( vertices[j] );
      if ( temp[0] > threshold )
        x.append( temp );
    }

  x.sort();
  y.setSize( x.getSize() );

  CDFTools<CoordType> cdfTools;
  y = cdfTools.cdf( x );
}

template class SpatialDescriptorFunctionLRD<float>;
template class SpatialDescriptorFunctionLRD<double>;
template class SpatialDescriptorFunctionLRD<long double>;
