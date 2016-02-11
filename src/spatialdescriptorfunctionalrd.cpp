/*!
 * \class  SpatialDescriptorFunctionAssymetricLongRangeDistance
 * \author Javier Arpon (ja), INRA
 * \date   2016.02.10 - creation (ja)
 * \brief  Assymetric-Longe-Range-Distance-(ALRD)-function statistics for spatial point processes
 * It studies the interacion between objects at long range distances taking as distance threshold
 * the maximum distance between neighbors in each pattern
****************************************************************/

#include "spatialdescriptorfunctionalrd.h"

#include <cdftools.h>

#include <cmath>

#define TRACE
#include <trace.h>

template<class CoordType>
SpatialDescriptorFunctionALRD<CoordType>::SpatialDescriptorFunctionALRD() : SpatialDescriptor<CoordType>()
{
}

template<class CoordType>
void SpatialDescriptorFunctionALRD<CoordType>::eval(
  const Vertices<CoordType>& vertices,
  Vector<CoordType>& x,
  Vector<CoordType>& y)
{
  const int numVertices = vertices.getSize();
  int i, j = 0;
  Vector<CoordType> temp, gFunction;

  x.setSize( 0 );
  temp.setSize( 1 );

  CoordType distanceThreshold;
  gFunction = vertices.squareNearestNeighborDistances( );
  gFunction.apply( sqrt );
  distanceThreshold = gFunction.max();

  for (i = 0; i < numVertices; ++i)
    for (j = i+1; j < numVertices; ++j)
    {
      temp[0] = vertices[i].distance( vertices[j] );
      if ( temp[0] > distanceThreshold )
        x.append( temp );
    }

  x.sort();
  y.setSize( x.getSize() );

  CDFTools<CoordType> cdfTools;
  y = cdfTools.cdf( x );
}

template class SpatialDescriptorFunctionALRD<float>;
template class SpatialDescriptorFunctionALRD<double>;
template class SpatialDescriptorFunctionALRD<long double>;
