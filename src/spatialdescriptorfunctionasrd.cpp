/*!
 * \class  SpatialDescriptorFunctionAssymetricShortRangeDistance
 * \author Javier Arpon (ja), INRA
 * \date   2016.02.18 - creation (ja)
 * \brief  Assymetric-Longe-Range-Distance-(ASRD)-function statistics for spatial point processes
 * It studies the interacion between objects at short range distances taking as distance threshold
 * the maximum distance between neighbors in each pattern
****************************************************************/

#include "spatialdescriptorfunctionasrd.h"

#include <cdftools.h>

#include <cmath>

#define TRACE
#include <trace.h>

template<class CoordType>
SpatialDescriptorFunctionASRD<CoordType>::SpatialDescriptorFunctionASRD() : SpatialDescriptor<CoordType>()
{
}

template<class CoordType>
void SpatialDescriptorFunctionASRD<CoordType>::eval(
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

template class SpatialDescriptorFunctionASRD<float>;
template class SpatialDescriptorFunctionASRD<double>;
template class SpatialDescriptorFunctionASRD<long double>;

