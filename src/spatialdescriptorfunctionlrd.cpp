/*!
 * \class  SpatialDescriptorFunctionLongRangeDistance
 * \author Javier Arpon (ja), INRA
 * \date   2016.02.10 - creation (ja)
 * \brief  Longe-Range-Distance-(LRD)-function statistics for spatial point processes
 * It studies the interacion between objects at long range distances taking as distance threshold
 * the maximum distance between neighbors in the observed pattern
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
void SpatialDescriptorFunctionLRD<CoordType>::setDistanceThreshold(const CoordType& distanceThreshold )
{
  _distanceThreshold = distanceThreshold;
}

template<class CoordType>
const CoordType& SpatialDescriptorFunctionLRD<CoordType>::getDistanceThreshold() const
{
  return _distanceThreshold;
}

template<class CoordType>
void SpatialDescriptorFunctionLRD<CoordType>::eval(
  const Vertices<CoordType>& vertices,
  Vector<CoordType>& x,
  Vector<CoordType>& y)
{
  const int numVertices = vertices.getSize();
  int i, j = 0;
  Vector<CoordType> temp;

  x.setSize( 0 );
  temp.setSize( 1 );

  for (i = 0; i < numVertices; ++i)
    for (j = i+1; j < numVertices; ++j)
    {
      temp[0] = vertices[i].distance( vertices[j] );
      if ( temp[0] > _distanceThreshold )
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
