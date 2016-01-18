/*!
 * \class  SpatialDescriptorFunctionZ
 * \author Javier Arpon (ja), INRA
 * \date   2016.01.18 - creation (ja)
 * \brief  Z-function statistics for spatial point processes
 * It studies the interacion between objects at long range distances
****************************************************************/

#include "spatialdescriptorfunctionz.h"

#include <cdftools.h>

#include <cmath>

template<class CoordType>
void SpatialDescriptorFunctionZ<CoordType>::eval(
  const Vertices<CoordType>& vertices,
  Vector<CoordType>& x,
  Vector<CoordType>& y)
{
  const int numVertices = vertices.getSize();
  int i, j, k = 0;
  Vector<CoordType> gFunction;
  CoordType threshold, temp;

  gFunction = vertices.squareNearestNeighborDistances();
  gFunction.apply( sqrt );
  threshold = gFunction.max();

  x.setSize( (numVertices*(numVertices-1))/2 );
  y.setSize( (numVertices*(numVertices-1))/2 );

  for (i = 0; i < numVertices; ++i)
    for (j = i+1; j < numVertices; ++j)
    {
      temp = vertices[i].distance( vertices[j] );
      if ( temp >= threshold )
        x[k++] = temp;
    }

  x.sort();

  CDFTools<CoordType> cdfTools;
  y = cdfTools.cdf( x );
}

template class SpatialDescriptorFunctionZ<float>;
template class SpatialDescriptorFunctionZ<double>;
template class SpatialDescriptorFunctionZ<long double>;
