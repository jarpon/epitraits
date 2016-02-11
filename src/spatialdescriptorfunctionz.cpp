/*!
 * \class  SpatialDescriptorFunctionZ
 * \author Javier Arpon (ja), INRA
 * \date   2016.01.18 - creation (ja)
 * \brief  Z-function statistics for spatial point processes
 * It studies the interacion between each objects and its furthest one
****************************************************************/

#include "spatialdescriptorfunctionz.h"

#include <cdftools.h>

#include <cmath>

#define TRACE
#include <trace.h>

template<class CoordType>
SpatialDescriptorFunctionZ<CoordType>::SpatialDescriptorFunctionZ() : SpatialDescriptor<CoordType>()
{
}

template<class CoordType>
void SpatialDescriptorFunctionZ<CoordType>::eval(
  const Vertices<CoordType>& vertices,
  Vector<CoordType>& x,
  Vector<CoordType>& y)
{
  const int numVertices = vertices.getSize();
  int i, j = 0;
  x.setSize( numVertices );
  x.setZeros();
  CoordType temp;

  for (i = 0; i < numVertices; ++i)
    for (j = 0; j < numVertices ; ++j)
      if ( i != j )
      {
        temp = vertices[i].distance( vertices[j] );
        if ( temp > x[i] )
          x[i] = temp;
      }

  x.sort();
  y.setSize( x.getSize() );

  CDFTools<CoordType> cdfTools;
  y = cdfTools.cdf( x );
}

template class SpatialDescriptorFunctionZ<float>;
template class SpatialDescriptorFunctionZ<double>;
template class SpatialDescriptorFunctionZ<long double>;
