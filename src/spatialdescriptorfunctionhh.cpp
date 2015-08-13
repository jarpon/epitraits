/*!
 * \class  SpatialDescriptorFunctionH - for more than 1 kind of object
 * \author Javier Arpon (ja), INRA
 * \date   2015.08.10 - creation (ja)
 * \brief  H'-function statistics for spatial point processes
****************************************************************/

#include "spatialdescriptorfunctionhh.h"
#include <programerror.h>

#include <cdftools.h>

template<class CoordType>
SpatialDescriptorFunctionHH<CoordType>::SpatialDescriptorFunctionHH() : SpatialDescriptor<CoordType>()
{
}

template<class CoordType>
void SpatialDescriptorFunctionHH<CoordType>::setVertices(
  const Vertices<CoordType>& vertices)
{
    _allVertices = vertices;
}

template<class CoordType>
void SpatialDescriptorFunctionHH<CoordType>::setOrder(
  const Vector<CoordType>& verticesOrder)
{
    _allVerticesOrder = verticesOrder;
}

//template<class CoordType>
//void SpatialDescriptorFunctionHH<CoordType>::setOrder(
//  const Vector<string>& verticesOrder)
//{
//    _allVerticesLabelOrder = verticesOrder;
//}


template<class CoordType>
void SpatialDescriptorFunctionHH<CoordType>::setVertices(
  const Vertices<CoordType>& vertices1, const Vertices<CoordType>& vertices2)
{
    _verticesKind1 = vertices1;
    _verticesKind2 = vertices2;
}

template<class CoordType>
void SpatialDescriptorFunctionHH<CoordType>::setVerticesKind1(
  const Vertices<CoordType>& vertices)
{
    _verticesKind1 = vertices;
}


template<class CoordType>
void SpatialDescriptorFunctionHH<CoordType>::setVerticesKind2(
  const Vertices<CoordType>& vertices)
{
    _verticesKind2 = vertices;
}


template<class CoordType>
void SpatialDescriptorFunctionHH<CoordType>::eval(
  const Vertices<CoordType>& vertices,
  Vector<CoordType>& x,
  Vector<CoordType>& y)
{
  if ( _verticesKind2.getNumVertices() == 0 )
  {
    ProgramError error;
    error.setWhat( "Error calling the program" );
    error.setWhat( "There is not other kind of objects set" );
    error.setWhat( "Are you looking for the normal H-Function?s" );
  }

  const int numVertices = vertices.getSize();
  int i, j, k = 0;

  x.setSize( (numVertices*(_verticesKind2.getSize())) );
  y.setSize( (numVertices*(_verticesKind2.getSize())) );

  for (i = 0; i < numVertices; ++i)
  {
    for (j = i+1; j < _verticesKind2.getSize(); ++j)
    {
      x[k++] = vertices[i].distance( _verticesKind2[j] );
    }
  }
  x.sort();

  CDFTools<CoordType> cdfTools;
  y = cdfTools.cdf( x );
}

template class SpatialDescriptorFunctionHH<float>;
template class SpatialDescriptorFunctionHH<double>;
template class SpatialDescriptorFunctionHH<long double>;

