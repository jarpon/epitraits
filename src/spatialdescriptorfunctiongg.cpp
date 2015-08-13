/*!
 * \class  SpatialDescriptorFunctionG - for more than 1 kind of object
 * \author Javier Arpon (ja), INRA
 * \date   2015.08.10 - creation (ja)
 * \brief  G'-function statistics for spatial point processes
****************************************************************/

#include "spatialdescriptorfunctiongg.h"
#include <programerror.h>
#include <cmath>

#include <cdftools.h>

template<class CoordType>
SpatialDescriptorFunctionGG<CoordType>::SpatialDescriptorFunctionGG() : SpatialDescriptor<CoordType>()
{
}

template<class CoordType>
void SpatialDescriptorFunctionGG<CoordType>::setVertices(
  const Vertices<CoordType>& vertices)
{
    _allVertices = vertices;
}

template<class CoordType>
void SpatialDescriptorFunctionGG<CoordType>::setOrder(
  const Vector<CoordType>& verticesOrder)
{
    _allVerticesOrder = verticesOrder;
}

//template<class CoordType>
//void SpatialDescriptorFunctionGG<CoordType>::setOrder(
//  const Vector<string>& verticesOrder)
//{
//    _allVerticesLabelOrder = verticesOrder;
//}


template<class CoordType>
void SpatialDescriptorFunctionGG<CoordType>::setVertices(
  const Vertices<CoordType>& vertices1, const Vertices<CoordType>& vertices2)
{
    _verticesKind1 = vertices1;
    _verticesKind2 = vertices2;
}

template<class CoordType>
void SpatialDescriptorFunctionGG<CoordType>::setVerticesKind1(
  const Vertices<CoordType>& vertices)
{
    _verticesKind1 = vertices;
}


template<class CoordType>
void SpatialDescriptorFunctionGG<CoordType>::setVerticesKind2(
  const Vertices<CoordType>& vertices)
{
    _verticesKind2 = vertices;
}

//template<class CoordType>
//Vector<Coordtype> SpatialDescriptorFunctionGG<CoordType>::getClosestNeighbors() const
//{
////  const int n = getNumVertices();

////  if ( n == 0 )
////  {
////    Vector<T> vector;
////    return vector;
////  }

////  if ( n == 1 )
////  {
////    Vector<T> vector( 1 );
////    vector[0] = 0;
////    return vector;
////  }

////  ConstVertexIterator<T> vi( *this );
////  ConstVertexIterator<T> vj( *this );
////  Vector<T> squareDistances( n-1 );
////  Vector<T> snnd( n );
////  int i, j, k;

////  for (i = 0; i < n; ++i, ++vi)
////  {
////    const Vector<T>& rvi = *vi.current();
////    vj.rewind();
////    for (j = 0, k = 0; j < n; ++j, ++vj)
////      if ( i != j )
////        squareDistances[k++] = rvi.sdistance( *vj.current() );
////    snnd[i] = squareDistances.min();
////  }

////  return snnd;
//}


template<class CoordType>
void SpatialDescriptorFunctionGG<CoordType>::eval(
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

//  x = getClosestNeighbors();
  x.apply( sqrt );
  x.sort();


  CDFTools<CoordType> cdfTools;
  y = cdfTools.cdf( x );
}

template class SpatialDescriptorFunctionGG<float>;
template class SpatialDescriptorFunctionGG<double>;
template class SpatialDescriptorFunctionGG<long double>;

