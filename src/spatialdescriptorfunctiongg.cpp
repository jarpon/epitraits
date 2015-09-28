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

#define TRACE
#include <trace.h>

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
//Vector<CoordType> SpatialDescriptorFunctionGG<CoordType>::getClosestNeighbor()
//{
//  const int n1 = _verticesKind1.getNumVertices();
//  const int n2 = _verticesKind2.getNumVertices();

//  if ( ( n1 == 0 ) || ( n2 == 0 ) )
//  {
//    Vector<CoordType> vector;
//    return vector;
//  }

//  if ( ( n1 == 1 ) || ( n2 == 1 ) )
//  {
//    Vector<CoordType> vector( 1 );
//    vector[0] = 0;
//    return vector;
//  }

////  ConstVertexIterator<CoordType> vi( *this );
////  ConstVertexIterator<CoordType> vj( *this );
////  ConstVertexIterator<CoordType> vi( _verticesKind1.vertexIterator() ); //if this does not work, move it to shape.cpp and use lines above
////  ConstVertexIterator<CoordType> vj();

//  Vector<CoordType> squareDistances( n1+n2-1 );
//  Vector<CoordType> snnd( n1+n2 );
//  int i, j, k;

//  for (i = 0; i < n1; ++i, ++vi)
//  {
//    const Vector<CoordType>& rvi = *vi.current();
//    vj.rewind();
//    for (j = 0, k = 0; j < n2; ++j, ++vj)
//      //if ( i != j )
//        squareDistances[k++] = rvi.sdistance( *vj.current() );
//    snnd[i] = squareDistances.min();
//  }

//  for (i = 0; i < n2; ++i, ++vi)
//  {
//    const Vector<CoordType>& rvi = *vi.current();
//    vj.rewind();
//    for (j = 0, k = squareDistances.getSize(); j < n1; ++j, ++vj)
//      //if ( i != j )
//        squareDistances[k++] = rvi.sdistance( *vj.current() );
//    snnd[i+n1] = squareDistances.min();
//  }

//  return snnd;
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
    error.setWhat( "There is not other kind of objects set up" );
    error.setWhat( "Are you looking for the normal G-Function?" );
  }

  EVAL("inside G'");
  x = _verticesKind1.nearestNeighborDistances( _verticesKind2 );
  x.apply( sqrt );
  x.sort();


  CDFTools<CoordType> cdfTools;
  y = cdfTools.cdf( x );
}

template class SpatialDescriptorFunctionGG<float>;
template class SpatialDescriptorFunctionGG<double>;
template class SpatialDescriptorFunctionGG<long double>;

