/*!
 * \class  SpatialModel2
 * \author Philippe Andrey (pa), INRA
 * \date   2012.01.16 - creation (pa)
 * \brief  Base class for spatial point process models
****************************************************************/

#include "spatialmodel2.h"

#include <programerror.h>

template<class CoordType,class PixelType>
SpatialModel2<CoordType,PixelType>::SpatialModel2()
{
  _randomGenerator = 0;
  _pixelSize.setOnes( 2 );
  _mask = 0;
}

template<class CoordType,class PixelType>
SpatialModel2<CoordType,PixelType>::~SpatialModel2()
{
}

/*! Sets the random number generator to use.
****************************************************************/
template <class CoordType,class PixelType>
void SpatialModel2<CoordType,PixelType>::setRandomGenerator(RandomGenerator& randomGenerator)
{
  _randomGenerator = &randomGenerator;
}

/*! Returns the random number generator.
****************************************************************/
template <class CoordType,class PixelType>
RandomGenerator& SpatialModel2<CoordType,PixelType>::getRandomGenerator()
{
  return *_randomGenerator;
}

/*! Sets the physical size of pixels.
****************************************************************/
template<class CoordType,class PixelType>
void SpatialModel2<CoordType,PixelType>::setPixelSize(const Vector<CoordType>& pixelSize)
{
  _pixelSize = pixelSize;
}

/*! Returns the physical size of pixels.
****************************************************************/
template<class CoordType,class PixelType>
const Vector<CoordType>& SpatialModel2<CoordType,PixelType>::getPixelSize() const
{
  return _pixelSize;
}

/*! Sets the mask defining the domain.
****************************************************************/
template<class CoordType,class PixelType>
void SpatialModel2<CoordType,PixelType>::setMask(const PixelMatrix<PixelType>& mask)
{
  _mask = &mask;
}

/*! Returns the mask defining the domain.
****************************************************************/
template<class CoordType,class PixelType>
const PixelMatrix<PixelType>& SpatialModel2<CoordType,PixelType>::getMask() const
{
  if ( _mask == 0 )
  {
    ProgramError programError;
    programError.setWhere( "const PixelMatrix<PixelType>& SpatialModel2<CoordType,PixelType>::getMask() const" );
    programError.setWhat( "No mask is set" );
    throw programError;
  }

  return *_mask;
}

///*! Generates \c numSamples of \c numVertices each.
//****************************************************************/
//template<class CoordType,class PixelType>
//ShapeSet<CoordType> SpatialModel2<CoordType,PixelType>::drawSamples(
//  const int numSamples,
//  const int numVertices)
//{
//  ShapeSet<CoordType> shapeSet;

//  for (int i = 0; i < numSamples; ++i)
//  {
//    shapeSet.addShape( new Vertices<CoordType>(drawSample(numVertices)) );
//  }

//  return shapeSet;
//}

///*! Generates \c numSamples of \c numVertices each.
//****************************************************************/
//template<class CoordType,class PixelType>
//ShapeSet<CoordType> SpatialModel2<CoordType,PixelType>::drawSamples(
//  const int numSamples,
//  const int numVerticesDist1,
//  const int numVerticesDist2)
//{
//  ShapeSet<CoordType> shapeSet;

//  for (int i = 0; i < numSamples; ++i)
//  {
//    shapeSet.addShape( new Vertices<CoordType>(drawSample(numVerticesDist1,numVerticesDist2)) );
//  }

//  return shapeSet;
//}

template class SpatialModel2<float,float>;
template class SpatialModel2<double,double>;
template class SpatialModel2<long double,long double>;

template class SpatialModel2<double,float>;
template class SpatialModel2<long double,float>;
