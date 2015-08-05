/*!
 * \class  SpatialModelBoundaryInteraction
 * \author Philippe Andrey (pa), INRA
 * \date   2012.01.16 - creation (pa)
 * \brief  Boundary interaction model
 * \details
 * This class implements a 2D spatial model with an interaction
 * to the boundary.
 * A margin is defined at the periphery by an arbitrary proportion.
 * With a specified probability, each point is taken in this margin.
 * Conditioned on the domain (peripheral or central), each point
 * is uniformly distributed.
 *
 * \code
 * SpatialModelBoundaryInteraction<float> spatialModel;
 * spatialModel.setBoundary( yourBoundary ); // defines the domain
 * spatialModel.setMargin( 0.1 ); // margin is 10% of domain size
 * spatialModel.setMarginProb( 0.5 ); // prob of being in the margin
 * spatialModel.initialize();
 * yourVertices = spatialModel.drawVertices( 10 ); // simulates the model
 * \endcode
****************************************************************/

#include "spatialmodelboundaryinteraction.h"

//#define TRACE
#include <trace.h>

template<class CoordType>
SpatialModelBoundaryInteraction<CoordType>::SpatialModelBoundaryInteraction()
  : SpatialModel2D<CoordType,CoordType>()
{
  _margin = 0.0;
  _marginProb = 0.0;
}

template<class CoordType>
SpatialModelBoundaryInteraction<CoordType>::~SpatialModelBoundaryInteraction()
{
}

template<class CoordType>
void SpatialModelBoundaryInteraction<CoordType>::setMargin(const CoordType margin)
{
  _margin = margin;
}

template<class CoordType>
CoordType SpatialModelBoundaryInteraction<CoordType>::getMargin() const
{
  return _margin;
}

template<class CoordType>
void SpatialModelBoundaryInteraction<CoordType>::setMarginProb(const float marginProb)
{
  _marginProb = marginProb;
}

template<class CoordType>
float SpatialModelBoundaryInteraction<CoordType>::getMarginProb() const
{
  return _marginProb;
}

template<class CoordType>
void SpatialModelBoundaryInteraction<CoordType>::initialize()
{
  const Curve<CoordType>& boundary = this->getBoundary();
  _innerBoundary = boundary;
  _innerBoundary.detach();
  _innerBoundary.center();
  _innerBoundary.scale( 1.0-2.0*_margin );
  _innerBoundary.translate( boundary.cog() );
  _innerBoundingBox = _innerBoundary.boundingBox();
}

template<class CoordType>
const Curve<CoordType>& SpatialModelBoundaryInteraction<CoordType>::getInnerBoundary() const
{
  return _innerBoundary;
}

template<class CoordType>
Vertices<CoordType> SpatialModelBoundaryInteraction<CoordType>::drawSample(const int numVertices)
{
  Vertices<CoordType> vertices( 2, numVertices );
  RandomGenerator& randomGenerator = this->getRandomGenerator();
  const BoundingBox<CoordType>& outerBoundingBox = this->getBoundingBox();
  const Curve<CoordType>& outerBoundary = this->getBoundary();

  for (int v = 0; v < numVertices; ++v)
  {
    if ( randomGenerator.uniformLF() < _marginProb )
    {
      // draw a point at the periphery
      do {
        vertices[v][X] = randomGenerator.uniformLF( outerBoundingBox.min(X), outerBoundingBox.max(X) );
        vertices[v][Y] = randomGenerator.uniformLF( outerBoundingBox.min(Y), outerBoundingBox.max(Y) );
      } while ( !outerBoundary.contains(vertices[v]) || _innerBoundary.contains(vertices[v]) );
    }
    else
    {
      // draw a point in the central domain
      do {
        vertices[v][X] = randomGenerator.uniformLF( _innerBoundingBox.min(X), _innerBoundingBox.max(X) );
        vertices[v][Y] = randomGenerator.uniformLF( _innerBoundingBox.min(Y), _innerBoundingBox.max(Y) );
      } while ( _innerBoundary.contains(vertices[v],X,Y) == false );
    }
  }

  return vertices;
}

template class SpatialModelBoundaryInteraction<float>;
template class SpatialModelBoundaryInteraction<double>;
template class SpatialModelBoundaryInteraction<long double>;
