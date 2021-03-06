/*!
 * \class  SpatialModelMaximalRepulsion3D
 * \author Javier Arpón (ja), INRA
 * \author Philippe Andrey (pa), INRA
 * \date   XXXX.XX.XX - creation (ja)
 * \date   2015.10.12 - integration (pa)
 * \brief  3D hardcore point process with maximal repulsion
****************************************************************/

#include "spatialmodelmaximalrepulsion3d.h"

#include <convergencetest.h>
#include <dataset.h>
#include <exception.h>
#include <stringtools.h>
#include <cmath>

#define TRACE
#include <trace.h>

/*! Constructor.
****************************************************************/
template<class CoordType>
SpatialModelMaximalRepulsion3D<CoordType>::SpatialModelMaximalRepulsion3D() : SpatialModelHardcoreDistance3D<CoordType>()
{
  _numMonteCarloCycles = 1000;
}

/*! Sets the maximum number of Monte Carlo cycles.
 *
 * The default value is 1000.
****************************************************************/
template<class CoordType>
void SpatialModelMaximalRepulsion3D<CoordType>::setNumMonteCarloCycles(const int numMonteCarloCycles)
{
  _numMonteCarloCycles = numMonteCarloCycles;
}

/*! Returns the maximum number of Monte Carlo cycles.
****************************************************************/
template<class CoordType>
int SpatialModelMaximalRepulsion3D<CoordType>::getNumMonteCarloCycles() const
{
  return _numMonteCarloCycles;
}

/*! Returns the energy of the specified configuration.
 *
 * The default implementation returns the inverse sum of distances
 * to nearest neighbours.
 *
 * This method can be reimplemented in subclasses for alternative
 * energy functions.
****************************************************************/
template<class CoordType>
CoordType SpatialModelMaximalRepulsion3D<CoordType>::energy(const Vertices<CoordType>& vertices) const
{
  Vector<CoordType> minDistances = vertices.squareNearestNeighborDistances();
  minDistances.apply( sqrt );
  return 1.0 / minDistances.mean();
}

#if 0
/*! Gets the current energy of the system
****************************************************************/
template<class CoordType>
CoordType SpatialModelMaximalRepulsion3D<CoordType>::getEnergy1(const Vertices<CoordType>& vertices)
{
  const int numPoints = vertices.getNumVertices();
  CoordType s = 0.0;

  for (int i = 0; i < numPoints; ++i)
    for (int j = i+1; j < numPoints; ++j)
       s += vertices[i].distance(vertices[j]);

  EVAL( s );
  EVAL( 1/s );
  return 1.0/s;
  //return numPoints/distBetweenComp;
}

/*! Gets the current energy of the system
****************************************************************/
template<class CoordType>
CoordType SpatialModelMaximalRepulsion3D<CoordType>::getEnergy3(const Vertices<CoordType>& vertices)
{
  CoordType energy;
  int numPoints = vertices.getNumVertices();
  T beta = _hardcoreDistances.mean()/(_nucleusVolume/numPoints);

  for ( int i = 0; i < numPoints; ++i )
    for ( int j = 0; j < numPoints; ++j )
      if ( j!= i)
        energy += exp(-(pow(vertices[i].distance(vertices[j]),2))/(2*pow(beta,2)));
  return energy;
}

/*! Gets the current energy of the system
****************************************************************/
template<class CoordType>
CoordType SpatialModelMaximalRepulsion3D<CoordType>::getEnergy4(const Vertices<CoordType>& vertices)
{
  T partialEnergy;
  int numPoints = vertices.getNumVertices();
  T beta = _hardcoreDistances.mean()/(_nucleusVolume/numPoints);

  for ( int i = 0; i < numPoints; ++i )
    for ( int j = 0; j < numPoints; ++j )
      if ( j!= i)
        partialEnergy += pow(vertices[i].distance(vertices[j]),2);

  return exp(-(partialEnergy/(2*pow(beta,2))));

}
#endif

/*! Generates a sample according to this model.
****************************************************************/
template<class CoordType>
Vertices<CoordType> SpatialModelMaximalRepulsion3D<CoordType>::drawSample(const int numPoints)
{
  ENTER("Vertices<CoordType> SpatialModelMaximalRepulsion3D<CoordType>::drawSample(const int)");

  const int outerSteps = getNumMonteCarloCycles();
  const int innerSteps = numPoints;
  const CoordType beta = computeBeta();
  const CoordType maxRadius = this->getTriMesh().equivalentRadius() / 50.0;
  Vertices<CoordType> vertices = SpatialModelHardcoreDistance3D<CoordType>::drawSample( numPoints );
  RandomGenerator& randomGenerator = this->getRandomGenerator();
  CoordType currentEnergy, newEnergy, deltaEnergy;
  Vector<CoordType> vertex;
  bool acceptTransition;
  int i, j, v;

  ConvergenceTest<CoordType> convergenceTest;  
  CoordType sumEnergy;
  bool converged = false;
  _energyProfile.setZeros( outerSteps );
  convergenceTest.setData( _energyProfile );
  convergenceTest.setRange( 100 );

  for (i = 0; !converged && i < outerSteps; ++i)
  {
#if 0
    if ( (i%20)==0 )
    {
      string filename = "max-repulsion-vertices-" + StringTools::toString(i,4,'0') + ".vx";
      vertices.save( filename, true );
    }
#endif
    for (j = 0, sumEnergy = 0.0; j < innerSteps; ++j)
    {
      v = randomGenerator.uniformL( numPoints );
      vertex = vertices[v];
      if ( j == 0 ) currentEnergy = energy( vertices );
      moveVertex( vertices, v, maxRadius );
      newEnergy = energy( vertices );
      deltaEnergy = newEnergy - currentEnergy;
      acceptTransition = deltaEnergy <= .0 || randomGenerator.uniformLF() < exp( -beta*deltaEnergy );

      if ( acceptTransition )
        currentEnergy = newEnergy;
      else vertices[v] = vertex;
      sumEnergy += currentEnergy;
    }

    _energyProfile[i] = sumEnergy / innerSteps;
    converged = convergenceTest.isPositive( i );
  }

  if ( !converged )
  {
    Exception exception;
    exception.setWhere( "Vertices<CoordType> SpatialModelMaximalRepulsion3D<CoordType>::drawSample(const int)" );
    exception.setWhat( "Convergence not reached after " + StringTools::toString(outerSteps) + " iterations" );
    throw exception;
  }

  LEAVE();

  return vertices;
}

/*! Returns the energy profile of the current or last call to
 * \c drawSample.
 *
 * The energy profile gives the energy value after each Monte Carlo step.
****************************************************************/
template<class CoordType>
const Vector<CoordType>& SpatialModelMaximalRepulsion3D<CoordType>::getEnergyProfile() const
{
  return _energyProfile;
}

/*! Randomly moves vertex \c v in \c vertices around its current
 * position within a cube of half size \c radius.
 *
 * The new position is constrained to lie within the associated
 * triMesh and to respect the hardcore distances.
****************************************************************/
template<class CoordType>
void SpatialModelMaximalRepulsion3D<CoordType>::moveVertex(
  Vertices<CoordType>& vertices,
  const int v,
  const CoordType radius)
{
  RandomGenerator& randomGenerator = this->getRandomGenerator();
  const TriMeshQuery<CoordType>& triMeshQuery = this->getTriMeshQuery();
  Vector<CoordType> vertex;

  do
  {
    EVAL( triMeshQuery.contains(vertices[v]) );
    EVAL( checkHardcoreDistances(vertices[v],v,vertices) );
    EVAL (radius);
    vertex = vertices[v];
    vertex[X] += randomGenerator.uniformLF( -radius, radius );
    vertex[Y] += randomGenerator.uniformLF( -radius, radius );
    vertex[Z] += randomGenerator.uniformLF( -radius, radius );
  } while ( !triMeshQuery.contains(vertex) || !checkHardcoreDistances(vertex,v,vertices) );

  vertices[v] = vertex;
}

/*! Checks if the hardcore distances between a new (moved) vertex and the others are respected.
****************************************************************/
template<class CoordType>
bool SpatialModelMaximalRepulsion3D<CoordType>::checkHardcoreDistances(
  const Vector<CoordType>& vertex,
  const int v,
  const Vertices<CoordType>& vertices) const
{
  const Vector<CoordType>& hardcoreDistances = this->getHardcoreDistances();
  Vector<CoordType> triMeshVertex;
  if ( this->getTriMeshQuery().closestPoint(vertex,triMeshVertex) < hardcoreDistances[v] )
    return false;

  for (int i = 0; i < vertices.getNumVertices(); ++i)
    if ( i != v && vertices[i].distance(vertex) < hardcoreDistances[i] + hardcoreDistances[v] )
      return false;

  return true;
}

/*! Compute the value of beta for the present context (mesh + objects).
 *
 * Beta is set to the value that leads to 5% acceptance of a
 * transition between a completely random configuration and the
 * same configuration with increased energy after one vertex has
 * been moved.
 *
 * The returned value is an estimation of beta as defined above.
 * The estimation is based on 100 random configurations and their
 * local perturbations.
****************************************************************/
template<class CoordType>
CoordType SpatialModelMaximalRepulsion3D<CoordType>::computeBeta()
{
  const int n = 100;
  const int numPoints = this->getHardcoreDistances().getSize();
  const CoordType maxRadius = this->getTriMesh().equivalentRadius() / 50.0;
  RandomGenerator& randomGenerator = this->getRandomGenerator();
  Vertices<CoordType> vertices1;
  Vertices<CoordType> vertices2;
  CoordType sumDelta = 0.0;

  for (int i = 0; i < n; ++i)
  {
    vertices1 = SpatialModelHardcoreDistance3D<CoordType>::drawSample( numPoints );
    vertices2 = vertices1;
    vertices2.detach();
    int v = randomGenerator.uniformL( numPoints );
    moveVertex( vertices2, v, maxRadius );
    sumDelta += fabs( energy(vertices1) - energy(vertices2) );
  }

  return log( 20.0 ) / ( sumDelta/n );
}

template class SpatialModelMaximalRepulsion3D<float>;
template class SpatialModelMaximalRepulsion3D<double>;
template class SpatialModelMaximalRepulsion3D<long double>;
