/*!
 * \class  SpatialModelMaximalRepulsion3D2
 * \author Javier Arpón (ja), INRA
 * \author Philippe Andrey (pa), INRA
 * \date   XXXX.XX.XX - creation (ja)
 * \date   2015.10.12 - integration (pa)
 * \brief  3D hardcore point process with maximal repulsion
****************************************************************/

#include "spatialmodelmaximalrepulsion3d2.h"
#include <convergencetest.h>
#include <dataset.h>
#include <exception.h>
#include <stringtools.h>

#include <cmath>
#include <cstdlib>

#define TRACE
#include <trace.h>

/*! Constructor.
****************************************************************/
template<class CoordType>
SpatialModelMaximalRepulsion3D2<CoordType>::SpatialModelMaximalRepulsion3D2() : SpatialModelHardcoreDistance3D<CoordType>()
{
  _numMonteCarloCycles = 1000;
  _beta = -1.0;
}

/*! Sets the maximum number of Monte Carlo cycles.
 *
 * The default value is 1000.
****************************************************************/
template<class CoordType>
void SpatialModelMaximalRepulsion3D2<CoordType>::setNumMonteCarloCycles(const int numMonteCarloCycles)
{
  _numMonteCarloCycles = numMonteCarloCycles;
}

/*! Returns the maximum number of Monte Carlo cycles.
****************************************************************/
template<class CoordType>
int SpatialModelMaximalRepulsion3D2<CoordType>::getNumMonteCarloCycles() const
{
  return _numMonteCarloCycles;
}

/*! Returns the current value of the beta parameter.
 *
 * \sa initializeBeta
****************************************************************/
template<class CoordType>
CoordType SpatialModelMaximalRepulsion3D2<CoordType>::getBeta() const
{
  return _beta;
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
 *
 * \sa getBeta
****************************************************************/
template<class CoordType>
void SpatialModelMaximalRepulsion3D2<CoordType>::initializeBeta(const int numObjects)
{
  const int n = 100;
  const CoordType maxRadius = this->getTriMesh().equivalentRadius() / 50.0;
  RandomGenerator& randomGenerator = this->getRandomGenerator();
  Vertices<CoordType> vertices1;
  Vertices<CoordType> vertices2;
  CoordType sumDelta = 0.0;

  for (int i = 0; i < n; ++i)
  {
    vertices1 = SpatialModelHardcoreDistance3D<CoordType>::drawSample( numObjects );
    vertices2 = vertices1;
    vertices2.detach();
    int v = randomGenerator.uniformL( numObjects );
    moveVertex( vertices2, v, maxRadius );
    sumDelta += fabs( energy(vertices1) - energy(vertices2) );
  }

  _beta = log( 20.0 ) / ( sumDelta/n );
}



/*! Returns the energy of the specified configuration.
 *
 * The default implementation returns the inverse sum of distances
 * to furthest neighbours.
****************************************************************/
template<class CoordType>
CoordType SpatialModelMaximalRepulsion3D2<CoordType>::energy(const Vertices<CoordType>& vertices) const
{
  const int numVertices = vertices.getNumVertices();
  CoordType temp;
  Vector<CoordType> x (numVertices);

  for ( int i = 0; i < numVertices; ++i)
  {
    x[i] = FLT_MAX;
    for ( int j = 0; j < numVertices ; ++j)
      if ( i != j )
      {
        temp = vertices[i].distance( vertices[j] );
        if ( temp < x[i] )
          x[i] = 1.0/temp;
      }
  }
//  EVAL(x.sum());
  return x.sum();
}


#if 0
/*! Returns the energy of the specified configuration.
 *
 * The default implementation returns the inverse sum of distances
 * to all chromocenters.
****************************************************************/
template<class CoordType>
CoordType SpatialModelMaximalRepulsion3D2<CoordType>::energy(const Vertices<CoordType>& vertices) const
{
  const int numVertices = vertices.getNumVertices();
  Vector<CoordType> sumInterDistances ( numVertices );

  for (int i = 0; i < numVertices; ++i)
    //for (int j = i+1; j < numVertices; ++j)
    for (int j = 0; j < numVertices, j != i; ++j)
//      sumInterDistances[i] += 1.0/vertices[i].distance(vertices[j]);
      sumInterDistances[i] += vertices[i].distance(vertices[j]);

//  return sumInterDistances;
  return 1.0/sumInterDistances.mean();
}


/*! Returns the energy of the specified configuration.
 *
 * The default implementation returns the inverse sum of distances
 * to furthest neighbours.
****************************************************************/
template<class CoordType>
CoordType SpatialModelMaximalRepulsion3D2<CoordType>::energy(const Vertices<CoordType>& vertices) const
{
  const int numVertices = vertices.getNumVertices();
  CoordType temp;
  Vector<CoordType> x (numVertices);

  for ( int i = 0; i < numVertices; ++i)
    for ( int j = 0; j < numVertices ; ++j)
      if ( i != j )
      {
        temp = vertices[i].distance( vertices[j] );
        if ( temp < x[i] )
          x[i] = 1.0/temp;
      }
  return x.mean();
}

/*! Returns the energy of the specified configuration.
 *
 * The default implementation returns the inverse sum of distances
 * to furthest neighbours.
****************************************************************/
template<class CoordType>
CoordType SpatialModelMaximalRepulsion3D2<CoordType>::energy(const Vertices<CoordType>& vertices) const
{
  const int numVertices = vertices.getNumVertices();
  CoordType temp;
  Vector<CoordType> x (numVertices);

  for ( int i = 0; i < numVertices; ++i)
    for ( int j = 0; j < numVertices ; ++j)
      if ( i != j )
      {
        temp = vertices[i].distance( vertices[j] );
        if ( temp > x[i] )
          x[i] = temp;
      }
  return 1.0 / x.mean();
}


/*! Returns the energy of the specified configuration.
 *
 * The default implementation returns the inverse sum of distances
 * to all chromocenters.
****************************************************************/
template<class CoordType>
CoordType SpatialModelMaximalRepulsion3D2<CoordType>::energy(const Vertices<CoordType>& vertices) const
{
  const int numChromocenters = vertices.getNumVertices();
  CoordType sumInterDistances = 0.0;
  for (int i = 0; i < numChromocenters; ++i)
    for (int j = i+1; j < numChromocenters; ++j)
       sumInterDistances += vertices[i].distance(vertices[j]);

  return 1.0/sumInterDistances;
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
CoordType SpatialModelMaximalRepulsion3D2<CoordType>::energy(const Vertices<CoordType>& vertices) const
{
  Vector<CoordType> minDistances = vertices.squareNearestNeighborDistances();
  minDistances.apply( sqrt );
  return 1.0 / minDistances.mean();
}



/*! Gets the current energy of the system
****************************************************************/
template<class CoordType>
CoordType SpatialModelMaximalRepulsion3D2<CoordType>::getEnergy3(const Vertices<CoordType>& vertices)
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
CoordType SpatialModelMaximalRepulsion3D2<CoordType>::getEnergy4(const Vertices<CoordType>& vertices)
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
Vertices<CoordType> SpatialModelMaximalRepulsion3D2<CoordType>::drawSample(const int numPoints)
{
  ENTER("Vertices<CoordType> SpatialModelMaximalRepulsion3D2<CoordType>::drawSample(const int)");

  const int outerSteps = getNumMonteCarloCycles();
  const int innerSteps = numPoints;
  const CoordType beta = getBeta();
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
    exception.setWhere( "Vertices<CoordType> SpatialModelMaximalRepulsion3D2<CoordType>::drawSample(const int)" );
    exception.setWhat( "Convergence not reached after " + StringTools::toString(outerSteps) + " iterations" );
    throw exception;
  }

  if ( !this->valid(vertices) )
  {
    cerr << "!!!!!!!!!!! ERROR !!!!!!!!!!" << endl;
    cerr << "In Vertices<CoordType> SpatialModelMaximalRepulsion3D2<CoordType>::drawSample(const int)" << endl;
    exit( EXIT_FAILURE );
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
const Vector<CoordType>& SpatialModelMaximalRepulsion3D2<CoordType>::getEnergyProfile() const
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
void SpatialModelMaximalRepulsion3D2<CoordType>::moveVertex(
  Vertices<CoordType>& vertices,
  const int v,
  const CoordType radius)
{
  RandomGenerator& randomGenerator = this->getRandomGenerator();
  const TriMeshQuery<CoordType>& triMeshQuery = this->getTriMeshQuery();
  Vector<CoordType> vertex;

  do
  {
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
bool SpatialModelMaximalRepulsion3D2<CoordType>::checkHardcoreDistances(
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

template class SpatialModelMaximalRepulsion3D2<float>;
template class SpatialModelMaximalRepulsion3D2<double>;
template class SpatialModelMaximalRepulsion3D2<long double>;
