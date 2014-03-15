#include "trimeshspatialmodel.h"

#include <randomgenerator.h>
#include <boundingbox.h>
#include <trimesh.h>
#include <alinevector.h>
#include <cmath>

#include <trimeshquery.h>

#include <programerror.h>

//#define TRACE
#include <trace.h>

using namespace std;



///*! Sets the random number generator to use for generating vertices.
//****************************************************************/
//template<class CoordType>
//void TriMeshSpatialModel<CoordType>::setRandomGenerator( RandomGenerator& randomGenerator )
//{
//  _randomGenerator = &randomGenerator;
//}

/*! Initializes all the parameters.
****************************************************************/
template<class CoordType>
TriMeshSpatialModel<CoordType>::TriMeshSpatialModel()
{
  _triMesh = 0;
  _numCompartments = 0;
//  _volumeRadius = 0;
//  _volumeRadiusRange.setZeros(2);
  _hardcoreDistance.setZeros(1);
  _hardcoreDistanceRange.setZeros(2);
  _distanceToBorder.setZeros(1);
  _distanceToBorderRange.setZeros(2);
}

/*! Destroys it.
****************************************************************/
template<class CoordType>
TriMeshSpatialModel<CoordType>::~TriMeshSpatialModel()
{
}

/*! Sets the number of compartments to generate
****************************************************************/
template<class CoordType>
void TriMeshSpatialModel<CoordType>::setNumberOfCompartments(const int numCompartments)
{
  _numCompartments = numCompartments;
}

/*! Sets the triMesh to work with.
****************************************************************/
template<class CoordType>
void TriMeshSpatialModel<CoordType>::setTriMesh(const TriMesh<CoordType>& triMesh)
{
  //_triMesh = &triMesh;
  _triMeshQuery.setTriMesh( triMesh );
  //_triMesh.computeNormals(); // see in the future if this is needed or need to be changed
}

/*! Initializes the class.
****************************************************************/
template<class CoordType>
void TriMeshSpatialModel<CoordType>::initialize()
{
  //computes the bounding box of the triMesh and saves it
  const BoundingBox<CoordType> boundingBox = _triMeshQuery.getTriMesh().boundingBox();
  //boundingBox.vertices().save( _outputDir + "Nucleus-BB", true );

  //_triMeshQuery.setTriMesh( _triMesh );

  xMinBB = boundingBox.min(0);
  xMaxBB = boundingBox.max(0);
  yMinBB = boundingBox.min(1);
  yMaxBB = boundingBox.max(1);
  zMinBB = boundingBox.min(2);
  zMaxBB = boundingBox.max(2);

  //shows limits of the boundingBox
  EVAL(xMinBB);
  EVAL(xMaxBB);
  EVAL(yMinBB);
  EVAL(yMaxBB);
  EVAL(zMinBB);
  EVAL(zMaxBB);

  //shows the maximum distance to the border that can be used
  float maxDistance = ( boundingBox.max(0) - boundingBox.min(0) )/2;
  for ( int i = 1; i < 3; ++i)
  {
    if ( ( boundingBox.max(i) - boundingBox.min(i) )/2 < maxDistance )
      maxDistance = ( boundingBox.max(i) - boundingBox.min(i) )/2;
  }
  EVAL(maxDistance)
}

/*! Sets a unique hardcore distance for all compartments.
****************************************************************/
template<class CoordType>
void TriMeshSpatialModel<CoordType>::setHardcoreDistance(const CoordType hardcoreDistance)
{
  _hardcoreDistance[0] = hardcoreDistance;
}

/*! Sets a vector of hardcore distances.
****************************************************************/
template<class CoordType>
void TriMeshSpatialModel<CoordType>::setHardcoreDistance(const Vector<CoordType>& hardcoreDistances)
{
  _hardcoreDistance = hardcoreDistances;
}

/*! Gets the correspondent hardcore distances of element n.
****************************************************************/
template<class CoordType>
const CoordType& TriMeshSpatialModel<CoordType>::getHardcoreDistance(const int numCompartment) const
{
  return _hardcoreDistance[numCompartment-1];
}

/*! Gets a vector of hardcore distances.
****************************************************************/
template<class CoordType>
Vector<CoordType>& TriMeshSpatialModel<CoordType>::getHardcoreDistance()
{
  return _hardcoreDistance;
}

/*! Sets a unique distance to the border for all compartments.
****************************************************************/
template<class CoordType>
void TriMeshSpatialModel<CoordType>::setDistanceToBorder(const CoordType distanceToBorder)
{
  _distanceToBorder[0] = distanceToBorder;
}

/*! Sets a vector of distances to the border.
****************************************************************/
template<class CoordType>
void TriMeshSpatialModel<CoordType>::setDistanceToBorder(const Vector<CoordType>& distancesToTheBorder)
{
  _distanceToBorder = distancesToTheBorder;
}

///*! Generates \c numSamples of \c numVertices each.
//****************************************************************/
//template<class CoordType>
//ShapeSet<CoordType> TriMeshSpatialModel<CoordType>::drawSamples(
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

/*! Diverts to the correct function depending on the conditions.
****************************************************************/
template<class CoordType>
Vertices<CoordType> TriMeshSpatialModel<CoordType>::drawSample(const int numPoints)
{
  _numCompartments = numPoints;

  if ( ( _hardcoreDistance[0] == 0 ) && ( _hardcoreDistanceRange[1] == 0 ) &&
       ( _distanceToBorder[0] == 0 ) && ( _distanceToBorderRange[1] == 0 ) )

    return randomVertices();

  else if ( ( _hardcoreDistance[0] == 0 ) && ( _hardcoreDistanceRange[1] == 0 ) &&
          ( ( _distanceToBorder[0] != 0 ) || ( _distanceToBorderRange[1] != 0 ) ) )

    return distanceToTheBorder();

  else if ( ( ( _hardcoreDistance[0] != 0 ) || ( _hardcoreDistanceRange[1] != 0 ) ) &&
              ( _distanceToBorder[0] == 0 ) && ( _distanceToBorderRange[1] == 0 ) )

    return hardcoreDistance();

  else //if ( ( ( _hardcoreDistance[0] != 0 ) || ( _hardcoreDistanceRange[1] != 0 ) ) &&
       //     ( ( _distanceToBorder[0] != 0 ) || ( _distanceToBorderRange[1] != 0 ) ) )

    return hardcoreAndToTheBorderDistances();

}

/*! Generates a random vertex into the trimesh.
****************************************************************/
template<class CoordType>
void TriMeshSpatialModel<CoordType>::drawPosition(Vector<CoordType>& position)
{
  ENTER( "void TriMeshSpatialModel<CoordType>::drawPosition(Vector<CoordType>&)" );

  RandomGenerator& randomGenerator = this->getRandomGenerator();

  do
  {
    position[X] = randomGenerator.uniformLF(xMinBB,xMaxBB);
    position[Y] = randomGenerator.uniformLF(yMinBB,yMaxBB);
    position[Z] = randomGenerator.uniformLF(zMinBB,zMaxBB);
  //} while ( _triMesh->contains(position) == false );
    } while ( _triMeshQuery.contains(position) == false );

  LEAVE();
}

/*! Generates a random vertex into the trimesh taking into account a distance to the border.
****************************************************************/
template<class CoordType>
void TriMeshSpatialModel<CoordType>::drawPositionFromBorder(Vector<CoordType>& position, CoordType distanceToBorder)
{
  //ENTER("void TriMeshSpatialModel<CoordType>::drawPositionFromBorder(Vector<CoordType,PixelType>, T)");
  Vector<CoordType> triMeshVertex, normalVector;
  int attempts = 0;
  const int maxAttempts = 10000;
  bool found = false;

  do
  {
    drawPosition( position );

    //from the random poin we find the closest point to it in the triMesh
    //the segment between them will be the normal vector
    //from this normal vector and using the distance to the border we obtain
    //the "position" vertex which will be our wanted vertex

    //_triMesh->closestPoint( position, triMeshVertex );
    _triMeshQuery.getTriMesh().closestPoint( position, triMeshVertex );

    normalVector = position - triMeshVertex;
    normalVector /= position.distance( triMeshVertex );
    normalVector *= distanceToBorder;
    position = triMeshVertex + normalVector;

    //we check if the obtained point is inside of the triMesh
    //and if it keeps the condition of being the minimum distance to the border
    //the distance to the border introduced manually before

    //if ( _triMesh->contains(position) )
    if ( _triMeshQuery.contains(position) )
    {
      //_triMesh->closestPoint( position, triMeshVertex );
      _triMeshQuery.getTriMesh().closestPoint( position, triMeshVertex );
      found = fabs( distanceToBorder - position.distance(triMeshVertex) ) < position.epsilon();
    }
  } while ( found == false && ++attempts <= maxAttempts );

  if ( found == false || attempts > maxAttempts )
  {
    Exception exception;
    exception.setWhere( "void TriMeshSpatialModel<CoordType>::drawPositionFromBorder()" );
    exception.setWhat( "The number of attempts trying to generate a wanted point was too high" );
  }

  //LEAVE();
}


/*! Generates random vertices into the trimesh without any condition.
****************************************************************/
template<class CoordType>
Vertices<CoordType> TriMeshSpatialModel<CoordType>::randomVertices()
{
  ENTER("void TriMeshSpatialModel<CoordType>::randomVertices()");

  Vertices<CoordType> vertices( 3, _numCompartments, 0, 0 );

  for (int i = 0; i < _numCompartments; ++i)
  {
    drawPosition( vertices[i] );
  }

  LEAVE();

  return vertices;
}

/*! Generates vertices into the trimesh with a fixed distance to the border.
****************************************************************/
template<class CoordType>
Vertices<CoordType> TriMeshSpatialModel<CoordType>::distanceToTheBorder()
{
  //ENTER( "void TriMeshSpatialModel<CoordType>::distanceToTheBorder()" );
  Vertices<CoordType> vertices( 3, _numCompartments, 0, 0 );

  if ( _numCompartments != 1 && _distanceToBorder.getSize() == 1 )
  {
    float temp = _distanceToBorder[0];
    _distanceToBorder.setSize( _numCompartments );
    _distanceToBorder.fill( temp );
  }
  else if ( _numCompartments != _distanceToBorder.getSize() )
  {
    ProgramError error;
    error.setWhere( "void TriMeshSpatialModel<CoordType>::distanceToTheBorder()" );
    error.setWhat( "The number of compartments and the vector length of distances to the border are not the same" );
    return vertices;
    //LEAVE();
  }

  //calls the random drawing function as many times as the number of compartments we want
  for (int i = 0; i < _numCompartments; ++i)
    drawPositionFromBorder( vertices[i] , _distanceToBorder[i] );

  //LEAVE();
  return vertices;
}

/*! Generates vertices into the trimesh with a fixed hardcore distance between compartments.
****************************************************************/
template<class CoordType>
Vertices<CoordType> TriMeshSpatialModel<CoordType>::hardcoreDistance()
{
  //ENTER("void TriMeshSpatialModel<CoordType>::hardcoreDistance()");

  Vertices<CoordType> vertices( 3, 0, 0, 0 );
  Vector<CoordType> vertex( 3 );
  int attempts;
  const int maxAttempts = 10000 * _numCompartments;

  if ( _numCompartments != 1 && _hardcoreDistance.getSize() == 1 )
  {
    float temp = _hardcoreDistance[0];
    _hardcoreDistance.setSize( _numCompartments );
    _hardcoreDistance.fill( temp );
  }
  else if ( _numCompartments != _hardcoreDistance.getSize() )
  {
    ProgramError error;
    error.setWhere( "void TriMeshSpatialModel<CoordType>::hardcoreDistance()" );
    error.setWhat( "The number of compartments and the vector length of hardcore distances are not the same" );
    //LEAVE();
    return vertices;
  }

  for (int i = 0; i < _numCompartments; ++i)
  {
    attempts = 0;

    do
    { //checks the current vertex with previous vertices if hardcores distances are respected
      drawPosition( vertex );
    } while ( checkHardcoreDistances(vertex,vertices) == false && ++attempts < maxAttempts );

    //if the checking got true then the current vertex is included into the list
    if ( attempts < maxAttempts )
      vertices.append( vertex );
  }

  if ( attempts > maxAttempts )
  {
    Exception exception;
    exception.setWhere( "void TriMeshSpatialModel<CoordType>::hardcoreDistance()" );
    exception.setWhat( "The number of attempts trying to generate a wanted point was too high" );
  }

  //LEAVE();

  return vertices;
}

/*! Generates vertices into the trimesh with a fixed hardcore distance between compartments
 * and fixed distances to the border.
****************************************************************/
template<class CoordType>
Vertices<CoordType> TriMeshSpatialModel<CoordType>::hardcoreAndToTheBorderDistances()
{
  //ENTER("void TriMeshSpatialModel<CoordType>::hardcoreAndToTheBorderDistances()");
  Vertices<CoordType> vertices( 3, 0, 0, 0 );
  Vector<CoordType> vertex(3);

  int attempts = 0;
  const int maxAttempts = 10000;

  if ( _numCompartments != 1 && _hardcoreDistance.getSize() == 1 && _distanceToBorder.getSize() == 1 )
  {
    float tempHD = _hardcoreDistance[0];
    float tempDB = _distanceToBorder[0];
    _hardcoreDistance.setSize( _numCompartments );
    _distanceToBorder.setSize( _numCompartments );
    _hardcoreDistance.fill( tempHD );
    _distanceToBorder.fill( tempDB );
  }
  else if (  ( _numCompartments != _hardcoreDistance.getSize() && _hardcoreDistance.getSize() != 1 ) ||
             ( _numCompartments != _distanceToBorder.getSize() && _distanceToBorder.getSize() != 1 ) )
  {
    ProgramError error;
    error.setWhere( "void TriMeshSpatialModel<CoordType>::hardcoreAndToTheBorderDistances()" );
    error.setWhat( "The number of compartments and vectors' lengths of hardcore and to the border distances are not the same" );
    return vertices;
    //LEAVE();
  }

  for (int i = 0; i < _numCompartments; ++i)
  {
    attempts = 0;
    do
    {
      //checks if both conditions are true at the same time
      drawPositionFromBorder( vertex , _distanceToBorder[i] );
    } while ( checkHardcoreDistances(vertex,vertices) == false && ++attempts < maxAttempts );

  vertices.append( vertex );

  }

  if ( attempts > maxAttempts )
  {
    Exception exception;
    exception.setWhere( "void TriMeshSpatialModel<CoordType>::hardcoreAndToTheBorderDistances()" );
    exception.setWhat( "The number of attempts trying to generate a wanted point was too high" );
  }

  return vertices;
  //LEAVE();
}

/*! UPDATE: this check is already included during the "random generator with distance to the border" call.
 * Checks if the distance between a vertex and its closest point of the triMesh is correct.
****************************************************************/
template<class CoordType>
bool TriMeshSpatialModel<CoordType>::checkDistancesToBorder(const Vector<CoordType>& vertex, CoordType distanceToBorder)
{
  Vector<CoordType> triMeshVertex = vertex;
//  _triMesh->closestPoint( vertex, triMeshVertex );
  _triMeshQuery.getTriMesh().closestPoint( vertex, triMeshVertex );

  if ( abs( abs(triMeshVertex.distance(vertex)) - distanceToBorder ) < vertex.epsilon() ) return true;
  else return false;
}

/*! Checks if the hardcore distances from a vertex to previous generated vertices are correct.
****************************************************************/
template<class CoordType>
bool TriMeshSpatialModel<CoordType>::checkHardcoreDistances(const Vector<CoordType>& vertex, const Vertices<CoordType>& vertices)
{
  bool checkDistance = true;
  const int k = vertices.getNumVertices();
  Vector<float> triMeshVertex;

  for (int i = 0; checkDistance == true && i < vertices.getNumVertices(); ++i)
    if ( vertices[i].distance( vertex ) > _hardcoreDistance[i] + _hardcoreDistance[k] )
    {
      checkDistance = true;
      if ( checkDistance == true )
      {
        //_triMesh->closestPoint( vertex, triMeshVertex );
        _triMeshQuery.getTriMesh().closestPoint( vertex, triMeshVertex );
        if ( vertex.distance( triMeshVertex ) > _hardcoreDistance[k] )
          checkDistance = true;
        else checkDistance = false;
      }
    }
    else checkDistance = false;

  return checkDistance;

#if 0
  bool checkDistance = true;

  for (int i = 0; checkDistance && i < vertices.getNumVertices(); ++i)
    if ( ( _hardcoreDistance[0] != 0 ) && ( _hardcoreDistanceRange[1] == 0 ) )
    {
      if ( vertices[i].distance( vertex ) > _hardcoreDistance[i] + _hardcoreDistance[vertices.getNumVertices()] )
        checkDistance = true;
      else checkDistance = false;
    }
    else if ( ( _hardcoreDistance[0] == 0) && ( _hardcoreDistanceRange[1] != 0 ) )
    {
      if ( ( vertex.distance( vertices[jj] ) > 2 * _hardcoreDistanceRange[0] ) &&
           ( vertex.distance( vertices[jj] ) < 2 * _hardcoreDistanceRange[1] ) )
        checkDistance = true;
      else checkDistance = false;
    }
  return checkDistance;
#endif
}

/*! Saves generated vertices.
****************************************************************/
template<class CoordType>
void TriMeshSpatialModel<CoordType>::save(const string outputDir)
{
    _vertices.save( outputDir + "00", true );
}

//template class TriMeshSpatialModel<double>;
template class TriMeshSpatialModel<float>;
//template class TriMeshSpatialModel<int>;
