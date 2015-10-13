#include "trimeshspatialmodel.h"
#include "maximalrepulsion.h"

#include <randomgenerator.h>
#include <boundingbox.h>
#include <trimesh.h>
#include <alinevector.h>
#include <cmath>

#include <trimeshquery.h>

#include <programerror.h>
#include <dataset.h>

//#define TRACE
#include <trace.h>
#include <sstream>

#include <stopwatch.h>

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
  _hardcoreDistances.setZeros(1);
  _hardcoreDistancesRange.setZeros(2);
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
void TriMeshSpatialModel<CoordType>::setHardcoreDistances(const CoordType hardcoreDistances)
{
  _hardcoreDistances[0] = hardcoreDistances;
}

/*! Sets a vector of hardcore distances.
****************************************************************/
template<class CoordType>
void TriMeshSpatialModel<CoordType>::setHardcoreDistances(const Vector<CoordType>& hardcoreDistancess)
{
  _hardcoreDistances = hardcoreDistancess;
}

/*! Gets the correspondent hardcore distances of element n.
****************************************************************/
template<class CoordType>
const CoordType& TriMeshSpatialModel<CoordType>::getHardcoreDistances(const int numCompartment) const
{
  return _hardcoreDistances[numCompartment-1];
}

/*! Gets a vector of hardcore distances.
****************************************************************/
template<class CoordType>
Vector<CoordType>& TriMeshSpatialModel<CoordType>::getHardcoreDistances()
{
  return _hardcoreDistances;
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

/*! Diverts to the correct function depending on the conditions.
****************************************************************/
template<class CoordType>
Vertices<CoordType> TriMeshSpatialModel<CoordType>::drawSample(const int numPoints)
{
  _numCompartments = numPoints;

  if ( ( _hardcoreDistances[0] == 0 ) && ( _hardcoreDistancesRange[1] == 0 ) &&
       ( _distanceToBorder[0] == 0 ) && ( _distanceToBorderRange[1] == 0 ) )

    return randomVertices();

  else if ( ( _hardcoreDistances[0] == 0 ) && ( _hardcoreDistancesRange[1] == 0 ) &&
          ( ( _distanceToBorder[0] != 0 ) || ( _distanceToBorderRange[1] != 0 ) ) )

    return distanceToTheBorder();

  else if ( ( ( _hardcoreDistances[0] != 0 ) || ( _hardcoreDistancesRange[1] != 0 ) ) &&
              ( _distanceToBorder[0] == 0 ) && ( _distanceToBorderRange[1] == 0 ) )

    return hardcoreDistances();

  else //if ( ( ( _hardcoreDistances[0] != 0 ) || ( _hardcoreDistancesRange[1] != 0 ) ) &&
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
  const int maxAttempts = 200;
  bool found = false;

  do
  {
    drawPosition( position );

    //from the random point we find the closest point to it in the triMesh
    //the segment between them will be the normal vector
    //from this normal vector and using the distance to the border we obtain
    //the "position" vertex which will be our wanted vertex

    _triMeshQuery.closestPoint( position, triMeshVertex );

    normalVector = position - triMeshVertex;
    normalVector /= position.distance( triMeshVertex );
    normalVector *= distanceToBorder;
    position = triMeshVertex + normalVector;

    //we check if the obtained point is inside of the triMesh
    //and if it keeps the condition of being the minimum distance to the border
    //the distance to the border introduced manually before

    if ( _triMeshQuery.contains(position) )
      found = fabs( distanceToBorder - _triMeshQuery.closestPoint( position, triMeshVertex ) ) < position.epsilon();

  } while ( found == false && ++attempts <= maxAttempts );

  //EVAL(attempts);
  if ( attempts > maxAttempts )
  {
    //EVAL(attempts);
    Exception exception;
    exception.setWhere( "void TriMeshSpatialModel<CoordType>::drawPositionFromBorder()" );
    exception.setWhat( "The number of attempts trying to generate a wanted point was too high" );
    throw exception;
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
  RandomGenerator& randomGenerator = this->getRandomGenerator();

  if ( _numCompartments != 1 && _distanceToBorder.getSize() == 1 && _distanceToBorderRange[1] == 0 )
  {
    float temp = _distanceToBorder[0];
    _distanceToBorder.setSize( _numCompartments );
    _distanceToBorder.fill( temp );
  }
  if ( _numCompartments != 0 && _distanceToBorder.getSize() != 0 && _distanceToBorderRange[1] != 0 )
  {
      _distanceToBorder.setSize( _numCompartments );
      for ( int i = 0; i < _numCompartments; ++i )
      {
          float temp = randomGenerator.uniformLF(_distanceToBorderRange[0],_distanceToBorderRange[1]);
          _distanceToBorder[i] = temp;
      }
  }
  else if ( _numCompartments != _distanceToBorder.getSize() )
  {
    ProgramError programError;
    programError.setWhere( "void TriMeshSpatialModel<CoordType>::hardcoreDistances()" );
    programError.setWhat( "The number of compartments and the vector length of hardcore distances are not the same" );
    throw programError;
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
Vertices<CoordType> TriMeshSpatialModel<CoordType>::hardcoreDistances()
{
  //ENTER("void TriMeshSpatialModel<CoordType>::hardcoreDistances()");

  Vertices<CoordType> vertices( 3, 0, 0, 0 );
  Vector<CoordType> vertex( 3 );
//  RandomGenerator& randomGenerator = this->getRandomGenerator();

  //Vector<CoordType> test( 3 );
  int attempts;
  const int maxAttempts = 150;

  if ( _numCompartments != 1 && _hardcoreDistances.getSize() == 1 )
  {
    float temp = _hardcoreDistances[0];
    _hardcoreDistances.setSize( _numCompartments );
    _hardcoreDistances.fill( temp );
  }
  else if ( _numCompartments != _hardcoreDistances.getSize() )
  {
    ProgramError programError;
    programError.setWhere( "void TriMeshSpatialModel<CoordType>::hardcoreDistances()" );
    programError.setWhat( "The number of compartments and the vector length of hardcore distances are not the same" );
    throw programError;
  }

  //randomizes the order of the compartments
//  randomizesOrder(randomGenerator);
//  EVAL(_hardcoreDistances);
  for (int i = 0; i < _numCompartments; ++i)
  {
    attempts = 0;

    do
    { //checks the current vertex with previous vertices if hardcores distances are respected
      drawPosition( vertex );
      //EVAL(vertex);
    } while ( ( checkInterObjectDistances(vertex,vertices) == false
              || checkObjectToBorderDistance(vertex, i) == false )
      && ++attempts < maxAttempts );

    //EVAL(_triMeshQuery.getTriMesh().closestPoint(vertex,test));
    //if the checking got true then the current vertex is included into the list
    if ( attempts < maxAttempts )
      vertices.append( vertex );
    else
    {
      Exception exception;
      exception.setWhere( "void TriMeshSpatialModel<CoordType>::hardcoreDistances()" );
      exception.setWhat( "The number of attempts trying to generate a wanted point was too high" );
      throw exception;
    }

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
  RandomGenerator& randomGenerator = this->getRandomGenerator();

//  Stopwatch stopWatch;
//  stopWatch.start( "One pattern" );

  int attempts = 0;
  const int maxAttempts = 200;

  if ( _numCompartments != 1 && _hardcoreDistances.getSize() == 1 && _distanceToBorder.getSize() == 1 )
  {
    float tempHD = _hardcoreDistances[0];
    float tempDB = _distanceToBorder[0];
    _hardcoreDistances.setSize( _numCompartments );
    _distanceToBorder.setSize( _numCompartments );
    _hardcoreDistances.fill( tempHD );
    _distanceToBorder.fill( tempDB );
  }
  else if (  ( _numCompartments != _hardcoreDistances.getSize() && _hardcoreDistances.getSize() != 1 ) ||
             ( _numCompartments != _distanceToBorder.getSize() && _distanceToBorder.getSize() != 1 ) )
  {
    ProgramError programError;
    programError.setWhere( "void TriMeshSpatialModel<CoordType>::hardcoreAndToTheBorderDistances()" );
    programError.setWhat( "The number of compartments and vectors' lengths of hardcore and to the border distances are not the same" );
    throw programError;
    return vertices;
    //LEAVE();
  }

  //randomizes the order of the compartments
  randomizesOrder(randomGenerator);
  EVAL(_hardcoreDistances);

  for (int i = 0; i < _numCompartments; ++i)
  {
    attempts = 0;

    do
    {
      //checks if both conditions are true at the same time
      drawPositionFromBorder( vertex , _distanceToBorder[i] );

    } while ( checkInterObjectDistances(vertex,vertices) == false && ++attempts < maxAttempts );

    vertices.append( vertex );

    EVAL(attempts);
    if ( attempts >= maxAttempts )
    {
      //EVAL(attempts);
      Exception exception;
      exception.setWhere( "void TriMeshSpatialModel<CoordType>::hardcoreAndToTheBorderDistances()" );
      exception.setWhat( "The number of attempts trying to generate a wanted point was too high" );
      throw exception;
    }

  }

//  stopWatch.stop( "Fonction vertices" );
//  stopWatch.print();
  //LEAVE();
  return vertices;
}

///*! UPDATE: this check is already included during the "random generator with distance to the border" call.
// * Checks if the distance between a vertex and its closest point of the triMesh is correct.
//****************************************************************/
//template<class CoordType>
//bool TriMeshSpatialModel<CoordType>::checkDistancesToBorder(
//  const Vector<CoordType>& vertex,
//  CoordType distanceToBorder)
//{
//  Vector<CoordType> triMeshVertex = vertex;
//  //const CoordType distance = _triMeshQuery.getTriMesh().closestPoint( vertex, triMeshVertex );
//  return abs(_triMeshQuery.getTriMesh().closestPoint( vertex, triMeshVertex )-distanceToBorder) < vertex.epsilon();
//}



/////*! Checks if the hardcore distances from a vertex to previous generated vertices are correct.
//****************************************************************/
//template<class CoordType>
//bool TriMeshSpatialModel<CoordType>::checkHardcoreDistances(const Vector<CoordType>& vertex, const Vertices<CoordType>& vertices)
//{
//  const int k = vertices.getNumVertices();

//  bool checkDistance = true;
//  Vector<float> triMeshVertex;

//  checkDistance = _triMeshQuery.getTriMesh().closestPoint(vertex, triMeshVertex) > _hardcoreDistances[k];

//  for (int i = 0; checkDistance == true && i < vertices.getNumVertices(); ++i)
//    if ( vertices[i].distance( vertex ) < _hardcoreDistances[i] + _hardcoreDistances[k] )
//      return false;

//  return checkDistance;

// //  return checkHardcoreDistances( vertex, k, vertices );


//}*/

/*! Randomizes the position of the objects
****************************************************************/
template<class CoordType>
void TriMeshSpatialModel<CoordType>::randomizesOrder(RandomGenerator& randomGenerator)
{
 //RandomGenerator& randomGenerator = this->getRandomGenerator();

  //creates copies of the distances to the border and volumes vectors
  Vector<CoordType> hardcoreDistancesCopy (_hardcoreDistances);
  Vector<CoordType>  distanceToBorderCopy (_distanceToBorder);

  //to shuffle the order of the compartments to be generated, first we stablish the normal raising order
  Vector<int> shuffledCompartmentsIndexes (_numCompartments);
  for (int order = 0; order < _numCompartments; ++order)
    shuffledCompartmentsIndexes[order] = order;

  //reorganizes the order of the compartments
  for ( int p = 0; p < _numCompartments; ++p)
  {
    int k1 = randomGenerator.uniformL( _numCompartments );
    int k2 = randomGenerator.uniformL( _numCompartments );
    int temp = shuffledCompartmentsIndexes[k1];
    shuffledCompartmentsIndexes[k1] = shuffledCompartmentsIndexes[k2];
    shuffledCompartmentsIndexes[k2] = temp;
  }

  for (int i = 0; i < _numCompartments; ++i)
  {
    _hardcoreDistances[i] = hardcoreDistancesCopy[shuffledCompartmentsIndexes[i]];
    _distanceToBorder[i]  =  distanceToBorderCopy[shuffledCompartmentsIndexes[i]];
  }
}

/*! Checks the distance from the compartment to the border
****************************************************************/
template<class CoordType>
bool TriMeshSpatialModel<CoordType>::checkObjectToBorderDistance(const Vector<CoordType>& vertex, const int& numObject)
{
  Vector<float> triMeshVertex;
  return _triMeshQuery.closestPoint(vertex, triMeshVertex) > _hardcoreDistances[numObject];
}

/*! Checks if the hardcore distances from a vertex to previous generated vertices are correct.
****************************************************************/
template<class CoordType>
bool TriMeshSpatialModel<CoordType>::checkInterObjectDistances(const Vector<CoordType>& vertex, const Vertices<CoordType>& vertices)
{
  const int k = vertices.getNumVertices();

  for (int i = 0; i < vertices.getNumVertices(); ++i)
    if ( vertices[i].distance( vertex ) < (_hardcoreDistances[i] + _hardcoreDistances[k]) )
      return false;

  return true;
}

// /*! Checks if the hardcore distances between a new (moved) vertex and the others are respected.
//****************************************************************/
//template<class CoordType>
//bool TriMeshSpatialModel<CoordType>::checkHardcoreDistances(const Vector<CoordType>& vertex, const int numVertexChanged, const Vertices<CoordType>& vertices) const
//{
//  bool checkDistance;
//  //const int k = vertices.getNumVertices();
//  Vector<float> triMeshVertex;

//  checkDistance = _triMeshQuery.getTriMesh().closestPoint(vertex, triMeshVertex) > _hardcoreDistances[numVertexChanged];

//  for (int i = 0; checkDistance == true && i < vertices.getNumVertices(); ++i)
//    if ( i!= numVertexChanged)
//        if ( vertices[i].distance( vertex ) < _hardcoreDistances[i] + _hardcoreDistances[numVertexChanged] )
//          return false;

//  return checkDistance;

//#if 0
//  bool checkDistance = true;

//  for (int i = 0; checkDistance && i < vertices.getNumVertices(); ++i)
//    if ( ( _hardcoreDistances[0] != 0 ) && ( _hardcoreDistancesRange[1] == 0 ) )
//    {
//      if ( vertices[i].distance( vertex ) > _hardcoreDistances[i] + _hardcoreDistances[vertices.getNumVertices()] )
//        checkDistance = true;
//      else checkDistance = false;
//    }
//    else if ( ( _hardcoreDistances[0] == 0) && ( _hardcoreDistancesRange[1] != 0 ) )
//    {
//      if ( ( vertex.distance( vertices[jj] ) > 2 * _hardcoreDistancesRange[0] ) &&
//           ( vertex.distance( vertices[jj] ) < 2 * _hardcoreDistancesRange[1] ) )
//        checkDistance = true;
//      else checkDistance = false;
//    }
//  return checkDistance;
//#endif
//}


//template class TriMeshSpatialModel<double>;
template class TriMeshSpatialModel<float>;
//template class TriMeshSpatialModel<int>;
