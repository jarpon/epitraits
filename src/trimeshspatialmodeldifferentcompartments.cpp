#include "trimeshspatialmodeldifferentcompartments.h"
#include "maximalrepulsion.h"

#include <randomgenerator.h>
#include <boundingbox.h>
#include <trimesh.h>
#include <alinevector.h>
#include <cmath>

#include <trimeshquery.h>

#include <programerror.h>
#include <dataset.h>

#define TRACE
#include <trace.h>
#include <sstream>

#include <stopwatch.h>

using namespace std;



/*! Initializes all the parameters.
****************************************************************/
template<class CoordType>
TriMeshSpatialModelDifferentCompartments<CoordType>::TriMeshSpatialModelDifferentCompartments()
{
}

/*! Destroys it.
****************************************************************/
template<class CoordType>
TriMeshSpatialModelDifferentCompartments<CoordType>::~TriMeshSpatialModelDifferentCompartments()
{
}

/*! Sets the triMesh to work with.
****************************************************************/
template<class CoordType>
void TriMeshSpatialModelDifferentCompartments<CoordType>::setTriMesh(const TriMesh<CoordType>& triMesh)
{
  _triMeshQuery.setTriMesh( triMesh );
}

template<class CoordType>
const TriMesh<CoordType>& TriMeshSpatialModelDifferentCompartments<CoordType>::getTriMesh() const
{
  return _triMeshQuery.getTriMesh();
}

template<class CoordType>
const TriMeshQuery<CoordType>& TriMeshSpatialModelDifferentCompartments<CoordType>::getTriMeshQuery() const
{
  return _triMeshQuery;
}


/*! Initializes the class.
****************************************************************/
template<class CoordType>
void TriMeshSpatialModelDifferentCompartments<CoordType>::initialize()
{
  _boundingBox = _triMeshQuery.getTriMesh().boundingBox();

#if 0
  //shows the maximum distance to the border that can be used
  float maxDistance = ( boundingBox.max(0) - boundingBox.min(0) )/2;
  for ( int i = 1; i < 3; ++i)
  {
    if ( ( boundingBox.max(i) - boundingBox.min(i) )/2 < maxDistance )
      maxDistance = ( boundingBox.max(i) - boundingBox.min(i) )/2;
  }
  EVAL(maxDistance);
#endif
}

/*! Generates a random vertex into the trimesh.
****************************************************************/
template<class CoordType>
void TriMeshSpatialModelDifferentCompartments<CoordType>::drawPosition(Vector<CoordType>& position)
{
  RandomGenerator& randomGenerator = this->getRandomGenerator();
  const Vector<CoordType> v1 = _boundingBox.getVertex1();
  const Vector<CoordType> v2 = _boundingBox.getVertex2();

  do
  {
    position[X] = randomGenerator.uniformLF( v1[X], v2[X] );
    position[Y] = randomGenerator.uniformLF( v1[Y], v2[Y] );
    position[Z] = randomGenerator.uniformLF( v1[Z], v2[Z] );
  } while ( _triMeshQuery.contains(position) == false );
}


#if 0

/*! Sets the number of compartments to generate
****************************************************************/
template<class CoordType>
void TriMeshSpatialModelDifferentCompartments<CoordType>::setNumberOfCompartments(const int numCompartments)
{
  _numCompartments = numCompartments;
}

/*! Sets a unique hardcore distance for all compartments.
****************************************************************/
template<class CoordType>
void TriMeshSpatialModelDifferentCompartments<CoordType>::setHardcoreDistances(const CoordType hardcoreDistances)
{
  _hardcoreDistances[0] = hardcoreDistances;
}

/*! Sets a vector of hardcore distances.
****************************************************************/
template<class CoordType>
void TriMeshSpatialModelDifferentCompartments<CoordType>::setHardcoreDistances(const Vector<CoordType>& hardcoreDistancess)
{
  _hardcoreDistances = hardcoreDistancess;
}

/*! Gets the correspondent hardcore distances of element n.
****************************************************************/
template<class CoordType>
const CoordType& TriMeshSpatialModelDifferentCompartments<CoordType>::getHardcoreDistances(const int numCompartment) const
{
  return _hardcoreDistances[numCompartment-1];
}

/*! Gets a vector of hardcore distances.
****************************************************************/
template<class CoordType>
Vector<CoordType>& TriMeshSpatialModelDifferentCompartments<CoordType>::getHardcoreDistances()
{
  return _hardcoreDistances;
}

/*! Gets the number of different compartments within the scenario.
****************************************************************/
template<class CoordType>
int TriMeshSpatialModelDifferentCompartments<CoordType>::getNumDistributions()
{
  return _numDistributions;
}

/*! Sets a unique distance to the border for all compartments.
****************************************************************/
template<class CoordType>
void TriMeshSpatialModelDifferentCompartments<CoordType>::setDistanceToBorder(const CoordType distanceToBorder)
{
  _distanceToBorder[0] = distanceToBorder;
}

/*! Sets a vector of distances to the border.
****************************************************************/
template<class CoordType>
void TriMeshSpatialModelDifferentCompartments<CoordType>::setDistanceToBorder(const Vector<CoordType>& distancesToTheBorder)
{
  _distanceToBorder = distancesToTheBorder;
}

/*! Diverts to the correct function depending on the conditions.
 *  First version.
 *  Currently, working as always with just one kind of distribution.
****************************************************************/
template<class CoordType>
Vertices<CoordType> TriMeshSpatialModelDifferentCompartments<CoordType>::drawSample(const int numPoints)
{
  ENTER("void TriMeshSpatialModelDifferentCompartments<CoordType>::drawSample( ... )")

  if ( _numDistributions == 1 )
  {
    if ( _numCompartments == 0 )
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
  else drawSample2(  );

  LEAVE();
}

/*! Diverts to the correct function depending on the conditions.
 *  New version:
 *  First, it randomizes the order of the objects, mixing all objects independtently of their kind.
 *  The randomization changes all the vectors information to keep the proper order.
 *  Then it calls the correct model.
 *  The output is an updated vertices variable with all objects and the corresponding vector with their objects association.
****************************************************************/
template<class CoordType>
void TriMeshSpatialModelDifferentCompartments<CoordType>::drawSample( Vertices<CoordType>& vertices, Vector<string>& order )
{
  RandomGenerator& randomGenerator = this->getRandomGenerator();
  randomizesOrder(randomGenerator);

  order = allDistributionObjectsNames;

  if ( ( _hardcoreDistances[0] == 0 ) && ( _hardcoreDistancesRange[1] == 0 ) &&
       ( _distanceToBorder[0] == 0 ) && ( _distanceToBorderRange[1] == 0 ) )

    vertices = randomVertices();

  else if ( ( _hardcoreDistances[0] == 0 ) && ( _hardcoreDistancesRange[1] == 0 ) &&
          ( ( _distanceToBorder[0] != 0 ) || ( _distanceToBorderRange[1] != 0 ) ) )

    vertices = distanceToTheBorder();

  else if ( ( ( _hardcoreDistances[0] != 0 ) || ( _hardcoreDistancesRange[1] != 0 ) ) &&
              ( _distanceToBorder[0] == 0 ) && ( _distanceToBorderRange[1] == 0 ) )

    vertices = hardcoreDistances();

  else //if ( ( ( _hardcoreDistances[0] != 0 ) || ( _hardcoreDistancesRange[1] != 0 ) ) &&
       //     ( ( _distanceToBorder[0] != 0 ) || ( _distanceToBorderRange[1] != 0 ) ) )

    vertices = hardcoreAndToTheBorderDistances();

}


/*! Diverts to the correct function depending on the conditions.
 *  New version for two different organizations:
 *  First, it randomizes the order of the objects, mixing all objects independtently of their kind.
 *  The randomization changes all the vectors information to keep the proper order.
 *  Then it calls the correct model.
 *  The output is two updated vertices with the corresponding objects of each class.
****************************************************************/
template<class CoordType>
Vertices<CoordType> TriMeshSpatialModelDifferentCompartments<CoordType>::drawSample2( )
{
//  RandomGenerator& randomGenerator = this->getRandomGenerator();
//  randomizesOrder(randomGenerator);

//  order = allDistributionObjectsNames;

  EVAL ("this is going well");
  EVAL (_allDistributionHardcoreDistances);
  EVAL (_allDistributionNumObjects.getSize());
  EVAL (_numDistributions);
  _hardcoreDistances = _allDistributionHardcoreDistances;
  _numCompartments = _allDistributionNumObjects.sum();
  EVAL(_hardcoreDistances);
  return hardcoreDistances();

  if ( ( _hardcoreDistances[0] == 0 ) && ( _hardcoreDistancesRange[1] == 0 ) &&
       ( _distanceToBorder[0] == 0 ) && ( _distanceToBorderRange[1] == 0 ) )

    //randomVertices(); //remains ordering

  {

      Vertices<CoordType> tempAllVertices = randomVertices();

//      for ( int i = 0; i < _allDistributionNumObjects.getSize(); i++ )
//      {
//          if ( allDistributionObjectsNames[i] == _objectsNames[1] )
//              dist1.append( tempAllVertices[i] );
//          else
//              dist2.append( tempAllVertices[i] );
//      }

  }

  else if ( ( _hardcoreDistances[0] == 0 ) && ( _hardcoreDistancesRange[1] == 0 ) &&
          ( ( _distanceToBorder[0] != 0 ) || ( _distanceToBorderRange[1] != 0 ) ) )
    {

        Vertices<CoordType> tempAllVertices = distanceToTheBorder();

//        for ( int i = 0; i < _allDistributionNumObjects.getSize(); i++ )
//        {
//            if ( allDistributionObjectsNames[i] == _objectsNames[1] )
//                dist1.append( tempAllVertices[i] );
//            else
//                dist2.append( tempAllVertices[i] );
//        }

    }

  else if ( ( ( _hardcoreDistances[0] != 0 ) || ( _hardcoreDistancesRange[1] != 0 ) ) &&
              ( _distanceToBorder[0] == 0 ) && ( _distanceToBorderRange[1] == 0 ) )

  {
    ENTER("void TriMeshSpatialModelDifferentCompartments<CoordType>::drawSample2( ... )")

      hardcoreDistances();

//      for ( int i = 0; i < _allDistributionNumObjects.getSize(); i++ )
//      {
//          if ( allDistributionObjectsNames[i] == _objectsNames[1] )
//              dist1.append( tempAllVertices[i] );
//          else
//              dist2.append( tempAllVertices[i] );
//      }

    LEAVE();
  }

  else //if ( ( ( _hardcoreDistances[0] != 0 ) || ( _hardcoreDistancesRange[1] != 0 ) ) &&
       //     ( ( _distanceToBorder[0] != 0 ) || ( _distanceToBorderRange[1] != 0 ) ) )
  {
      Vertices<CoordType> tempAllVertices = hardcoreAndToTheBorderDistances();

//      for ( int i = 0; i < _allDistributionNumObjects.getSize(); i++ )
//      {
//          if ( allDistributionObjectsNames[i] == _objectsNames[1] )
//              dist1.append( tempAllVertices[i] );
//          else
//              dist2.append( tempAllVertices[i] );
//      }
  }
}

/*! Generates a random vertex into the trimesh.
****************************************************************/
template<class CoordType>
void TriMeshSpatialModelDifferentCompartments<CoordType>::drawPosition(Vector<CoordType>& position)
{
  ENTER( "void TriMeshSpatialModelDifferentCompartments<CoordType>::drawPosition(Vector<CoordType>&)" );

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
void TriMeshSpatialModelDifferentCompartments<CoordType>::drawPositionFromBorder(Vector<CoordType>& position, CoordType distanceToBorder)
{
  //ENTER("void TriMeshSpatialModelDifferentCompartments<CoordType>::drawPositionFromBorder(Vector<CoordType,PixelType>, T)");
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
    exception.setWhere( "void TriMeshSpatialModelDifferentCompartments<CoordType>::drawPositionFromBorder()" );
    exception.setWhat( "The number of attempts trying to generate a wanted point was too high" );
    throw exception;
  }

  //LEAVE();
}

/*! Generates random vertices into the trimesh without any condition.
****************************************************************/
template<class CoordType>
Vertices<CoordType> TriMeshSpatialModelDifferentCompartments<CoordType>::randomVertices()
{
  ENTER("void TriMeshSpatialModelDifferentCompartments<CoordType>::randomVertices()");

  Vertices<CoordType> vertices( 3, _numCompartments, 0, 0 );

  for (int i = 0; i < _numCompartments; ++i)
  {
    drawPosition( vertices[i] );
  }

  return vertices;

  LEAVE();

}

/*! With the purpose of analysing more than one kind of objects distribution
 * this saves the current (latest) info of the objects distribution in global variables.
****************************************************************/
template<class CoordType>
void TriMeshSpatialModelDifferentCompartments<CoordType>::addDistribution()
{
  ENTER("void TriMeshSpatialModelDifferentCompartments<CoordType>::addDistribution()");

    _allDistributionNumObjects[_allDistributionNumObjects.getSize()+1] = _numCompartments;
    _allDistributionHardcoreDistances.append( _hardcoreDistances );
    _allDistributionDistancesToTheBorder.append( _distanceToBorder );
    _numDistributions++;

  LEAVE();
}

/*! With the purpose of analysing more than one kind of objects distribution
 * this saves the current (latest) info of the objects distribution in global variables.
 * Adding the name of the specific kind of compartments/objects
****************************************************************/
template<class CoordType>
void TriMeshSpatialModelDifferentCompartments<CoordType>::addDistribution( const string& objectsName )
{
  ENTER("void TriMeshSpatialModelDifferentCompartments<CoordType>::addDistribution()");

//    _allDistributionNumObjects[_allDistributionNumObjects.getSize()+1] = _numCompartments;
    _allDistributionNumObjects[_allDistributionNumObjects.getSize()] = _numCompartments;
    _allDistributionHardcoreDistances.append( _hardcoreDistances );
    _allDistributionDistancesToTheBorder.append( _distanceToBorder );
//    _objectsNames[_objectsNames.getSize()+1] = objectsName;
    _objectsNames[_objectsNames.getSize()] = objectsName;

    for ( int i = 0; i < _numCompartments; i++ )
        //allDistributionObjectsNames[allDistributionObjectsNames.getSize()+1] = objectsName;
        allDistributionObjectsNames[allDistributionObjectsNames.getSize()] = objectsName;

    _numDistributions++;

  LEAVE();
}

/*! Adds the information of a new objects distribution:
 * number of objects and their radii
****************************************************************/
template<class CoordType>
void TriMeshSpatialModelDifferentCompartments<CoordType>::addDistribution( const int numObjects, const Vector<CoordType>& hardcoreDistances )
{
  ENTER("void TriMeshSpatialModelDifferentCompartments<CoordType>::addDistribution(...)");

  EVAL (_allDistributionNumObjects);
    //not the best way
    //_allDistributionNumObjects.setSize( _allDistributionNumObjects.getSize() + 1 );
    EVAL (_allDistributionNumObjects);
    Vector<int> temp (1);
    temp[0] = numObjects;

    _allDistributionNumObjects.append( temp );
    _allDistributionHardcoreDistances.append( hardcoreDistances );
    _numCompartments += numObjects;
    //_allDistributionDistancesToTheBorder.append( distancesToTheBorder );

    ++_numDistributions;

    EVAL (_allDistributionHardcoreDistances);
    EVAL (_allDistributionNumObjects);
    EVAL( _numCompartments );
    EVAL (_numDistributions);


  LEAVE();
}


/*! Adds the information of a new objects distribution:
 * number of objects, their radii and their distances to the border
****************************************************************/
template<class CoordType>
void TriMeshSpatialModelDifferentCompartments<CoordType>::addDistribution( const int numObjects, const string& objectsName, const Vector<CoordType>& hardcoreDistances, const Vector<CoordType>& distancesToTheBorder )
{
  ENTER("void TriMeshSpatialModelDifferentCompartments<CoordType>::addDistribution(...)");

//    _allDistributionNumObjects[_allDistributionNumObjects.getSize()+1] = _numCompartments;
    _allDistributionNumObjects[_allDistributionNumObjects.getSize()] = _numCompartments;
    _allDistributionHardcoreDistances.append( hardcoreDistances );
    _allDistributionDistancesToTheBorder.append( distancesToTheBorder );

    //_objectsNames[_objectsNames.getSize()+1] = objectsName;
    _objectsNames[_objectsNames.getSize()] = objectsName;

    for ( int i = 1; i <= numObjects; i++ )
//        allDistributionObjectsNames[allDistributionObjectsNames.getSize()+1] = objectsName;
    allDistributionObjectsNames[allDistributionObjectsNames.getSize()] = objectsName;

    _numDistributions++;

  LEAVE();
}

/*! Sets the names of the different kind of objects
****************************************************************/
template<class CoordType>
void TriMeshSpatialModelDifferentCompartments<CoordType>::addObjectsNames( const Vector<string>& objectsNames )
{
  ENTER("void TriMeshSpatialModelDifferentCompartments<CoordType>::addObjectsNames(...)");

    _objectsNames = objectsNames;

    for ( int j = 1; j <= _numDistributions; j++ )
        for ( int i = 1; i <= _numCompartments; i++ )
            //allDistributionObjectsNames[allDistributionObjectsNames.getSize()+1] = objectsNames[j];
            allDistributionObjectsNames[allDistributionObjectsNames.getSize()] = objectsNames[j];

  LEAVE();
}

/*! Sets the name of the selected kind of objects
****************************************************************/
template<class CoordType>
void TriMeshSpatialModelDifferentCompartments<CoordType>::addObjectsName( const int kindObject, const string& objectsName )
{
  ENTER("void TriMeshSpatialModelDifferentCompartments<CoordType>::addObjectsName(...)");

    _objectsNames[kindObject-1] = objectsName;

  LEAVE();
}

/*! Returns a vector with the names of the different kind of objects
****************************************************************/
template<class CoordType>
Vector<string> TriMeshSpatialModelDifferentCompartments<CoordType>::getObjectsNames()
{
  ENTER("void TriMeshSpatialModelDifferentCompartments<CoordType>::getObjectsNames()");

    return _objectsNames;

  LEAVE();
}

/*! Returns the name of the selected kind of objects
****************************************************************/
template<class CoordType>
string TriMeshSpatialModelDifferentCompartments<CoordType>::getObjectsName( const int kindObject )
{
  ENTER("void TriMeshSpatialModelDifferentCompartments<CoordType>::getObjectsName(...)");

    return _objectsNames[kindObject-1];

  LEAVE();
}

/*! Returns the name of the selected kind of objects
****************************************************************/
template<class CoordType>
int TriMeshSpatialModelDifferentCompartments<CoordType>::getObjectsNumber( const string& kindObject )
{
  ENTER("void TriMeshSpatialModelDifferentCompartments<CoordType>::getObjectsNumber(...)");

    for ( int i = 0; i < _objectsNames.getSize(); i++ )
        if ( _objectsNames[i] == kindObject )
            return i;

  LEAVE();
}

/*! Gets the information of the requested objects distribution.
****************************************************************/
template<class CoordType>
void TriMeshSpatialModelDifferentCompartments<CoordType>::getDistribution( const int kindObjects, int numObjects, Vector<CoordType> hardcoreDistances, Vector<CoordType> distancesToTheBorder )
{
  ENTER("void TriMeshSpatialModelDifferentCompartments<CoordType>::getDistribution(...)");

  if ( kindObjects == 1 || _numDistributions == 1 )
    {
        numObjects = _numCompartments;
        hardcoreDistances = _hardcoreDistances;
        distancesToTheBorder = _distanceToBorder;
    }
    else
    {
        int temp = 0;
        for (int i = 0; i < kindObjects; i++)
          temp += _allDistributionNumObjects[i];

        numObjects = _allDistributionNumObjects[kindObjects-1];

        for (int j = 0; j < kindObjects; j++)
        {
            hardcoreDistances[j] = _allDistributionHardcoreDistances[temp+j];
            distancesToTheBorder[j] = _allDistributionDistancesToTheBorder[temp+j];
        }

    }

  LEAVE();
}

/*! Generates vertices into the trimesh with a fixed distance to the border.
****************************************************************/
template<class CoordType>
Vertices<CoordType> TriMeshSpatialModelDifferentCompartments<CoordType>::distanceToTheBorder()
{
  //ENTER( "void TriMeshSpatialModelDifferentCompartments<CoordType>::distanceToTheBorder()" );
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
    programError.setWhere( "void TriMeshSpatialModelDifferentCompartments<CoordType>::hardcoreDistances()" );
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
Vertices<CoordType> TriMeshSpatialModelDifferentCompartments<CoordType>::hardcoreDistances()
{
  //ENTER("void TriMeshSpatialModelDifferentCompartments<CoordType>::hardcoreDistances()");
  EVAL(_hardcoreDistances);
  EVAL(_numCompartments);
  Vertices<CoordType> vertices( 3, 0, 0, 0 );
  Vector<CoordType> vertex( 3 );
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
    programError.setWhere( "void TriMeshSpatialModelDifferentCompartments<CoordType>::hardcoreDistances()" );
    programError.setWhat( "The number of compartments and the vector length of hardcore distances are not the same" );
    throw programError;
  }

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
      exception.setWhere( "void TriMeshSpatialModelDifferentCompartments<CoordType>::hardcoreDistances()" );
      exception.setWhat( "The number of attempts trying to generate a wanted point was too high" );
      throw exception;
    }

  }

  //LEAVE();

  for ( int i = 0; i < _numCompartments; ++i)
  {
    EVAL(vertices[i]);
  }

  return vertices;
}

/*! Generates vertices into the trimesh with a fixed hardcore distance between compartments
 * and fixed distances to the border.
****************************************************************/
template<class CoordType>
Vertices<CoordType> TriMeshSpatialModelDifferentCompartments<CoordType>::hardcoreAndToTheBorderDistances()
{
  //ENTER("void TriMeshSpatialModelDifferentCompartments<CoordType>::hardcoreAndToTheBorderDistances()");
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
    programError.setWhere( "void TriMeshSpatialModelDifferentCompartments<CoordType>::hardcoreAndToTheBorderDistances()" );
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
    if ( attempts > maxAttempts )
    {
      //EVAL(attempts);
      Exception exception;
      exception.setWhere( "void TriMeshSpatialModelDifferentCompartments<CoordType>::hardcoreAndToTheBorderDistances()" );
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
//bool TriMeshSpatialModelDifferentCompartments<CoordType>::checkDistancesToBorder(
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
//bool TriMeshSpatialModelDifferentCompartments<CoordType>::checkHardcoreDistances(const Vector<CoordType>& vertex, const Vertices<CoordType>& vertices)
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
void TriMeshSpatialModelDifferentCompartments<CoordType>::randomizesOrder(RandomGenerator& randomGenerator)
{
 //RandomGenerator& randomGenerator = this->getRandomGenerator();

  //creates copies of the distances to the border and volumes vectors
  Vector<CoordType> hardcoreDistancesCopy ( _hardcoreDistances );
  Vector<CoordType>  distanceToBorderCopy ( _distanceToBorder );
  Vector<string>  allDistributionObjectsNamesCopy ( allDistributionObjectsNames );

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
    allDistributionObjectsNames[i] = allDistributionObjectsNamesCopy[shuffledCompartmentsIndexes[i]];
  }
}

/*! Checks the distance from the compartment to the border
****************************************************************/
template<class CoordType>
bool TriMeshSpatialModelDifferentCompartments<CoordType>::checkObjectToBorderDistance(const Vector<CoordType>& vertex, const int& numObject)
{
  Vector<float> triMeshVertex;
  return _triMeshQuery.closestPoint(vertex, triMeshVertex) > _hardcoreDistances[numObject];
}

/*! Checks if the hardcore distances from a vertex to previous generated vertices are correct.
****************************************************************/
template<class CoordType>
bool TriMeshSpatialModelDifferentCompartments<CoordType>::checkInterObjectDistances(const Vector<CoordType>& vertex, const Vertices<CoordType>& vertices)
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
//bool TriMeshSpatialModelDifferentCompartments<CoordType>::checkHardcoreDistances(const Vector<CoordType>& vertex, const int numVertexChanged, const Vertices<CoordType>& vertices) const
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

#endif


template class TriMeshSpatialModelDifferentCompartments<double>;
template class TriMeshSpatialModelDifferentCompartments<float>;
template class TriMeshSpatialModelDifferentCompartments<long double>;

