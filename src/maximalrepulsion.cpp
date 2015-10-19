#include <trimeshspatialmodel.h>
#include "maximalrepulsion.h"

#include <randomgenerator.h>
#include <trimesh.h>
#include <alinevector.h>
#include <cmath>

#include <trimeshquery.h>

#include <programerror.h>
#include <dataset.h>

#include <curvestack.h>
#include <vertexstack.h>

#define TRACE
#include <trace.h>
#include <sstream>

#include <limits>       // std::numeric_limits

using namespace std;

/*! Initializes all the parameters.
****************************************************************/
template<class CoordType>
MaximalRepulsionTriMeshSpatialModel<CoordType>::MaximalRepulsionTriMeshSpatialModel()
{
}

/*! Destroys it.
****************************************************************/
template<class CoordType>
MaximalRepulsionTriMeshSpatialModel<CoordType>::~MaximalRepulsionTriMeshSpatialModel()
{
}

/*! Sets a vector of hardcore distances.
****************************************************************/
template<class CoordType>
void MaximalRepulsionTriMeshSpatialModel<CoordType>::setHardcoreDistance(const Vector<CoordType>& hardcoreDistances)
{
  _hardcoreDistances = hardcoreDistances;
}

// /*! Sets the triMeshSpatialModel to work with ------ (not needed)
//****************************************************************/
//template<class CoordType>
//void MaximalRepulsionTriMeshSpatialModel<CoordType>::setTriMeshSpatialModel(TriMeshSpatialModel<CoordType>& triMeshSpatialModel)
//{
//  _triMeshSpatialModel = &triMeshSpatialModel;
//}

 /*! Sets the triMesh to work with.
****************************************************************/
template<class CoordType>
void MaximalRepulsionTriMeshSpatialModel<CoordType>::setTriMesh(const TriMesh<CoordType>& triMesh)
{
  _triMeshQuery.setTriMesh( triMesh );
  _nucleusVolume = abs(_triMeshQuery.getTriMesh().volume());
}

/*! Gets the current energy of the system
 * Minimizing the inverse of the sum of all the interdistances
****************************************************************/
template<class CoordType>
CoordType MaximalRepulsionTriMeshSpatialModel<CoordType>::getEnergy1(const Vertices<CoordType>& vertices)
{
  CoordType distBetweenComp;
  int numPoints = vertices.getNumVertices();

  for ( int i = 0; i < numPoints; ++i )
    for ( int j = 0; j < numPoints; ++j )
      if ( j!= i)
          distBetweenComp += vertices[i].distance(vertices[j]);

  return 1/distBetweenComp;
  //return numPoints/distBetweenComp;
}

/*! Gets the current energy of the system
 * Minimizing the inverse of the sum of the distances between closest objects
****************************************************************/
template<class CoordType>
CoordType MaximalRepulsionTriMeshSpatialModel<CoordType>::getEnergy2(const Vertices<CoordType>& vertices)
{
  int numPoints = vertices.getNumVertices();
  Vector<CoordType> distancesCompartment(numPoints), minDistances(numPoints);

  minDistances = vertices.squareNearestNeighborDistances();
//  for ( int i = 0; i < numPoints; ++i )
//  {
//    for ( int j = 0; j < numPoints; ++j )
//    {
//      if ( j!= i)
//        distancesCompartment[j] = vertices[i].distance(vertices[j]);
//      else
//        distancesCompartment[j] = numeric_limits<float>::max();
//    }
////    EVAL(distancesCompartment);
//    minDistances[i] = distancesCompartment.min();
//  }

  EVAL(minDistances);

  return 1/minDistances.sum();

}


/*! Gets the current energy of the system
****************************************************************/
template<class CoordType>
CoordType MaximalRepulsionTriMeshSpatialModel<CoordType>::getEnergy3(const Vertices<CoordType>& vertices)
{
  CoordType energy;
  int numPoints = vertices.getNumVertices();
  CoordType beta = _hardcoreDistances.mean()/(_nucleusVolume/numPoints);

  for ( int i = 0; i < numPoints; ++i )
    for ( int j = 0; j < numPoints; ++j )
      if ( j!= i)
        energy += exp(-(pow(vertices[i].distance(vertices[j]),2))/(2*pow(beta,2)));

  return energy;
}

/*! Gets the current energy of the system
****************************************************************/
template<class CoordType>
CoordType MaximalRepulsionTriMeshSpatialModel<CoordType>::getEnergy4(const Vertices<CoordType>& vertices)
{
  CoordType partialEnergy;
  int numPoints = vertices.getNumVertices();
  CoordType beta = _hardcoreDistances.mean()/(_nucleusVolume/numPoints);

  for ( int i = 0; i < numPoints; ++i )
    for ( int j = 0; j < numPoints; ++j )
      if ( j!= i)
        partialEnergy += pow(vertices[i].distance(vertices[j]),2);

  return exp(-(partialEnergy/(2*pow(beta,2))));
}

/*! Moves a compartment inside the container
****************************************************************/
template<class CoordType>
Vector<CoordType> MaximalRepulsionTriMeshSpatialModel<CoordType>::moveCompartment(const int numCompartment, const Vertices<CoordType>& vertices, const CoordType moveLimit)
//Vector<float> MaximalRepulsionTriMeshSpatialModel<CoordType>::moveCompartment(const int numCompartment, const Vertices<CoordType>& vertices, const CoordType moveLimit)
{
  RandomGenerator& randomGenerator = this->getRandomGenerator();
  Vector<CoordType> movedVertex(3);
  CoordType distX, distY, distZ;
  bool checkCondition = false;

  //moves the compartment in a random direction a distance smaller than twice the current distance to the border
  do
  {
    distX = randomGenerator.uniformLF(-moveLimit,moveLimit);//(globalCycles+1));
    distY = randomGenerator.uniformLF(-moveLimit,moveLimit);//(globalCycles+1));
    distZ = randomGenerator.uniformLF(-moveLimit,moveLimit);//(globalCycles+1));

    movedVertex[X] = vertices[numCompartment][X]+distX;
    movedVertex[Y] = vertices[numCompartment][Y]+distY;
    movedVertex[Z] = vertices[numCompartment][Z]+distZ;

    //checks that the new point's location respects the real volume and that is located inside the nucleus
    //if ( _triMeshQuery.contains(movedVertex) == true )
      if ( checkHardcoreDistances( movedVertex, numCompartment, vertices) == true )
        checkCondition = true;

  } while ( checkCondition == false );

  return movedVertex;
}

 /*! Checks if the hardcore distances between a new (moved) vertex and the others are respected.
****************************************************************/
template<class CoordType>
bool MaximalRepulsionTriMeshSpatialModel<CoordType>::checkHardcoreDistances(const Vector<CoordType>& vertex, const int numVertexChanged, const Vertices<CoordType>& vertices) const
{
  bool checkDistance;
  //const int k = vertices.getNumVertices();
  Vector<float> triMeshVertex;

  checkDistance = _triMeshQuery.closestPoint(vertex, triMeshVertex) > _hardcoreDistances[numVertexChanged];

  for (int i = 0; checkDistance == true && i < vertices.getNumVertices(); ++i)
    if ( i!= numVertexChanged)
        if ( vertices[i].distance( vertex ) < _hardcoreDistances[i] + _hardcoreDistances[numVertexChanged] )
          checkDistance = false;

  return checkDistance;
}

// /*! Processes the maxima repulsion
//****************************************************************/
//template<class CoordType>
//void MaximalRepulsionTriMeshSpatialModel<CoordType>::setMaximalRepulsion()
//{


//}

 /*! Find beta for this nucleus (used with the energy)
****************************************************************/
template<class CoordType>
CoordType MaximalRepulsionTriMeshSpatialModel<CoordType>::findBeta()
{
  EVAL(_vertices.getSize());
  EVAL(_numCompartments);
  EVAL(_hardcoreDistances.getSize());
  CoordType beta, epsilon = 1.0;
  CoordType temp = getEnergy(_vertices, 1)*(1/((_numCompartments-1)*(2*abs(_triMeshQuery.getTriMesh().equivalentRadius())-2*_hardcoreDistances.mean())));
  for (  beta = 0; beta < 1000000 && epsilon > 0.001; beta += 100 )
    epsilon = exp(-beta*temp);

  EVAL(epsilon)
  EVAL(beta);
  return beta;
}

/*! Processes the maxima repulsion
****************************************************************/
template<class CoordType>
Vertices<CoordType> MaximalRepulsionTriMeshSpatialModel<CoordType>::drawSample(const int numPoints)
{
  ENTER("Vertices<CoordType> MaximalRepulsionTriMeshSpatialModel<CoordType>::drawSample(const int)");

  bool info = true;
  RandomGenerator& randomGenerator = this->getRandomGenerator();
  Vector<CoordType> triMeshVertex(3), movedVertex(3);

  //calls the function to generate random points with volume inside the nucleus trimesh
  _numCompartments = numPoints;
  _vertices = hardcoreDistances();

  Vertices<CoordType> currentVertices = _vertices, minEnergyVertices = _vertices;
  Vector<CoordType> distancesToBorder(numPoints), numberMovements(numPoints), numberAttempts(numPoints);

  int i, globalMovements = 0;
  int numberLoopAttempts = 250;
  int maxTotalLoopAttempts = 200;

  DataSet dataset, dataset2;

  CoordType beta = findBeta();

  //CurveStack<CoordType> curveStack(3,_numCompartments,0,0,0,true);

  //sets up the minimum energy (over the initial real)
  Vector<CoordType> minEnergy(2);
  minEnergy[0]=1000;
  minEnergy[1]=1;

  //creates a vector of distances to the border
  for ( int c = 0; c < _numCompartments; ++c)
  {
    distancesToBorder[c] = _triMeshQuery.closestPoint( _vertices[c], triMeshVertex );
    EVAL(_vertices[c]);
  }

  //two loops
  for ( int globalCycles = 0; globalCycles < maxTotalLoopAttempts && numberMovements.sum() !=0 && globalMovements-minEnergy[1]<100; ++globalCycles )
  {
    numberMovements.setZeros();
    numberAttempts.setZeros();
    int numMovements = 0;

//    curveStack.setOpen(true);

    for ( int oneCycle = 0; oneCycle < numberLoopAttempts && globalMovements-minEnergy[1]<100 ; ++oneCycle )
    {
      //chooses a random compartment to move
      i = randomGenerator.uniformL( numPoints );

      const Vector<CoordType> vertex = currentVertices[i];
      //float currentDistanceToBorder = distancesToBorder[i];

//      //calculates old energy ------- using 1st method -- all interdistances
      float oldEnergy = getEnergy(currentVertices, 2);
      //calculates new energy ------- using 2nd method -- distance to the closest one
      //float oldEnergy = getEnergy(currentVertices, 1);

      //moves the compartment in a random direction a distance smaller than twice the current distance to the border
      //CoordType radius = 3*currentDistanceToBorder;
      //CoordType radius = (distancesToBorder[i]-_hardcoreDistances[i])*.99;
      CoordType radius = (_hardcoreDistances[i])*.1;
      movedVertex = moveCompartment(i, currentVertices, radius);
      currentVertices[i] = movedVertex;

      ++numberAttempts[i];
//      //calculates new energy ------- using 1st method -- all interdistances
      float newEnergy = getEnergy(currentVertices, 2);
      //calculates new energy ------- using 2nd method -- distance to the closest one
//      float newEnergy = getEnergy(currentVertices, 1);

      bool secondChance = false;
      //gives "an opportunity" to accept a "bad" energy change depending on its probability
      if ( newEnergy-oldEnergy > 0 )
      {
        float tempEnergy = exp(-beta*(newEnergy-oldEnergy));
        float coef = randomGenerator.uniformLF(0,1);

        if ( coef < tempEnergy ) secondChance = true;
//        EVAL(coef);
//        EVAL(tempEnergy);
//        EVAL(newEnergy);
//        EVAL(oldEnergy );
      }

      //checks whether the object is within the confined space or not
      //**to remember: distances to other objects are not tested because they do are in getEnergy
      if ( ( _triMeshQuery.contains(movedVertex) ) )
      {
        float tempDistanceToTheBorder = _triMeshQuery.closestPoint( movedVertex, triMeshVertex );

        if  ( tempDistanceToTheBorder > _hardcoreDistances[i] )
        {
          //check if the conditions and therefore the new position are accepted
          if ( (newEnergy-oldEnergy < 0) || secondChance == true )
          {
            distancesToBorder[i] = tempDistanceToTheBorder;

            if (info == true )
            {
              dataset.setValue("probs",numMovements,exp(-(newEnergy-oldEnergy)));
              dataset.setValue("deltaEnergy",numMovements,newEnergy-oldEnergy);
              dataset.setValue("newEnergy",numMovements,newEnergy);
              dataset.setValue("numMovements",numMovements,numMovements);

              dataset2.setValue("numMovements",globalMovements,globalMovements);
              dataset2.setValue("deltaEnergy",globalMovements,newEnergy-oldEnergy);
              dataset2.setValue("newEnergy",globalMovements,newEnergy);
            }

            if ( newEnergy < minEnergy[0] )
            {
              minEnergy[0] = newEnergy;
              minEnergy[1] = globalMovements;
              minEnergyVertices = currentVertices;
              EVAL(minEnergy[0]);
              EVAL(minEnergy[1]);
            }

    //        EVAL(newEnergy);
            ++numMovements;
            ++numberMovements[i];
            ++globalMovements;
          }
        }
      }
      else currentVertices[i] = vertex; //if the movement is not accepted we return to the previous position
      //EVAL(numberMovements);

//      for ( int c = 0; c < _numCompartments; ++c)
//        curveStack[c].append(currentVertices[c]);

    }

//    if (info == true)
//    {
//      ostringstream iss;
//      iss << globalCycles;
//      currentVertices.save( "/home/jarpon/Desktop/new/maxRepulsion" + iss.str() + ".vx", true );
//      dataset.save("/home/jarpon/Desktop/new/data" + iss.str() + ".csv",true);
//      dataset2.save("/home/jarpon/Desktop/new/dataEnergy.csv",true);
//    }

//    EVAL(globalMovements);
//    EVAL(numberMovements);
//    EVAL(numberAttempts);

  }

//  EVAL(curveStack.getNumVertices());

//  EVAL(minEnergy[0]);
//  EVAL(globalMovements);
  for (int cc = 0; cc < _numCompartments; ++cc)
  {
    EVAL(minEnergyVertices[cc]);
  }
  LEAVE();
  return minEnergyVertices;
}

template<class CoordType>
CoordType MaximalRepulsionTriMeshSpatialModel<CoordType>::getEnergy(const Vertices<CoordType>& vertices, const int& method)
{
  switch( method )
  {
    case 1:
      return getEnergy1( vertices );
      break;
    case 2:
      return getEnergy2( vertices );
      break;
    case 3:
      return getEnergy3( vertices );
      break;
    case 4:
      return getEnergy4( vertices );
      break;
  }
}

//template class MaximalRepulsionTriMeshSpatialModel<double>;
template class MaximalRepulsionTriMeshSpatialModel<float>;
//template class MaximalRepulsionTriMeshSpatialModel<int>;
