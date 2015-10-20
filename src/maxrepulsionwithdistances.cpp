//#include "trimeshspatialmodel.h"
//#include "maxrepulsionwithdistances.h"

//#include <randomgenerator.h>
//#include <trimesh.h>
//#include <alinevector.h>
//#include <cmath>

//#include <trimeshquery.h>

//#include <programerror.h>
//#include <dataset.h>

//#include <curvestack.h>
//#include <vertexstack.h>

////#define TRACE
//#include <trace.h>
//#include <sstream>

//using namespace std;

///*! Initializes all the parameters.
//****************************************************************/
//template<class CoordType>
//MaxRepulsionWithDistancesTriMeshSpatialModel<CoordType>::MaxRepulsionWithDistancesTriMeshSpatialModel()
//{
//}

///*! Destroys it.
//****************************************************************/
//template<class CoordType>
//MaxRepulsionWithDistancesTriMeshSpatialModel<CoordType>::~MaxRepulsionWithDistancesTriMeshSpatialModel()
//{
//}

///*! Sets a vector of hardcore distances.
//****************************************************************/
//template<class CoordType>
//void MaxRepulsionWithDistancesTriMeshSpatialModel<CoordType>::setHardcoreDistance(const Vector<CoordType>& hardcoreDistances)
//{
//  _hardcoreDistances = hardcoreDistances;
//}

//// /*! Sets the triMeshSpatialModel to work with ------ (not needed)
////****************************************************************/
////template<class CoordType>
////void MaxRepulsionWithDistancesTriMeshSpatialModel<CoordType>::setTriMeshSpatialModel(TriMeshSpatialModel<CoordType>& triMeshSpatialModel)
////{
////  _triMeshSpatialModel = &triMeshSpatialModel;
////}

// /*! Sets the triMesh to work with.
//****************************************************************/
//template<class CoordType>
//void MaxRepulsionWithDistancesTriMeshSpatialModel<CoordType>::setTriMesh(const TriMesh<CoordType>& triMesh)
//{
//  _triMeshQuery.setTriMesh( triMesh );
//  _nucleusVolume = abs(_triMeshQuery.getTriMesh().volume());
//}

///*! Gets the current energy of the system
//****************************************************************/
//template<class CoordType>
//CoordType MaxRepulsionWithDistancesTriMeshSpatialModel<CoordType>::getEnergy1(const Vertices<CoordType>& vertices)
//{
//  CoordType distBetweenComp;
//  int numPoints = vertices.getNumVertices();

//  for ( int i = 0; i < numPoints; ++i )
//    for ( int j = 0; j < numPoints; ++j )
//      if ( j!= i)
//          distBetweenComp += vertices[i].distance(vertices[j]);

//  return 1/distBetweenComp;
//  //return numPoints/distBetweenComp;
//}

///*! Gets the current energy of the system
//****************************************************************/
//template<class CoordType>
//CoordType MaxRepulsionWithDistancesTriMeshSpatialModel<CoordType>::getEnergy2(const Vertices<CoordType>& vertices)
//{
//  int numPoints = vertices.getNumVertices();
//  Vector<CoordType> distancesCompartment(numPoints), minDistances(numPoints);

//  for ( int i = 0; i < numPoints; ++i )
//  {
//    for ( int j = 0; j < numPoints; ++j )
//      if ( j!= i)
//        distancesCompartment[j] = vertices[i].distance(vertices[j]);
//      else
//        distancesCompartment[j] = 10000;
//    minDistances[i] = distancesCompartment.min();
//  }

//  return 1/minDistances.sum();
//}


///*! Gets the current energy of the system
//****************************************************************/
//template<class CoordType>
//CoordType MaxRepulsionWithDistancesTriMeshSpatialModel<CoordType>::getEnergy3(const Vertices<CoordType>& vertices)
//{
//  CoordType energy;
//  int numPoints = vertices.getNumVertices();
//  CoordType beta = _hardcoreDistances.mean()/(_nucleusVolume/numPoints);

//  for ( int i = 0; i < numPoints; ++i )
//    for ( int j = 0; j < numPoints; ++j )
//      if ( j!= i)
//        energy += exp(-(pow(vertices[i].distance(vertices[j]),2))/(2*pow(beta,2)));

//  return energy;
//}

///*! Gets the current energy of the system
//****************************************************************/
//template<class CoordType>
//CoordType MaxRepulsionWithDistancesTriMeshSpatialModel<CoordType>::getEnergy4(const Vertices<CoordType>& vertices)
//{
//  CoordType partialEnergy;
//  int numPoints = vertices.getNumVertices();
//  CoordType beta = _hardcoreDistances.mean()/(_nucleusVolume/numPoints);

//  for ( int i = 0; i < numPoints; ++i )
//    for ( int j = 0; j < numPoints; ++j )
//      if ( j!= i)
//        partialEnergy += pow(vertices[i].distance(vertices[j]),2);

//  return exp(-(partialEnergy/(2*pow(beta,2))));
//}

///*! Moves a compartment inside the container
//****************************************************************/
//template<class CoordType>
//Vector<CoordType> MaxRepulsionWithDistancesTriMeshSpatialModel<CoordType>::moveCompartment(const int numCompartment, const Vertices<CoordType>& vertices, const CoordType moveLimit)
////Vector<float> MaxRepulsionWithDistancesTriMeshSpatialModel<CoordType>::moveCompartment(const int numCompartment, const Vertices<CoordType>& vertices, const CoordType moveLimit)
//{
//  RandomGenerator& randomGenerator = this->getRandomGenerator();
//  Vector<CoordType> movedVertex(3);
//  CoordType distX, distY, distZ;
//  bool checkCondition = false;

//  //moves the compartment in a random direction a distance smaller than twice the current distance to the border
//  do
//  {
//    distX = randomGenerator.uniformLF(-moveLimit,moveLimit);//(globalCycles+1));
//    distY = randomGenerator.uniformLF(-moveLimit,moveLimit);//(globalCycles+1));
//    distZ = randomGenerator.uniformLF(-moveLimit,moveLimit);//(globalCycles+1));

//    movedVertex[X] = vertices[numCompartment][X]+distX;
//    movedVertex[Y] = vertices[numCompartment][Y]+distY;
//    movedVertex[Z] = vertices[numCompartment][Z]+distZ;

//    //checks that the new point's location respects the real volume and that is located inside the nucleus
//    //if ( _triMeshQuery.contains(movedVertex) == true )
//      if ( checkHardcoreDistances( movedVertex, numCompartment, vertices) == true )
//        checkCondition = true;

//  } while ( checkCondition == false );

//  return movedVertex;
//}


///*! Moves a compartment inside the container setting/checking the distance to the border
//****************************************************************/
//template<class CoordType>
//Vector<CoordType> MaxRepulsionWithDistancesTriMeshSpatialModel<CoordType>::moveCompartmentRespectingDistance(const int numCompartment, const Vertices<CoordType>& vertices, const CoordType moveLimit)
////Vector<float> MaxRepulsionWithDistancesTriMeshSpatialModel<CoordType>::moveCompartment(const int numCompartment, const Vertices<CoordType>& vertices, const CoordType moveLimit)
//{
//  RandomGenerator& randomGenerator = this->getRandomGenerator();
//  Vector<CoordType> movedVertex(3);
//  CoordType distX, distY, distZ;
//  bool checkCondition, found = false;
//  Vector<CoordType>  triMeshVertex, normalVector;
//  int numAttempts = 0;

//  //moves the compartment in a random direction a distance smaller than twice the current distance to the border
//  do
//  {
//    distX = randomGenerator.uniformLF(-moveLimit,moveLimit);//(globalCycles+1));
//    distY = randomGenerator.uniformLF(-moveLimit,moveLimit);//(globalCycles+1));
//    distZ = randomGenerator.uniformLF(-moveLimit,moveLimit);//(globalCycles+1));

//    movedVertex[X] = vertices[numCompartment][X]+distX;
//    movedVertex[Y] = vertices[numCompartment][Y]+distY;
//    movedVertex[Z] = vertices[numCompartment][Z]+distZ;

//    _triMeshQuery.closestPoint( movedVertex, triMeshVertex );
//    normalVector = movedVertex - triMeshVertex;
//    normalVector /= movedVertex.distance( triMeshVertex );
//    normalVector *= _distanceToBorder[numCompartment];
//    movedVertex = triMeshVertex + normalVector;

//    //checks that the new point's location respects the real volume and that is located inside the nucleus
//    //if ( _triMeshQuery.contains(movedVertex) == true )
//    if ( _triMeshQuery.contains(movedVertex) )
//      found = fabs( _distanceToBorder[numCompartment] - _triMeshQuery.closestPoint( movedVertex, triMeshVertex ) ) < movedVertex.epsilon();

//    if ( found == true && checkHardcoreDistances( movedVertex, numCompartment, vertices) == true )
//        checkCondition = true;

//  } while ( ( ++numAttempts < 50) || ( checkCondition == false ) );

//  return movedVertex;
//}


// /*! Checks if the hardcore distances between a new (moved) vertex and the others are respected.
//****************************************************************/
//template<class CoordType>
//bool MaxRepulsionWithDistancesTriMeshSpatialModel<CoordType>::checkHardcoreDistances(const Vector<CoordType>& vertex, const int numVertexChanged, const Vertices<CoordType>& vertices) const
//{
//  bool checkDistance;
//  //const int k = vertices.getNumVertices();
//  Vector<float> triMeshVertex;

//  checkDistance = _triMeshQuery.closestPoint(vertex, triMeshVertex) > _hardcoreDistances[numVertexChanged];

//  for (int i = 0; checkDistance == true && i < vertices.getNumVertices(); ++i)
//    if ( i!= numVertexChanged)
//        if ( vertices[i].distance( vertex ) < _hardcoreDistances[i] + _hardcoreDistances[numVertexChanged] )
//          checkDistance = false;

//  return checkDistance;
//}

//// /*! Processes the maxima repulsion
////****************************************************************/
////template<class CoordType>
////void MaxRepulsionWithDistancesTriMeshSpatialModel<CoordType>::setMaximalRepulsion()
////{


////}

// /*! Find beta for this nucleus (used with the energy)
//****************************************************************/
//template<class CoordType>
//CoordType MaxRepulsionWithDistancesTriMeshSpatialModel<CoordType>::findBeta()
//{
//  EVAL(_vertices.getSize());
//  EVAL(_numCompartments);
//  EVAL(_hardcoreDistances.getSize());
//  CoordType beta, epsilon = 1.0;
//  CoordType temp = getEnergy(_vertices, 1)*(1/((_numCompartments-1)*(2*abs(_triMeshQuery.getTriMesh().equivalentRadius())-2*_hardcoreDistances.mean())));
//  for (  beta = 0; beta < 1000000 && epsilon > 0.001; beta += 100 )
//    epsilon = exp(-beta*temp);

//  EVAL(epsilon)
//  EVAL(beta);
//  return beta;
//}

///*! Sets a vector of distances to the border.
//****************************************************************/
//template<class CoordType>
//void MaxRepulsionWithDistancesTriMeshSpatialModel<CoordType>::setDistanceToBorder(const Vector<CoordType>& distancesToTheBorder)
//{
//  _distanceToBorder = distancesToTheBorder;
//}

///*! Checks if the distance between a vertex and its closest point of the triMesh is correct.
//****************************************************************/
//template<class CoordType>
//bool MaxRepulsionWithDistancesTriMeshSpatialModel<CoordType>::checkDistancesToBorder(
//  const Vector<CoordType>& vertex,
//  CoordType position)
//{
//  Vector<CoordType> triMeshVertex = vertex;
//  //const CoordType distance = _triMeshQuery.getTriMesh().closestPoint( vertex, triMeshVertex );
//  return abs(_triMeshQuery.getTriMesh().closestPoint( vertex, triMeshVertex )-_distanceToBorder[position]) < vertex.epsilon();
//}

///*! Generates vertices into the trimesh with a fixed hardcore distance between compartments
// * and fixed distances to the border.
//****************************************************************/
//template<class CoordType>
//Vertices<CoordType> MaxRepulsionWithDistancesTriMeshSpatialModel<CoordType>::hardcoreAndToTheBorderDistances()
//{
//  //ENTER("void TriMeshSpatialModel<CoordType>::hardcoreAndToTheBorderDistances()");
//  Vertices<CoordType> vertices( 3, 0, 0, 0 );
//  Vector<CoordType> vertex(3);
//  RandomGenerator& randomGenerator = this->getRandomGenerator();

////  Stopwatch stopWatch;
////  stopWatch.start( "One pattern" );

//  int attempts = 0;
//  const int maxAttempts = 200;

//  if ( _numCompartments != 1 && _hardcoreDistances.getSize() == 1 && _distanceToBorder.getSize() == 1 )
//  {
//    float tempHD = _hardcoreDistances[0];
//    float tempDB = _distanceToBorder[0];
//    _hardcoreDistances.setSize( _numCompartments );
//    _distanceToBorder.setSize( _numCompartments );
//    _hardcoreDistances.fill( tempHD );
//    _distanceToBorder.fill( tempDB );
//  }
//  else if (  ( _numCompartments != _hardcoreDistances.getSize() && _hardcoreDistances.getSize() != 1 ) ||
//             ( _numCompartments != _distanceToBorder.getSize() && _distanceToBorder.getSize() != 1 ) )
//  {
//    ProgramError programError;
//    programError.setWhere( "void TriMeshSpatialModel<CoordType>::hardcoreAndToTheBorderDistances()" );
//    programError.setWhat( "The number of compartments and vectors' lengths of hardcore and to the border distances are not the same" );
//    throw programError;
//    return vertices;
//    //LEAVE();
//  }

//  //randomizes the order of the compartments
//  randomizesOrder(randomGenerator);
//  EVAL(_hardcoreDistances);

//  for (int i = 0; i < _numCompartments; ++i)
//  {
//    attempts = 0;

//    do
//    {
//      //checks if both conditions are true at the same time
//      drawPositionFromBorder( vertex , _distanceToBorder[i] );

//    } while ( checkInterObjectDistances(vertex,vertices) == false && ++attempts < maxAttempts );

//    vertices.append( vertex );

//    EVAL(attempts);
//    if ( attempts > maxAttempts )
//    {
//      //EVAL(attempts);
//      Exception exception;
//      exception.setWhere( "void TriMeshSpatialModel<CoordType>::hardcoreAndToTheBorderDistances()" );
//      exception.setWhat( "The number of attempts trying to generate a wanted point was too high" );
//      throw exception;
//    }

//  }

////  stopWatch.stop( "Fonction vertices" );
////  stopWatch.print();
//  //LEAVE();
//  return vertices;
//}

/////*! Checks if the hardcore distances from a vertex to previous generated vertices are correct.
////****************************************************************/
////template<class CoordType>
////bool MaxRepulsionWithDistancesTriMeshSpatialModel<CoordType>::checkInterObjectDistances(const Vector<CoordType>& vertex, const Vertices<CoordType>& vertices)
////{
////  const int k = vertices.getNumVertices();

////  for (int i = 0; i < vertices.getNumVertices(); ++i)
////    if ( vertices[i].distance( vertex ) < (_hardcoreDistances[i] + _hardcoreDistances[k]) )
////      return false;

////  return true;
////}

///*! Processes the maxima repulsion
//****************************************************************/
//template<class CoordType>
//Vertices<CoordType> MaxRepulsionWithDistancesTriMeshSpatialModel<CoordType>::drawSample(const int numPoints)
//{
//  ENTER("Vertices<CoordType> MaxRepulsionWithDistancesTriMeshSpatialModel<CoordType>::drawSample(const int)");

//  bool info = true;
//  RandomGenerator& randomGenerator = this->getRandomGenerator();
//  Vector<CoordType> triMeshVertex(3), movedVertex(3);

//  //calls the function to generate random points with volume inside the nucleus trimesh
//  _numCompartments = numPoints;
//  _vertices = hardcoreAndToTheBorderDistances();

//  Vertices<CoordType> currentVertices = _vertices, minEnergyVertices = _vertices;
//  //Vector<CoordType> distancesToBorder(numPoints), numberMovements(numPoints), numberAttempts(numPoints);
//  Vector<CoordType> numberMovements(numPoints), numberAttempts(numPoints);

//  int i, globalMovements = 0;
//  int numberLoopAttempts = 20;
//  int maxTotalLoopAttempts = 200;

//  DataSet dataset, dataset2;

//  CoordType beta = findBeta();
//  bool secondChance = false;

//  CurveStack<CoordType> curveStack(3,_numCompartments,0,0,0,true);

//  //sets up the minimum energy (over the initial real)
//  Vector<CoordType> minEnergy(2);
//  minEnergy[0]=1000;
//  minEnergy[1]=1;

////  //creates a vector of distances to the border
////  for ( int c = 0; c < _numCompartments; ++c)
////  {
////    distancesToBorder[c] = _triMeshQuery.closestPoint( _vertices[c], triMeshVertex );
////    EVAL(_vertices[c]);
////  }

//  //two loops
//  for ( int globalCycles = 0; globalCycles < maxTotalLoopAttempts && numberMovements.sum() !=0 && globalMovements-minEnergy[1]<100; ++globalCycles )
//  {
//    numberMovements.setZeros();
//    numberAttempts.setZeros();
//    int numMovements = 0;

////    curveStack.setOpen(true);

//    for ( int oneCycle = 0; oneCycle < numberLoopAttempts && globalMovements-minEnergy[1]<100 ; ++oneCycle )
//    {
//      //chooses a random compartment to move
//      i = randomGenerator.uniformL( numPoints );

//      const Vector<CoordType> vertex = currentVertices[i];
//      //float currentDistanceToBorder = distancesToBorder[i];

////      //calculates old energy ------- using 1st method -- all interdistances
////      float oldEnergy = getEnergy(currentVertices, 2);
//      //calculates new energy ------- using 2nd method -- distance to the closest one
//      float oldEnergy = getEnergy(currentVertices, 2);

//      //moves the compartment in a random direction a distance smaller than twice the current distance to the border
//      //CoordType radius = 3*currentDistanceToBorder;
//      CoordType radius = ( _distanceToBorder[i]-_hardcoreDistances[i] ) * 2;
////      movedVertex = moveCompartment(i, currentVertices, radius);
//      movedVertex = moveCompartmentRespectingDistance(i, currentVertices, radius);
//      currentVertices[i] = movedVertex;

//      ++numberAttempts[i];
////      //calculates new energy ------- using 1st method -- all interdistances
////      float newEnergy = getEnergy(currentVertices, 2);
//      //calculates new energy ------- using 2nd method -- distance to the closest one
//      float newEnergy = getEnergy(currentVertices, 2);

//      secondChance = false;
//      //gives "an opportunity" to accept a "bad" energy change depending on its probability
//      if ( newEnergy-oldEnergy > 0 )
//      {
//        float tempEnergy = exp(-beta*(newEnergy-oldEnergy));
//        float coef = randomGenerator.uniformLF(0,1);

//        if ( coef < tempEnergy ) secondChance = true;
//      }

//      //check if the conditions and therefore the new position are accepted
//      if ( ( (newEnergy-oldEnergy < 0) || ( secondChance == true ) ) )
//      {
//        //distancesToBorder[i] = _triMeshQuery.closestPoint( movedVertex, triMeshVertex );

//        if (info == true )
//        {
//          dataset.setValue("probs",numMovements,exp(-(newEnergy-oldEnergy)));
//          dataset.setValue("deltaEnergy",numMovements,newEnergy-oldEnergy);
//          dataset.setValue("newEnergy",numMovements,newEnergy);
//          dataset.setValue("numMovements",numMovements,numMovements);

//          dataset2.setValue("numMovements",globalMovements,globalMovements);
//          dataset2.setValue("deltaEnergy",globalMovements,newEnergy-oldEnergy);
//          dataset2.setValue("newEnergy",globalMovements,newEnergy);
//        }

//        if ( newEnergy < minEnergy[0] )
//        {
//          minEnergy[0] = newEnergy;
//          minEnergy[1] = globalMovements;
//          minEnergyVertices = currentVertices;
//          EVAL(minEnergy[0]);
//          EVAL(minEnergy[1]);
//        }



////        EVAL(newEnergy);
//        ++numMovements;
//        ++numberMovements[i];
//        ++globalMovements;
//      }
//      else currentVertices[i] = vertex; //if the movement is not accepted we return to the previous position
//      //EVAL(numberMovements);


//      for ( int c = 0; c < _numCompartments; ++c)
//        curveStack[c].append(currentVertices[c]);


//    }

//    curveStack.save( "/home/jarpon/Desktop/new/trace", true );

//    if (info == true)
//    {
//      ostringstream iss;
//      iss << globalCycles;
//      currentVertices.save( "/home/jarpon/Desktop/max+dist/maxRepulsion" + iss.str() + ".vx", true );
//      dataset.save("/home/jarpon/Desktop/max+dist/data" + iss.str() + ".csv",true);
//      dataset2.save("/home/jarpon/Desktop/max+dist/dataEnergy.csv",true);
//    }

////    EVAL(globalMovements);
////    EVAL(numberMovements);
////    EVAL(numberAttempts);

//  }

////  EVAL(curveStack.getNumVertices());

////  if (info == true) minEnergyVertices.save( "/home/jarpon/Desktop/new/minEnergyVertices.vx", true );
////  EVAL(minEnergy[0]);
////  EVAL(globalMovements);
//  for (int cc = 0; cc < _numCompartments; ++cc)
//  {
//    EVAL(minEnergyVertices[cc]);
//  }
//  LEAVE();
//  return minEnergyVertices;
//}

//template<class CoordType>
//CoordType MaxRepulsionWithDistancesTriMeshSpatialModel<CoordType>::getEnergy(const Vertices<CoordType>& vertices, const int& method)
//{
//  switch( method )
//  {
//    case 1:
//      return getEnergy1( vertices );
//      break;
//    case 2:
//      return getEnergy2( vertices );
//      break;
//    case 3:
//      return getEnergy3( vertices );
//      break;
//    case 4:
//      return getEnergy4( vertices );
//      break;
//  }
//}

////template class MaxRepulsionWithDistancesTriMeshSpatialModel<double>;
//template class MaxRepulsionWithDistancesTriMeshSpatialModel<float>;
////template class MaxRepulsionWithDistancesTriMeshSpatialModel<int>;
