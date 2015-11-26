///*!
// * \class  SpatialModelHardcoreBorderDistance3D
// * \author Javier Arpón (ja), INRA
// * \author Philippe Andrey (pa), INRA
// * \date   XXXX.XX.XX - creation (ja)
// * \date   2015.10.12 - integration (pa)
// * \brief  3D hardcore point process with distance-to-border constraint
//****************************************************************/

//#include "spatialmodelhardcoreborderdistance3ddifferentcompartments.h"

//#include <programerror.h>

///*! Constructor.
//****************************************************************/
//template<class CoordType>
//SpatialModelHardcoreBorderDistance3DDifferentCompartments<CoordType>::SpatialModelHardcoreBorderDistance3DDifferentCompartments() :
//  SpatialModelBorderDistance3DDifferentCompartments<CoordType>(),
//  SpatialModelHardcoreDistance3DDifferentCompartments<CoordType>()
//{
//}

///*! Generates vertices into the trimesh with a fixed hardcore distance between compartments
// * and fixed distances to the border.
//****************************************************************/
//template<class CoordType>
//Vertices<CoordType> SpatialModelHardcoreBorderDistance3DDifferentCompartments<CoordType>::drawSample(const int numVertices, const int numVertices2)
//{
//  const Vector<CoordType>& hardcoreDistances = this->getHardcoreDistances();
//  const Vector<CoordType>& distancesToBorder = this->getDistancesToBorder();

//  if ( numVertices != hardcoreDistances.getSize() || numVertices != distancesToBorder.getSize() )
//  {
//    ProgramError programError;
//    programError.setWhere( "void SpatialModelHardcoreBorderDistance3D<CoordType>::drawSample(const int)" );
//    programError.setWhat( "The number of vertices differs from the number of distances (hc or border)" );
//    throw programError;
//  }

//  const int numPermutations = numVertices;
//  const int maxAttempts = 200;
//  Vector<CoordType> shuffledHardcoreDistances = hardcoreDistances;
//  Vector<CoordType> shuffledDistancesToBorder = distancesToBorder;
//  Vertices<CoordType> vertices( 3, 0 );
//  Vector<CoordType> vertex( 3 );
//  int attempts;
//  bool success = false;

//  for (int p = 0; p < numPermutations && !success; ++p)
//  {
//    shuffleDistances( shuffledHardcoreDistances, shuffledDistancesToBorder );
//    vertices.setSize( 0 );
//    success = true;

//    for (int v = 0; v < numVertices; ++v)
//    {
//      attempts = 0;

//      do {
//        this->drawPositionFromBorder( vertex, shuffledDistancesToBorder[v] );
//      } while ( !this->validInterObjectDistances(vertex,vertices,shuffledHardcoreDistances)
//                && ++attempts < maxAttempts );

//      if ( attempts < maxAttempts )
//        vertices.append( vertex );
//      else
//      {
//        success = false;
//        break;
//      }
//    }
//  }

//  if ( !success )
//  {
//    Exception exception;
//    exception.setWhere( "void SpatialModelHardcoreBorderDistance3D<CoordType>::drawSample(const int)" );
//    exception.setWhat( "Too many unsuccessful attemps to generate a vertex" );
//    throw exception;
//  }

//  return vertices;
//}

//template<class CoordType>
//void SpatialModelHardcoreBorderDistance3DDifferentCompartments<CoordType>::shuffleDistances(
//  Vector<CoordType>& hardcoreDistances,
//  Vector<CoordType>& distancesToBorder)
//{
//  const int n = hardcoreDistances.getSize();
//  int i, i1, i2;
//  CoordType tmp;

//  for (i = 0; i < n; ++i)
//  {
//    i1 = this->getRandomGenerator().uniformL( n );
//    i2 = this->getRandomGenerator().uniformL( n );
//    if ( i1 != i2 )
//    {
//      tmp = hardcoreDistances[i1]; hardcoreDistances[i1] = hardcoreDistances[i2]; hardcoreDistances[i2] = tmp;
//      tmp = distancesToBorder[i1]; distancesToBorder[i1] = distancesToBorder[i2]; distancesToBorder[i2] = tmp;
//    }
//  }
//}

//template class SpatialModelHardcoreBorderDistance3DDifferentCompartments<float>;
//template class SpatialModelHardcoreBorderDistance3DDifferentCompartments<double>;
//template class SpatialModelHardcoreBorderDistance3DDifferentCompartments<long double>;
