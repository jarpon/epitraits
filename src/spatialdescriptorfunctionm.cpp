///*!
// * \class  SpatialDescriptorFunctionM
// * \author Javier Arpon (ja), INRA
// * \author Thotino Gobin-Gansou (tg-g), INRA
// * \date   2015.10.15 - creation (ja)
// * \brief  M-function statistics for spatial point patterns
// * \details This class implements the distance-to-border cumulative
// * distribution. This is one of the available spatial descriptors
// * for point patterns.
// * The class handles the 2D and the 3D cases.
// * \brief  Study the peripheral positioning of volumes respect an envelope property
//****************************************************************/

//#include "spatialdescriptorfunctionm.h"

//#include "spatialdescriptor.h"

//#include <cdftools.h>
//#include <trimesh.h>
//#include <cmath>
//#include <dataset.h>
//#include <trimeshquery.h>
//#include <sstream>

//#define TRACE
//#include <trace.h>

//template<class CoordType>
//SpatialDescriptorFunctionM<CoordType>::SpatialDescriptorFunctionM() : SpatialDescriptor<CoordType>()
//{
//    //_triMesh = 0;
//}

//template<class CoordType>
//void SpatialDescriptorFunctionM<CoordType>::setTriMesh(const TriMesh<CoordType>& triMesh)
//{
//  _triMesh = &triMesh;
//  _externalProperties = false;
//}

//template<class CoordType>
//void SpatialDescriptorFunctionM<CoordType>::setPropertiesTriMesh(const Shape<CoordType>& properties )
//{
//  _propertiesMap = &properties;
//  _externalProperties = true;

//  // in this example the property is the intensity map of the envelope
////  TriMesh<CoordType> triMeshPropertyMap( _propertiesTriMesh->getSpace(), _propertiesTriMesh->getSize() );

////  for ( int i = 0; i < _propertiesTriMesh->getSize(); ++i )
////  {
//////    triMeshPropertyMap[i][X] = int(_propertiesTriMesh[i][X]);
//////    triMeshPropertyMap[i][Y] = int(_propertiesTriMesh[i][Y]);
//////    triMeshPropertyMap[i][Z] = int(_propertiesTriMesh[i][Z]);
////      triMeshPropertyMap[i][X] = _propertiesTriMesh[i][X];
////      triMeshPropertyMap[i][Y] = _propertiesTriMesh[i][Y];
////      triMeshPropertyMap[i][Z] = _propertiesTriMesh[i][Z];

////  }

////  for ( int i = 0; i < _propertiesTriMesh->getNumTriangles(); ++i )
////  {
////    //EVAL(_propertiesTriMesh.getTriangleIndexes(i));
////    triMeshPropertyMap.addTriangle(
////          _propertiesTriMesh->getTriangleIndexes(i)[0],
////          _propertiesTriMesh->getTriangleIndexes(i)[1],
////          _propertiesTriMesh->getTriangleIndexes(i)[2] );
////  }

//  //_triMesh.setVertexColors( new TriMesh<CoordType>( triMeshPropertyMap ) );
////  _triMesh->setVertexColors( new TriMesh<CoordType>( _propertiesTriMesh ) );
////  _triMesh->setVertexColors( _propertiesTriMesh );
//}

//template<class CoordType>
//void SpatialDescriptorFunctionM<CoordType>::eval(
//    const Vertices<CoordType>& vertices,
//    Vector<CoordType>& xvalues,
//    Vector<CoordType>& yvalues
//    )
//{
//  int triangleIndex;

//    ENTER( "void SpatialDescriptorFunctionM<CoordType>::eval(...)" );
//    const int numVertices = vertices.getSize();

//    Vector<CoordType> vertex(3);
//    Vector<CoordType> triMeshVertex(3);
//    xvalues.setSize( numVertices );
//    TriMeshQuery<CoordType> triMeshQuery;

//    triMeshQuery.setTriMesh( _triMesh );
//    for (int i = 0; i < numVertices; ++i)
//    {
//      vertex = vertices[i];
//      triMeshQuery.closestPoint( vertex, triMeshVertex, triangleIndex );
//      xvalues[i]= getValue( vertex, triangleIndex );
//    }
//    xvalues.sort();

//    CDFTools<CoordType> cdftools;
//    yvalues = cdftools.cdf( xvalues );

////    DataSet dataCFDTools;
////    for(int i=0; i< numVertices; ++i)
////    {
////      dataCFDTools.setValue<CoordType>("yvalues", i, yvalues[i]);
////      dataCFDTools.setValue<CoordType>("xvalues", i, xvalues[i]);
////    }
////    dataCFDTools.save(realPointDir + "/CumulativeFonctionDistribution.csv", true);

//    LEAVE();
//}

///*!
// * Get the intensity at position \par position in the triangle \par triangle.
// *
// * The intensity is interpolated from the ones at the triangle vertices.
// */
//template<class CoordType>
//CoordType SpatialDescriptorFunctionM<CoordType>::triangleInterpolation(
//    const Triangle<CoordType>& triangle,
//    const Vector<CoordType>& properties,
//    const Vector<CoordType>& position)
//{

//  return ( ( ( properties[0] * position.distance( triangle[0] ) ) + ( properties[1] * position.distance( triangle[1] ) ) + ( properties[2] * position.distance( triangle[2] ) ) )
//      / ( position.distance( triangle[0] ) + position.distance(triangle[1]) + position.distance(triangle[2]) ) );

//}

///*!
// *
// */
//template<class CoordType>
//CoordType SpatialDescriptorFunctionM<CoordType>::getValue(
//    const Vector<CoordType>& position,
//    const int& triangleIndex )
//{


//  Triangle<CoordType> triangle;
//  triangle = _triMesh->triangle( triangleIndex );
//  Vector<unsigned int> vertices(3);
//  vertices = _triMesh->getTriangleIndexes( triangleIndex );
//  Vector<CoordType> properties(3);
//  TriMesh<int> colortm( 3,0 );
//  colortm = _triMesh->getVertexColors();

//  const int dim = _propertiesMap->getSpace();

////  for(int i = 0; i < triMesh.getNumTriangles() ; i++)
////  {
////    colortm.addTriangle(
////          triMesh.getTriangleIndexes(i)[0],
////          triMesh.getTriangleIndexes(i)[1],
////          triMesh.getTriangleIndexes(i)[2] );
////  }

//  if ( _externalProperties == false )
//  {
////    properties[X] = _triMesh->triangle(triangleIndex)[0][X];
////    properties[Y] = _triMesh->triangle(triangleIndex)[1][X];
////    properties[Z] = _triMesh->triangle(triangleIndex)[2][X];
//      properties[X] = triangle[0][X];
//      properties[Y] = triangle[1][X];
//      properties[Z] = triangle[2][X];
//  }
//  else if ( dim == 3 )
//  {
//    properties[X] = _externalProperties(triangleIndex)[0][X];
//    properties[Y] = _externalProperties(triangleIndex)[1][X];
//    properties[Z] = _externalProperties(triangleIndex)[2][X];
//  }
//  else if ( dim == 2 )
//  {
//    properties[X] = _externalProperties(triangleIndex[0]);
//    properties[Y] = _externalProperties(triangleIndex[1]);
//    properties[Z] = _externalProperties(triangleIndex[2]);
//  }


//  return triangleInterpolation( triangle, properties, position );
//}

////template<class CoordType>
////CoordType SpatialDescriptorFunctionM<CoordType>::getProperties(
////    int& triangleIndex,
////    Vector<CoordType>& position)
////{
////  Triangle<CoordType> triangle;
////  triangle = _triMesh.triangle( triangleIndex );
////  Vector<CoordType> vertices(3);
////  vertices = _triMesh->getTriangleIndexes( triangleIndex );
////  Vector<CoordType> properties(3);
////  TriMesh<int> colortm( 3,0 );
////  colortm = triMesh.getVertexColors();

////  for(int i = 0; i < triMesh.getNumTriangles() ; i++)
////  {
////    colortm.addTriangle(
////          triMesh.getTriangleIndexes(i)[0],
////          triMesh.getTriangleIndexes(i)[1],
////          triMesh.getTriangleIndexes(i)[2] );
////  }

//////  intensities[X] = CoordType ( colortm.triangle(triangleIndex)[0][X] );
//////  intensities[Y] = CoordType ( colortm.triangle(triangleIndex)[1][X] );
//////  intensities[Z] = CoordType ( colortm.triangle(triangleIndex)[2][X] );
//////  CoordType intensity;
//////  intensity = triangleInterpolation( triangle, intensities, position );
//////return intensity;

////  properties[X] = colortm.triangle(triangleIndex)[0][X];
////  properties[Y] = colortm.triangle(triangleIndex)[1][X];
////  properties[Z] = colortm.triangle(triangleIndex)[2][X];

////  return triangleInterpolation( triangle, properties, position );
////}


////template class SpatialDescriptorFunctionM<int>;
//template class SpatialDescriptorFunctionM<float>;
////template class SpatialDescriptorFunctionM<double>;
////template class SpatialDescriptorFunctionM<long double>;




