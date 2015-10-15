/*!
 * \class  SpatialDescriptorFunctionM
 * \author Javier Arpon (ja), INRA
 * \author Thotino Gobin-Gansou (tg-g), INRA
 * \date   2015.10.15 - creation (ja)
 * \brief  M-function statistics for spatial point patterns
 * \details This class implements the property-of-the-closest-point-of-the-boundary
 * cumulative distribution. This is one of the available spatial descriptors
 * for point patterns.
 * The class handles the 2D and the 3D cases.
 * \brief  Study the relationship of volumes
 * respect the property of the closest point of the boundary
****************************************************************/

#include "spatialdescriptorfunctionm.h"

#include "spatialdescriptor.h"

#include <cdftools.h>
#include <trimesh.h>
#include <cmath>
#include <dataset.h>
#include <trimeshquery.h>
#include <sstream>

#define TRACE
#include <trace.h>


/*! Constructor.
****************************************************************/
template<class CoordType>
SpatialDescriptorFunctionM<CoordType>::SpatialDescriptorFunctionM() : SpatialDescriptor<CoordType>()
{
  _curve = 0;
}

/*! Sets the boundary curve (used only in the 2D case).
****************************************************************/
template<class CoordType>
void SpatialDescriptorFunctionM<CoordType>::setCurve(const Curve<CoordType>& curve)
{
  _curve = &curve;
}

/*! Sets the boundary mesh (used only in the 3D case).
****************************************************************/
template<class CoordType>
void SpatialDescriptorFunctionM<CoordType>::setTriMesh(const TriMesh<CoordType>& triMesh)
{
  _triMeshQuery.setTriMesh( triMesh );
  _externalProperties = false;
}

/*! Computes the ECDF of the property of the closest point within the boundary
 *  on the pattern given by \c vertices.
 *
 * Returns the coordinates \c xEcdf, \c yEcdf of the points of the ECDF.
****************************************************************/
template<class CoordType>
void SpatialDescriptorFunctionM<CoordType>::eval(
  const Vertices<CoordType>& vertices,
  Vector<CoordType>& xEcdf,
  Vector<CoordType>& yEcdf)
{
  switch( vertices.getSpace() )
  {
    case 2: eval2D( vertices, xEcdf, yEcdf );
    case 3: eval3D( vertices, xEcdf, yEcdf );
  }
}

template<class CoordType>
void SpatialDescriptorFunctionM<CoordType>::eval2D(
  const Vertices<CoordType>& vertices,
  Vector<CoordType>& xEcdf,
  Vector<CoordType>& yEcdf)
{
  const int numVertices = vertices.getSize();
  Vector<CoordType> closestPoint;

  xEcdf.setSize( numVertices );
  for (int i = 0; i < numVertices; ++i)
    xEcdf[i] = _curve->closestPoint( vertices[i], closestPoint );
  xEcdf.sort();

  CDFTools<CoordType> cdftools;
  yEcdf = cdftools.cdf( xEcdf );
}

template<class CoordType>
void SpatialDescriptorFunctionM<CoordType>::eval3D(
  const Vertices<CoordType>& vertices,
  Vector<CoordType>& xEcdf,
  Vector<CoordType>& yEcdf)
{
  const int numVertices = vertices.getSize();
  int triangleIndex;
  Vector<CoordType> vertex(3), triMeshVertex(3);

  xEcdf.setSize( numVertices );
  for (int i = 0; i < numVertices; ++i)
  {
    vertex = vertices[i];
    _triMeshQuery.closestPoint( vertex, triMeshVertex, triangleIndex );
    xEcdf[i] = getValue( vertex, triangleIndex );
  }
  xEcdf.sort();

  CDFTools<CoordType> cdftools;
  yEcdf = cdftools.cdf( xEcdf );
}


template<class CoordType>
void SpatialDescriptorFunctionM<CoordType>::setPropertiesTriMesh(const Shape<CoordType>& properties )
{
  _propertiesMap = &properties;
  _externalProperties = true;

  // in this example the property is the intensity map of the envelope
//  TriMesh<CoordType> triMeshPropertyMap( _propertiesTriMesh->getSpace(), _propertiesTriMesh->getSize() );

//  for ( int i = 0; i < _propertiesTriMesh->getSize(); ++i )
//  {
////    triMeshPropertyMap[i][X] = int(_propertiesTriMesh[i][X]);
////    triMeshPropertyMap[i][Y] = int(_propertiesTriMesh[i][Y]);
////    triMeshPropertyMap[i][Z] = int(_propertiesTriMesh[i][Z]);
//      triMeshPropertyMap[i][X] = _propertiesTriMesh[i][X];
//      triMeshPropertyMap[i][Y] = _propertiesTriMesh[i][Y];
//      triMeshPropertyMap[i][Z] = _propertiesTriMesh[i][Z];

//  }

//  for ( int i = 0; i < _propertiesTriMesh->getNumTriangles(); ++i )
//  {
//    //EVAL(_propertiesTriMesh.getTriangleIndexes(i));
//    triMeshPropertyMap.addTriangle(
//          _propertiesTriMesh->getTriangleIndexes(i)[0],
//          _propertiesTriMesh->getTriangleIndexes(i)[1],
//          _propertiesTriMesh->getTriangleIndexes(i)[2] );
//  }

  //_triMesh.setVertexColors( new TriMesh<CoordType>( triMeshPropertyMap ) );
//  _triMesh->setVertexColors( new TriMesh<CoordType>( _propertiesTriMesh ) );
//  _triMesh->setVertexColors( _propertiesTriMesh );
}

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

/*!
 * Get the intensity at position \par position in the triangle \par triangle.
 *
 * The intensity is interpolated from the ones at the triangle vertices.
 */
template<class CoordType>
CoordType SpatialDescriptorFunctionM<CoordType>::triangleInterpolation(
    const Triangle<CoordType>& triangle,
    const Vector<CoordType>& position,
    bool externalInfo)
{
  Vector<CoordType> properties(3);
  Triangle<CoordType> triangle;

  if ( externalInfo == false )
  {
    triangle = _triMeshQuery->triangle( triangleIndex );

//    properties[X] = _triMesh->triangle(triangleIndex)[0][X];
//    properties[Y] = _triMesh->triangle(triangleIndex)[1][X];
//    properties[Z] = _triMesh->triangle(triangleIndex)[2][X];
    properties[X] = triangle[0][X];
    properties[Y] = triangle[1][X];
    properties[Z] = triangle[2][X];
  }
  else
  {
    triangle = _externalProperties->triangle( triangleIndex );
    properties[X] = triangle[0][X];
    properties[Y] = triangle[1][X];
    properties[Z] = triangle[2][X];
  }

  return ( ( ( properties[0] * position.distance( triangle[0] ) ) + ( properties[1] * position.distance( triangle[1] ) ) + ( properties[2] * position.distance( triangle[2] ) ) )
      / ( position.distance( triangle[0] ) + position.distance(triangle[1]) + position.distance(triangle[2]) ) );
}

/*!
 * Get the intensity at position \par position in the triangle \par triangle.
 *
 * The intensity is interpolated from the ones at the triangle vertices.
 */
template<class CoordType>
CoordType SpatialDescriptorFunctionM<CoordType>::verticesInterpolation(
    const Triangle<CoordType>& vertices,
    const Vector<CoordType>& position)
{
  Vector<CoordType> properties(3);

  properties[X] = _externalProperties[vertices[0]];
  properties[Y] = _externalProperties[vertices[1]];
  properties[Z] = _externalProperties[vertices[2]];

  return ( ( ( properties[0] * position.distance( triangle[0] ) ) + ( properties[1] * position.distance( triangle[1] ) ) + ( properties[2] * position.distance( triangle[2] ) ) )
      / ( position.distance( triangle[0] ) + position.distance(triangle[1]) + position.distance(triangle[2]) ) );
}

/*!
 *
 */
template<class CoordType>
CoordType SpatialDescriptorFunctionM<CoordType>::getValue(
    const Vector<CoordType>& position,
    const int& triangleIndex )
{
  Triangle<CoordType> triangle;
  triangle = _triMesh->triangle( triangleIndex );
  Vector<unsigned int> vertices(3);
  vertices = _triMesh->getTriangleIndexes( triangleIndex );
  Vector<CoordType> properties(3);
  TriMesh<int> colortm( 3,0 );
  colortm = _triMesh->getVertexColors();

  CoordType value;

  switch( _propertiesMap->getShapeType() )
  {
    case Vertices: value = verticesInterpolation( vertices, position );
    case TriMesh: value = triangleInterpolation( triangle, position, _externalProperties );
    default value = triangleInterpolation( triangle, position, _externalProperties );
  }

  return value;
  //return triangleInterpolation( triangle, properties, position );
}

//template<class CoordType>
//CoordType SpatialDescriptorFunctionM<CoordType>::getProperties(
//    int& triangleIndex,
//    Vector<CoordType>& position)
//{
//  Triangle<CoordType> triangle;
//  triangle = _triMesh.triangle( triangleIndex );
//  Vector<CoordType> vertices(3);
//  vertices = _triMesh->getTriangleIndexes( triangleIndex );
//  Vector<CoordType> properties(3);
//  TriMesh<int> colortm( 3,0 );
//  colortm = triMesh.getVertexColors();

//  for(int i = 0; i < triMesh.getNumTriangles() ; i++)
//  {
//    colortm.addTriangle(
//          triMesh.getTriangleIndexes(i)[0],
//          triMesh.getTriangleIndexes(i)[1],
//          triMesh.getTriangleIndexes(i)[2] );
//  }

////  intensities[X] = CoordType ( colortm.triangle(triangleIndex)[0][X] );
////  intensities[Y] = CoordType ( colortm.triangle(triangleIndex)[1][X] );
////  intensities[Z] = CoordType ( colortm.triangle(triangleIndex)[2][X] );
////  CoordType intensity;
////  intensity = triangleInterpolation( triangle, intensities, position );
////return intensity;

//  properties[X] = colortm.triangle(triangleIndex)[0][X];
//  properties[Y] = colortm.triangle(triangleIndex)[1][X];
//  properties[Z] = colortm.triangle(triangleIndex)[2][X];

//  return triangleInterpolation( triangle, properties, position );
//}


template class SpatialDescriptorFunctionM<float>;
template class SpatialDescriptorFunctionM<double>;
template class SpatialDescriptorFunctionM<long double>;




