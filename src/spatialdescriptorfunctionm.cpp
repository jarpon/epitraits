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
  _externalProperties = 0;
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
}

/*! Sets an external 3D mesh map with the properties of the boundary.
****************************************************************/
template<class CoordType>
void SpatialDescriptorFunctionM<CoordType>::setExternalProperties(const TriMesh<CoordType>& triMesh)
{
  _propertiesTriMesh = &triMesh;
  _externalProperties = 3;
}

/*! Sets an external curve map with the properties of the boundary.
****************************************************************/
template<class CoordType>
void SpatialDescriptorFunctionM<CoordType>::setExternalProperties(const Curve<CoordType>& curve)
{
  _propertiesCurve= &curve;
  _externalProperties = 4;
}

/*! Sets an external vertices map with the properties of the boundary.
****************************************************************/
template<class CoordType>
void SpatialDescriptorFunctionM<CoordType>::setExternalProperties(const Vertices<CoordType>& vertices)
{
  _propertiesVertices = &vertices;
  _externalProperties = 2;
}

/*! Sets an external vector with the properties of the boundary.
****************************************************************/
template<class CoordType>
void SpatialDescriptorFunctionM<CoordType>::setExternalProperties(const Vector<CoordType>& vector)
{
  _propertiesVector = &vector;
  _externalProperties = 1;
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
    //case 2: eval2D( vertices, xEcdf, yEcdf );
    case 3: eval3D( vertices, xEcdf, yEcdf );
  }
}

///*! Computes the ECDF of the property of the closest point within the boundary
// *  on the pattern given by \c vertices for 2D cases
// *
// * Returns the coordinates \c xEcdf, \c yEcdf of the points of the ECDF.
//****************************************************************/
//template<class CoordType>
//void SpatialDescriptorFunctionM<CoordType>::eval2D(
//  const Vertices<CoordType>& vertices,
//  Vector<CoordType>& xEcdf,
//  Vector<CoordType>& yEcdf)
//{
//  const int numVertices = vertices.getSize();
//  Vector<CoordType> closestPoint;

//  xEcdf.setSize( numVertices );
//  for (int i = 0; i < numVertices; ++i)
//  {
//    _curve->closestPoint( vertices[i], closestPoint );
//    xEcdf[i] = get2DValue( vertices[i], closestPoint );
//  }
//  xEcdf.sort();

//  CDFTools<CoordType> cdftools;
//  yEcdf = cdftools.cdf( xEcdf );
//}

/*! Computes the ECDF of the property of the closest point within the boundary
 *  on the pattern given by \c vertices for 3D cases
 *
 * Returns the coordinates \c xEcdf, \c yEcdf of the points of the ECDF.
****************************************************************/
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
    _triMeshQuery.closestPoint( vertices[i], triMeshVertex, triangleIndex );
    xEcdf[i] = get3DValue( vertices[i], triangleIndex );
  }
  xEcdf.sort();

  CDFTools<CoordType> cdftools;
  yEcdf = cdftools.cdf( xEcdf );
}

/*!
 * Get the intensity at position \par position in the triangle \par triangle.
 *
 * The intensity is interpolated from the ones at the triangle vertices.
 */
template<class CoordType>
Vector<CoordType> SpatialDescriptorFunctionM<CoordType>::useTriMesh(
//    const Vector<CoordType>& position,
    const int& triangleIndex)
{
  Vector<CoordType> properties(3);
//  Triangle<CoordType> triangle;
  Triangle<int> colorTriangle;
  TriMesh<int> tempTriMesh;

  if ( _externalProperties == 3 )
    tempTriMesh = _propertiesTriMesh->getVertexColors();
  else
    tempTriMesh = _triMeshQuery.getTriMesh().getVertexColors();

  colorTriangle = tempTriMesh.triangle( triangleIndex );
//  triangle = _triMeshQuery.getTriMesh().triangle( triangleIndex );

  properties[0] = colorTriangle[0][X];
  properties[1] = colorTriangle[1][X];
  properties[2] = colorTriangle[2][X];

  return properties;
//  return ( ( ( properties[0] * vertex.distance( triangle[0] ) ) + ( properties[1] * vertex.distance( triangle[1] ) ) + ( properties[2] * vertex.distance( triangle[2] ) ) )
//      / ( vertex.distance( triangle[0] ) + vertex.distance(triangle[1]) + vertex.distance(triangle[2]) ) );
}

///*!
// * Get the intensity at position \par position in the triangle \par triangle.
// *
// * The intensity is interpolated from the ones at the triangle vertices.
// */
//template<class CoordType>
//Vector<CoordType> SpatialDescriptorFunctionM<CoordType>::useCurve(
////    const Vector<CoordType>& position,
//    const int& triangleIndex)
//{
//  Vector<CoordType> properties(3);
////  Triangle<CoordType> triangle;
//  Triangle<int> colorTriangle;
//  TriMesh<int> tempTriMesh;

//  if ( _externalProperties == 3 )
//    tempTriMesh = _propertiesTriMesh->getVertexColors();
//  else
//    tempTriMesh = _triMeshQuery.getTriMesh().getVertexColors();

//  colorTriangle = tempTriMesh.triangle( triangleIndex );
////  triangle = _triMeshQuery.getTriMesh().triangle( triangleIndex );

//  properties[0] = colorTriangle[0][X];
//  properties[1] = colorTriangle[1][X];
//  properties[2] = colorTriangle[2][X];

//  return properties;
////  return ( ( ( properties[0] * vertex.distance( triangle[0] ) ) + ( properties[1] * vertex.distance( triangle[1] ) ) + ( properties[2] * vertex.distance( triangle[2] ) ) )
////      / ( vertex.distance( triangle[0] ) + vertex.distance(triangle[1]) + vertex.distance(triangle[2]) ) );
//}

/*!
 * Get the intensity at position \par position in the triangle \par triangle.
 *
 * The intensity is interpolated from the ones at the triangle vertices.
 */
template<class CoordType>
Vector<CoordType> SpatialDescriptorFunctionM<CoordType>::useVertices(
//    const Vector<CoordType>& vertex,
    const int& triangleIndex)
{
  Vector<CoordType> properties(3);
//  Triangle<CoordType> triangle;

  Vector<unsigned int> positions(3);
  positions = _triMeshQuery.getTriMesh().getTriangleIndexes( triangleIndex );

//  triangle = _triMeshQuery.getTriMesh( ).triangle( triangleIndex );

  properties[X] = _propertiesVertices[positions[0]][0][X];
  properties[Y] = _propertiesVertices[positions[1]][0][X];
  properties[Z] = _propertiesVertices[positions[2]][0][X];

  return properties;
//  return ( ( ( properties[0] * vertex.distance( triangle[0] ) ) + ( properties[1] * vertex.distance( triangle[1] ) ) + ( properties[2] * vertex.distance( triangle[2] ) ) )
//      / ( vertex.distance( triangle[0] ) + vertex.distance(triangle[1]) + vertex.distance(triangle[2]) ) );
}


/*!
 * Get the intensity at position \par position in the triangle \par triangle.
 *
 * The intensity is interpolated from the ones at the triangle vertices.
 */
template<class CoordType>
Vector<CoordType> SpatialDescriptorFunctionM<CoordType>::useVector(
//    const Vector<CoordType>& vertex,
    const int& triangleIndex)
{
  Vector<CoordType> properties(3);

  Vector<unsigned int> positions(3);
  positions = _triMeshQuery.getTriMesh().getTriangleIndexes( triangleIndex );



  properties[X] = _propertiesVector[positions[0]][X];
  properties[Y] = _propertiesVector[positions[1]][X];
  properties[Z] = _propertiesVector[positions[2]][X];

  return properties;
//  return ( ( ( properties[0] * vertex.distance( triangle[0] ) ) + ( properties[1] * vertex.distance( triangle[1] ) ) + ( properties[2] * vertex.distance( triangle[2] ) ) )
//      / ( vertex.distance( triangle[0] ) + vertex.distance(triangle[1]) + vertex.distance(triangle[2]) ) );}
}

/*!
 *
 */
template<class CoordType>
CoordType SpatialDescriptorFunctionM<CoordType>::get3DValue(
    const Vector<CoordType>& vertex,
    const int& triangleIndex )
{
  Vector<CoordType> propertyValues;
  Triangle<CoordType> triangle;

  triangle = _triMeshQuery.getTriMesh( ).triangle( triangleIndex );

  switch( _externalProperties )
  {
    case 0: propertyValues = useTriMesh( triangleIndex );
    case 1: propertyValues = useVector( triangleIndex );
    case 2: propertyValues = useVertices( triangleIndex );
    case 3: propertyValues = useTriMesh( triangleIndex );
  }

//  return value;
  return ( ( ( propertyValues[0] * vertex.distance( triangle[0] ) ) + ( propertyValues[1] * vertex.distance( triangle[1] ) ) + ( propertyValues[2] * vertex.distance( triangle[2] ) ) )
      / ( vertex.distance( triangle[0] ) + vertex.distance(triangle[1]) + vertex.distance(triangle[2]) ) );

}

///*!
// *
// */
//template<class CoordType>
//CoordType SpatialDescriptorFunctionM<CoordType>::get2DValue(
//    const Vector<CoordType>& vertex,
//    const int& triangleIndex )
//{
//  Vector<CoordType> propertyValues;
//  Triangle<CoordType> triangle;

//  switch( _externalProperties )
//  {
//    case 0: propertyValues = useCurve( triangleIndex );
//    case 1: propertyValues = useVector( triangleIndex );
//    case 4: propertyValues = useCurve( triangleIndex );
//  }

////  return value;
//  return ( ( ( propertyValues[0] * vertex.distance( triangle[0] ) ) + ( propertyValues[1] * vertex.distance( triangle[1] ) ) + ( propertyValues[2] * vertex.distance( triangle[2] ) ) )
//      / ( vertex.distance( triangle[0] ) + vertex.distance(triangle[1]) + vertex.distance(triangle[2]) ) );

//}


template class SpatialDescriptorFunctionM<float>;
template class SpatialDescriptorFunctionM<double>;
template class SpatialDescriptorFunctionM<long double>;




