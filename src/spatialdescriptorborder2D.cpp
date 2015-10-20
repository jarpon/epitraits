///*!
// * \class  SpatialDescriptorDistanceToBorder
// * \author Javier Arpón (ja), INRA
// * \date   2015.06.15 - creation (ja)
// * \brief  Study the peripheral positioning of points into a 2D contour.
//****************************************************************/

//#include "spatialdescriptorborder2D.h"

//#include <cdftools.h>
////#include <trimesh.h>
//#include <curve.h>

////#define TRACE
//#include <trace.h>

//template<class CoordType>
//SpatialDescriptorDistanceToBorder2D<CoordType>::SpatialDescriptorDistanceToBorder2D() : SpatialDescriptor<CoordType>()
//{
//  //_boundary = 0;
//}

//template<class CoordType>
//void SpatialDescriptorDistanceToBorder2D<CoordType>::setBoundary(const Curve<CoordType>& boundary)
//{
//  _boundary = &boundary;
//  //_boundary = boundary;
//}

//template<class CoordType>
//void SpatialDescriptorDistanceToBorder2D<CoordType>::eval(
//  const Vertices<CoordType>& vertices,
//  Vector<CoordType>& xvalues,
//  Vector<CoordType>& yvalues)
//{
//  ENTER( "void SpatialDescriptorDistanceToBorder2D<CoordType>::eval(...)" );
//    PRINT("here");
//  const int numVertices = vertices.getSize();

//  Vector<CoordType> testVertex(2);
//  Vector<CoordType> boundaryVertex(2);
//  xvalues.setSize( numVertices );

//  PRINT("here");
//  for (int i = 0; i < numVertices; ++i)
//  {
//    testVertex = vertices[i];
//    _boundary->closestPoint( testVertex, boundaryVertex ); //!!!
//    PRINT("here!!");

//    xvalues[i] = testVertex.distance( boundaryVertex );
//  }
//  xvalues.sort();
//  PRINT("here!!!!!!!");

//  CDFTools<CoordType> cdftools;
//  yvalues = cdftools.cdf( xvalues );
//  PRINT("here!!!!!!!!!!!!");

//  LEAVE();
//}

//template class SpatialDescriptorDistanceToBorder2D<float>;
////template class SpatialDescriptorDistanceToBorder2D<double>;
////template class SpatialDescriptorDistanceToBorder2D<long double>;
