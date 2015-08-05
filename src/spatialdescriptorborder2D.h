#ifndef SPATIALDESCRIPTORBORDER2D_H
#define SPATIALDESCRIPTORBORDER2D_H

#include "spatialdescriptor.h"
//#include <trimesh.h>
#include <curve.h>

template<class CoordType>
class SpatialDescriptorDistanceToBorder2D : public SpatialDescriptor<CoordType>
//class SpatialDescriptorDistanceToBorder2D : public SpatialDescriptor<float>
{
  public:

    SpatialDescriptorDistanceToBorder2D();

    //void setTriMesh(const TriMesh<CoordType>&);
//    void setBoundary(const Curve<CoordType>&);
    void setBoundary(const Curve<CoordType>&);

    void eval(const Vertices<CoordType>&,Vector<CoordType>&,Vector<CoordType>&);


private:

    //const TriMesh<CoordType>* _triMesh;
    const Curve<CoordType>* _boundary;
    //const Curve<CoordType> _boundary;

};


#endif // SPATIALDESCRIPTORBORDER2D_H
