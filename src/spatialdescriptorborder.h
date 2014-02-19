#ifndef SPATIALDESCRIPTORBORDER_H
#define SPATIALDESCRIPTORBORDER_H

#include "spatialdescriptor.h"
#include <trimesh.h>

template<class CoordType>
class SpatialDescriptorDistanceToBorder : public SpatialDescriptor<CoordType>
{
  public:

    SpatialDescriptorDistanceToBorder();

    void setTriMesh(const TriMesh<CoordType>&);

    void eval(const Vertices<CoordType>&,Vector<CoordType>&,Vector<CoordType>&);


private:

    const TriMesh<CoordType>* _triMesh;

};

#endif // SPATIALDESCRIPTORBORDER_H
