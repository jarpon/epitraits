#ifndef SPATIALDESCRIPTORCENTROID_H
#define SPATIALDESCRIPTORCENTROID_H

#include "spatialdescriptor.h"
#include <trimesh.h>

template<class CoordType>
class SpatialDescriptorDistanceToCentroid : public SpatialDescriptor<CoordType>
{
  public:

    SpatialDescriptorDistanceToCentroid();

    void setTriMesh(const TriMesh<CoordType>&);

    void eval(const Vertices<CoordType>&,Vector<CoordType>&,Vector<CoordType>&);


private:

    const TriMesh<CoordType>* _triMesh;

};


#endif // SPATIALDESCRIPTORCENTROID_H
