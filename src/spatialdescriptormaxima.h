#ifndef SPATIALDESCRIPTORMAXIMA_H
#define SPATIALDESCRIPTORMAXIMA_H

#include "spatialdescriptor.h"
#include <trimesh.h>

template<class CoordType>
class SpatialDescriptorMaximaRepulsion : public SpatialDescriptor<CoordType>
{
  public:

    SpatialDescriptorMaximaRepulsion();

    void setTriMesh(const TriMesh<CoordType>&);

    void eval(const Vertices<CoordType>&,Vector<CoordType>&,Vector<CoordType>&);


private:

    const TriMesh<CoordType>* _triMesh;

};

#endif // SPATIALDESCRIPTORMAXIMA_H
