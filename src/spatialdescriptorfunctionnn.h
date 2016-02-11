#ifndef SPATIALDESCRIPTORFUNCTIONNN_H
#define SPATIALDESCRIPTORFUNCTIONNN_H

#include "spatialdescriptor.h"

template<class CoordType>
class SpatialDescriptorFunctionNN : public SpatialDescriptor<CoordType>
{
  public:

    SpatialDescriptorFunctionNN();

    void eval(const Vertices<CoordType>&,Vector<CoordType>&,Vector<CoordType>&);
};

#endif // SPATIALDESCRIPTORFUNCTIONNN_H
