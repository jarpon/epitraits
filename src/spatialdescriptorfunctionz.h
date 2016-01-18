#ifndef SPATIALDESCRIPTORFUNCTIONZ_H
#define SPATIALDESCRIPTORFUNCTIONZ_H

#include <spatialdescriptor.h>

template<class CoordType>
class SpatialDescriptorFunctionZ : public SpatialDescriptor<CoordType>
{
  public:

    SpatialDescriptorFunctionZ();

    void eval(const Vertices<CoordType>&,Vector<CoordType>&,Vector<CoordType>&);
};

#endif // SPATIALDESCRIPTORFUNCTIONZ_H
