#ifndef SPATIALDESCRIPTORFUNCTIONLRD_H
#define SPATIALDESCRIPTORFUNCTIONLRD_H

#include <spatialdescriptor.h>

template<class CoordType>
class SpatialDescriptorFunctionLRD : public SpatialDescriptor<CoordType>
{
  public:

    SpatialDescriptorFunctionLRD();

    void eval(const Vertices<CoordType>&,Vector<CoordType>&,Vector<CoordType>&);
};

#endif // SPATIALDESCRIPTORFUNCTIONLRD_H
