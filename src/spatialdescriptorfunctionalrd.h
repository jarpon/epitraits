#ifndef SPATIALDESCRIPTORFUNCTIONALRD_H
#define SPATIALDESCRIPTORFUNCTIONALRD_H

#include <spatialdescriptor.h>

template<class CoordType>
class SpatialDescriptorFunctionALRD : public SpatialDescriptor<CoordType>
{
  public:

    SpatialDescriptorFunctionALRD();

    void eval(const Vertices<CoordType>&,Vector<CoordType>&,Vector<CoordType>&);

};

#endif // SPATIALDESCRIPTORFUNCTIONALRD_H
