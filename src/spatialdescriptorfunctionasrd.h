#ifndef SPATIALDESCRIPTORFUNCTIONASRD_H
#define SPATIALDESCRIPTORFUNCTIONASRD_H

#include <spatialdescriptor.h>

template<class CoordType>
class SpatialDescriptorFunctionASRD : public SpatialDescriptor<CoordType>
{
  public:

    SpatialDescriptorFunctionASRD();

    void eval(const Vertices<CoordType>&,Vector<CoordType>&,Vector<CoordType>&);

};

#endif // SPATIALDESCRIPTORFUNCTIONASRD_H

