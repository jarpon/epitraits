#ifndef SPATIALDESCRIPTORFUNCTIONSRD_H
#define SPATIALDESCRIPTORFUNCTIONSRD_H

#include <spatialdescriptor.h>

template<class CoordType>
class SpatialDescriptorFunctionSRD : public SpatialDescriptor<CoordType>
{
  public:

    SpatialDescriptorFunctionSRD();

    void setDistanceThreshold(const CoordType&);
    const CoordType& getDistanceThreshold() const;

    void eval(const Vertices<CoordType>&,Vector<CoordType>&,Vector<CoordType>&);

  private:

    CoordType _distanceThreshold;
};

#endif // SPATIALDESCRIPTORFUNCTIONSRD_H

