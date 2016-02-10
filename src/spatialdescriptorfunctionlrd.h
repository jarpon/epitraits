#ifndef SPATIALDESCRIPTORFUNCTIONLRD_H
#define SPATIALDESCRIPTORFUNCTIONLRD_H

#include <spatialdescriptor.h>

template<class CoordType>
class SpatialDescriptorFunctionLRD : public SpatialDescriptor<CoordType>
{
  public:

    SpatialDescriptorFunctionLRD();

    void setDistanceThreshold(const CoordType&);
    const CoordType& getDistanceThreshold() const;

    void eval(const Vertices<CoordType>&,Vector<CoordType>&,Vector<CoordType>&);

  private:

    CoordType _distanceThreshold;
};

#endif // SPATIALDESCRIPTORFUNCTIONLRD_H
