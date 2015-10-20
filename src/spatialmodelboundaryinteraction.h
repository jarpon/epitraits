#ifndef SPATIALMODELBOUNDARYINTERACTION_H
#define SPATIALMODELBOUNDARYINTERACTION_H

#include <spatialmodel2d.h>

template<class CoordType>
class SpatialModelBoundaryInteraction : public SpatialModel2D<CoordType,CoordType>
{
  public:

    SpatialModelBoundaryInteraction();
    ~SpatialModelBoundaryInteraction();

    void setMargin(const CoordType);
    CoordType getMargin() const;
    void setMarginProb(const float);
    float getMarginProb() const;

    void initialize();
    const Curve<CoordType>& getInnerBoundary() const;
    Vertices<CoordType> drawSample(const int);

  private:

    CoordType _margin;
    float _marginProb;
    Curve<CoordType> _innerBoundary;
    BoundingBox<CoordType> _innerBoundingBox;
};

#endif // SPATIALMODELBOUNDARYINTERACTION_H
