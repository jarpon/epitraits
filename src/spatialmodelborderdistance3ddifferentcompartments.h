#ifndef SPATIALMODELBORDERDISTANCE3DDIFFERENTCOMPARTMENTS_H
#define SPATIALMODELBORDERDISTANCE3DDIFFERENTCOMPARTMENTS_H

#include "trimeshspatialmodeldifferentcompartments.h"

template<class CoordType>
class SpatialModelBorderDistance3DDifferentCompartments : public virtual TriMeshSpatialModelDifferentCompartments<CoordType>
{
  public:

    SpatialModelBorderDistance3DDifferentCompartments();

    void setDistancesToBorder(const Vector<CoordType>&);
    const Vector<CoordType>& getDistancesToBorder() const;

    Vertices<CoordType> drawSample(const int);

  protected:

    void drawPositionFromBorder(Vector<CoordType>&, const CoordType);

  private:

    Vector<CoordType> _distancesToBorder;
};

#endif // SPATIALMODELBORDERDISTANCE3DDIFFERENTCOMPARTMENTS_H
