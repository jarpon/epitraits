#ifndef SPATIALMODELBORDERDISTANCE3DDIFFERENTCOMPARTMENTS_H
#define SPATIALMODELBORDERDISTANCE3DDIFFERENTCOMPARTMENTS_H

#include "trimeshspatialmodeldifferentcompartments.h"

template<class CoordType>
class SpatialModelBorderDistance3DDifferentCompartments : public virtual TriMeshSpatialModelDifferentCompartments<CoordType>
{
  public:

    SpatialModelBorderDistance3DDifferentCompartments();

    void setDistancesToBorder(const Vector<CoordType>&, const Vector<CoordType>&);
    const Vector<CoordType>& getDistancesToBorderDistribution1() const;
    const Vector<CoordType>& getDistancesToBorderDistribution2() const;

    //Vertices<CoordType> drawSample(const int);
    Vertices<CoordType> drawSample(const int, const int);


  protected:

    void drawPositionFromBorder(Vector<CoordType>&, const CoordType);

  private:

    Vector<CoordType> _distancesToBorderDistribution1, _distancesToBorderDistribution2, _distancesToBorder, _classBelongings;
};

#endif // SPATIALMODELBORDERDISTANCE3DDIFFERENTCOMPARTMENTS_H
