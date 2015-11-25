#ifndef SPATIALMODELHARDCOREDISTANCE3DDIFFERENTCOMPARTMENTS_H
#define SPATIALMODELHARDCOREDISTANCE3DDIFFERENTCOMPARTMENTS_H

#include "trimeshspatialmodeldifferentcompartments.h"

template<class CoordType>
class SpatialModelHardcoreDistance3DDifferentCompartments : public virtual TriMeshSpatialModelDifferentCompartments<CoordType>
{
  public:

    SpatialModelHardcoreDistance3DDifferentCompartments();

    void setHardcoreDistances(const Vector<CoordType>&);
    const Vector<CoordType>& getHardcoreDistances() const;

    Vertices<CoordType> drawSample(const int);

  protected:

    void shuffleDistances(Vector<CoordType>&);
    bool validObjectToBorderDistance(const Vector<CoordType>&, const CoordType&);
    bool validInterObjectDistances(const Vector<CoordType>&, const Vertices<CoordType>&,const Vector<CoordType>&);

  private:

    Vector<CoordType> _hardcoreDistances;
};

#endif // SPATIALMODELHARDCOREDISTANCE3DDIFFERENTCOMPARTMENTS_H
