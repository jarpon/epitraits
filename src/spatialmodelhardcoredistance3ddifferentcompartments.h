#ifndef SPATIALMODELHARDCOREDISTANCE3DDIFFERENTCOMPARTMENTS_H
#define SPATIALMODELHARDCOREDISTANCE3DDIFFERENTCOMPARTMENTS_H

#include "trimeshspatialmodeldifferentcompartments.h"
//#include <trimeshspatialmodel.h>
#include <alinematrix.h>
#include <vertices.h>

template<class CoordType>
class SpatialModelHardcoreDistance3DDifferentCompartments : public virtual TriMeshSpatialModelDifferentCompartments<CoordType>
{
  public:

    SpatialModelHardcoreDistance3DDifferentCompartments();

    void setHardcoreDistances(const Vector<CoordType>&, const Vector<CoordType>&);
    const Vector<CoordType>& getHardcoreDistancesDistribution1() const;
    const Vector<CoordType>& getHardcoreDistancesDistribution2() const;

    Vertices<CoordType> drawSample(const int, const int);

    const Vertices<CoordType>& getVerticesDistribution1() const;
    const Vertices<CoordType>& getVerticesDistribution2() const;

  protected:

    Vertices<CoordType> sortVertices(const Vertices<CoordType>&, const Vector<CoordType>&);
    void shuffleObjectsOrder(Vector<CoordType>&);
    bool validObjectToBorderDistance(const Vector<CoordType>&, const CoordType&);
    bool validInterObjectDistances(const Vector<CoordType>&, const Vertices<CoordType>&,const Vector<CoordType>&);

  private:

    Vector<CoordType> _hardcoreDistancesDistribution1, _hardcoreDistancesDistribution2, _hardcoreDistances, _classBelongings;
    int _numVerticesDist1, _numVerticesDist2;

    Matrix<int> _classCoordinates;
    Vertices<CoordType> _verticesDist1;
    Vertices<CoordType> _verticesDist2;
};

#endif // SPATIALMODELHARDCOREDISTANCE3DDIFFERENTCOMPARTMENTS_H
