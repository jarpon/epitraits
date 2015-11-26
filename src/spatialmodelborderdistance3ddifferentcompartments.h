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

    const Vertices<CoordType>& getVerticesDistribution1() const;
    const Vertices<CoordType>& getVerticesDistribution2() const;

  protected:

    Vertices<CoordType> sortVertices(const Vertices<CoordType>&, const Vector<CoordType>&);
    void shuffleObjects(Vector<CoordType>&);
    void drawPositionFromBorder(Vector<CoordType>&, const CoordType);

  private:

    Vector<CoordType> _distancesToBorderDistribution1, _distancesToBorderDistribution2, _distancesToBorder, _classBelongings;
    int _numVerticesDist1, _numVerticesDist2;

    Matrix<int> _classCoordinates;
    Vertices<CoordType> _verticesDist1;
    Vertices<CoordType> _verticesDist2;
};

#endif // SPATIALMODELBORDERDISTANCE3DDIFFERENTCOMPARTMENTS_H
