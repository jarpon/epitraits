#ifndef SPATIALMODELMAXIMALREPULSION3D_H
#define SPATIALMODELMAXIMALREPULSION3D_H

#include <spatialmodelhardcoredistance3d.h>

template<class CoordType>
class SpatialModelMaximalRepulsion3D : public SpatialModelHardcoreDistance3D<CoordType>
{
  public:

    SpatialModelMaximalRepulsion3D();

    void setNumMonteCarloCycles(const int);
    int getNumMonteCarloCycles() const;

    CoordType getBeta() const;
    void initializeBeta(const int);
    Vertices<CoordType> drawSample(const int);

    const Vector<CoordType>& getEnergyProfile() const;

  protected:

    virtual CoordType energy(const Vertices<CoordType>&) const;

  private:

    void moveVertex(Vertices<CoordType>&, const int, const CoordType);
    bool checkHardcoreDistances(const Vector<CoordType>&, const int, const Vertices<CoordType>&) const;

    CoordType _beta;
    Vector<CoordType> _energyProfile;
    int _numMonteCarloCycles;
};

#endif // SPATIALMODELMAXIMALREPULSION3D_H
