#ifndef MAXIMALREPULSION_H
#define MAXIMALREPULSION_H

#include <randomgenerator.h>
#include <vertices.h>
#include <trimeshquery.h>
#include <dataset.h>
#include "trimeshspatialmodel.h"

template<class CoordType>
class MaximalRepulsionTriMeshSpatialModel : public TriMeshSpatialModel<float>
{
  public:

    MaximalRepulsionTriMeshSpatialModel();
    ~MaximalRepulsionTriMeshSpatialModel();

    void setTriMesh(const TriMesh<CoordType>&);
//    void setTriMeshSpatialModel(TriMeshSpatialModel<CoordType> &);

    //void setMaximalRepulsion(const Vertices<CoordType>&);

    Vertices<CoordType> drawSample(const int);

    void setHardcoreDistance(const Vector<CoordType>&);

  private:

    CoordType getEnergy(const Vertices<CoordType>&, const int&);
    CoordType getEnergy1(const Vertices<CoordType>&);
    CoordType getEnergy2(const Vertices<CoordType>&);
    CoordType getEnergy3(const Vertices<CoordType>&);
    CoordType getEnergy4(const Vertices<CoordType>&);
    CoordType findBeta();

    Vector<CoordType> moveCompartment(const int, const Vertices<CoordType>&, const CoordType);

    bool checkHardcoreDistances(const Vector<CoordType>&, const int, const Vertices<CoordType>&) const;

//    const TriMesh<CoordType>* _triMesh;
//    TriMeshQuery<CoordType> _triMeshQuery;
//    Vertices<float> _vertices;
//    const TriMeshSpatialModel<CoordType>* _triMeshSpatialModel;

//    Vector<CoordType> _hardcoreDistances;
    CoordType _nucleusVolume;

};

#endif // MAXIMALREPULSION_H
