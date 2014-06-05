#ifndef TRIMESHSPATIALMODEL_H
#define TRIMESHSPATIALMODEL_H

#include <alinevector.h>
#include <trimesh.h>
#include <randomgenerator.h>
#include <vertices.h>
#include <spatialmodel.h>
#include <trimeshquery.h>
//#include "maximarepulsion.h"

//using namespace std;

template<class CoordType>
class TriMeshSpatialModel : public SpatialModel<CoordType,float>
{
  public:

    TriMeshSpatialModel();
    ~TriMeshSpatialModel();

    void setTriMesh(const TriMesh<CoordType>&);
    void setOutput(string output) { _outputDir = output ; }

    void setNumberOfCompartments(const int);

    void setHardcoreDistances(const Vector<CoordType>&);
    void setHardcoreDistances(const CoordType);
    void setHardcoreDistancesRange(const CoordType hardcoreDistancesMin, const CoordType hardcoreDistancesMax) { _hardcoreDistancesRange[0] = hardcoreDistancesMin; _hardcoreDistancesRange[1] = hardcoreDistancesMax; }

    void setDistanceToBorder(const CoordType);
    void setDistanceToBorder(const Vector<CoordType>&);
    void setDistanceToBorderRange( const CoordType distanceToBorderMin, const CoordType distanceToBorderMax ) { _distanceToBorderRange[0] = distanceToBorderMin; _distanceToBorderRange[1] = distanceToBorderMax; }

    void drawPosition(Vector<CoordType>&);
    void drawPositionFromBorder(Vector<CoordType>&, CoordType distanceToBorder);

//    bool checkHardcoreDistances(const Vector<CoordType>&, const Vertices<CoordType>&);
//    bool checkHardcoreDistances(const Vector<CoordType>&, const int, const Vertices<CoordType>&) const;
//    bool checkDistancesToBorder(const Vector<CoordType>&, CoordType distanceToBorder);

    bool checkObjectToBorderDistance(const Vector<CoordType>&, const int&);
    bool checkInterObjectDistances(const Vector<CoordType>&, const Vertices<CoordType>&);


    Vertices<CoordType> drawSample(const int);

    const CoordType& getHardcoreDistances(const int) const;
    Vector<CoordType>& getHardcoreDistances();
    //const CoordType& getDistanceToBorder();

    void initialize();

    Vertices<CoordType> hardcoreDistances();
    int _numCompartments;

    const TriMesh<CoordType>* _triMesh;
    TriMeshQuery<CoordType> _triMeshQuery;
    Vector<CoordType> _hardcoreDistances;
    Vertices<CoordType> _vertices;

  private:

    Vertices<CoordType> randomVertices();
    Vertices<CoordType> distanceToTheBorder();
    Vertices<CoordType> hardcoreAndToTheBorderDistances();

    void randomizesOrder(RandomGenerator&);



    string _outputDir;

    Vector<CoordType> _hardcoreDistancesRange;

    Vector<CoordType> _distanceToBorder;
    Vector<CoordType> _distanceToBorderRange;

    float xMinBB, xMaxBB, yMinBB, yMaxBB, zMinBB, zMaxBB;

};

#endif // TRIMESHSPATIALMODEL_H
