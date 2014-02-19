#ifndef TRIMESHSPATIALMODEL_H
#define TRIMESHSPATIALMODEL_H

#include <alinevector.h>
#include <trimesh.h>
#include <randomgenerator.h>
#include <vertices.h>
#include <spatialmodel.h>

//using namespace std;

template<class CoordType>
class TrimeshSpatialModel : public SpatialModel<CoordType,float>
{
  public:

    TrimeshSpatialModel();
    ~TrimeshSpatialModel();

//    void setRandomGenerator(RandomGenerator&);
    void setTriMesh(const TriMesh<CoordType>&);
    void setOutput(string output) { _outputDir = output ; }

    void setNumberOfCompartments(const int);

//    void setVolumeRadius(const T volumeRadius) { _volumeRadius = volumeRadius; }
//    void setVolumeRadius(Vector<T>&);
//    void setVolumeRadiusRange(const T volumeRadiusMin, const T volumeRadiusMax) { _volumeRadiusRange[0] = volumeRadiusMin; _volumeRadiusRange[1] = volumeRadiusMax; }

    void setHardcoreDistance(Vector<CoordType>&);
    void setHardcoreDistance(const CoordType);
    void setHardcoreDistanceRange(const CoordType hardcoreDistanceMin, const CoordType hardcoreDistanceMax) { _hardcoreDistanceRange[0] = hardcoreDistanceMin; _hardcoreDistanceRange[1] = hardcoreDistanceMax; }

    void setDistanceToBorder(const CoordType);
    void setDistanceToBorder(Vector<CoordType>&);
    void setDistanceToBorderRange( const CoordType distanceToBorderMin, const CoordType distanceToBorderMax ) { _distanceToBorderRange[0] = distanceToBorderMin; _distanceToBorderRange[1] = distanceToBorderMax; }

    void drawPosition(Vector<CoordType>&);
    void drawPositionFromBorder(Vector<CoordType>&, CoordType distanceToBorder);

    bool checkHardcoreDistances(const Vector<CoordType>&, const Vertices<CoordType>&);
    bool checkDistancesToBorder(const Vector<CoordType>&, CoordType distanceToBorder);

    Vertices<CoordType> drawSample(const int);
//    virtual Vertices<CoordType> drawSample(const int) = 0;
//    ShapeSet<CoordType> drawSamples(const int,const int);

    const CoordType& getHardcoreDistance(const int) const;
    Vector<CoordType>& getHardcoreDistance();
    //const CoordType& getDistanceToBorder();

    void initialize();

    void save(const string);

  private:

    Vertices<CoordType> randomVertices();
    Vertices<CoordType> hardcoreDistance();
    Vertices<CoordType> distanceToTheBorder();
    Vertices<CoordType> hardcoreAndToTheBorderDistances();

    const TriMesh<CoordType>* _triMesh;
    Vertices<CoordType> _vertices;

    string _outputDir;
    int _numCompartments;

    //CoordType _volumeRadius;
    //Vector<CoordType> _volumeRadiusRange;

    Vector<CoordType> _hardcoreDistance;
    Vector<CoordType> _hardcoreDistanceRange;

    Vector<CoordType> _distanceToBorder;
    Vector<CoordType> _distanceToBorderRange;

    float xMinBB, xMaxBB, yMinBB, yMaxBB, zMinBB, zMaxBB;

};

#endif // RANDOMPOINTSGENERATOR_H
