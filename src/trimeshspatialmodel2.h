#ifndef TriMeshSpatialModel22_H
#define TriMeshSpatialModel22_H

#include <alinevector.h>
#include <trimesh.h>
#include <randomgenerator.h>
#include <vertices.h>
#include <spatialmodel.h>
#include <trimeshquery.h>
//#include "maximarepulsion.h"

//using namespace std;

template<class CoordType>
class TriMeshSpatialModel2 : public SpatialModel<CoordType,float>
{
  public:

    TriMeshSpatialModel2();
    ~TriMeshSpatialModel2();

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
    Vertices<CoordType>  drawSample2( ); // vertices of two distributions

    const CoordType& getHardcoreDistances(const int) const;
    Vector<CoordType>& getHardcoreDistances();
    //const CoordType& getDistanceToBorder();

    void initialize();

    Vertices<CoordType> hardcoreDistances();




    // new version: working with more than one objects distribution
    Vector<string> allDistributionObjectsNames; // set in order to know whether there is an alteration in the order of the objects

    void addDistribution(); // adds current info (radii and/or distances to the border) to the model
    void addDistribution( const string& ); // adds current info (radii and/or distances to the border) to the model including the name of the objects kind
    void addDistribution( const int, const Vector<CoordType>& ); // adds the number of objects and radii of a new objects distribution
    void addDistribution( const int, const Vector<CoordType>&, const Vector<CoordType>& ); // adds the number of objects, radii and distances to the border of a new objects distribution
    void addDistribution( const int, const string&, const Vector<CoordType>&, const Vector<CoordType>& ); // gets the number of objects, kind of objects, radii and distances to the border of a new objects distribution
    void getDistribution( const int, int, Vector<CoordType>, Vector<CoordType> ); // gets the number of objects, radii and distances to the border of the X objects distribution
    //void getFirstDistribution( int, Vector<CoordType>, Vector<CoordType> ); // gets the number of objects, radii and distances to the border of the first objects distribution
    //void getSecondDistribution( int, Vector<CoordType>, Vector<CoordType> ); // gets the number of objects, radii and distances to the border of the second objects distribution

    void addObjectsName( const int, const string& );
    void addObjectsNames( const Vector<string>& );

    int getObjectsNumber (const string& );
    string getObjectsName( const int );
    Vector<string> getObjectsNames();

    int _numDistributions; // number of different kind of compartments/objects
    int _numCompartments; // number of compartments/objects current kind

    Vector<string> _objectsNames;
    int getNumDistributions();

    //void drawSample( Vertices<CoordType>&, Vertices<CoordType>& ); // vertices of two distributions
    void drawSample( Vertices<CoordType>&, Vector<string>& ); // vertices of any distribution and vector with the corresponding one

  private:

    Vertices<CoordType> randomVertices();
    Vertices<CoordType> distanceToTheBorder();
    Vertices<CoordType> hardcoreAndToTheBorderDistances();

//    void randomizesOrder(RandomGenerator&);

    const TriMesh<CoordType>* _triMesh;
    TriMeshQuery<CoordType> _triMeshQuery;
    Vector<CoordType> _hardcoreDistances;
    Vertices<CoordType> _vertices;
    void randomizesOrder(RandomGenerator&);

    string _outputDir;

    Vector<CoordType> _hardcoreDistancesRange;

    Vector<CoordType> _distanceToBorder;
    Vector<CoordType> _distanceToBorderRange;

    Vector<int> _allDistributionNumObjects;
    Vector<CoordType> _allDistributionHardcoreDistances;
    Vector<CoordType> _allDistributionDistancesToTheBorder;

    float xMinBB, xMaxBB, yMinBB, yMaxBB, zMinBB, zMaxBB;

};


#endif // TriMeshSpatialModel22_H

