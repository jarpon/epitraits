#ifndef TRIMESHSPATIALMODELDIFFERENTCOMPARTMENTS_H
#define TRIMESHSPATIALMODELDIFFERENTCOMPARTMENTS_H

#include <alinevector.h>
#include <trimesh.h>
#include <randomgenerator.h>
#include <vertices.h>
#include <spatialmodel.h>
#include <trimeshquery.h>

//using namespace std;

template<class CoordType>
class TriMeshSpatialModelDifferentCompartments : public SpatialModel<CoordType,float>
{
  public:

    TriMeshSpatialModelDifferentCompartments();
    ~TriMeshSpatialModelDifferentCompartments();

    void setTriMesh(const TriMesh<CoordType>&);
    const TriMesh<CoordType>& getTriMesh() const;
    const TriMeshQuery<CoordType>& getTriMeshQuery() const;

    void initialize();

  protected:

    void drawPosition(Vector<CoordType>&);

  private:

    TriMeshQuery<CoordType> _triMeshQuery;
    BoundingBox<CoordType> _boundingBox;

};

#endif // TRIMESHSPATIALMODELDIFFERENTCOMPARTMENTS_H
