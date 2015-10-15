#ifndef SPATIALDESCRIPTORFUNCTIONM_H
#define SPATIALDESCRIPTORFUNCTIONM_H

#include "spatialdescriptor.h"

#include <trimeshquery.h>

template<class CoordType>
class SpatialDescriptorFunctionM : public SpatialDescriptor<CoordType>
{
  public:

    SpatialDescriptorFunctionM();

    void setCurve(const Curve<CoordType>&);
    void setTriMesh(const TriMesh<CoordType>&);
    void setPropertiesTriMesh( const Shape<CoordType>& );

    void eval(const Vertices<CoordType>&, Vector<CoordType>&, Vector<CoordType>&);

private:

    bool _externalProperties;
    const TriMesh<CoordType>* _triMesh;
    const Shape<CoordType>* _propertiesMap;
    const Curve<CoordType>* _curve;
    TriMeshQuery<CoordType> _triMeshQuery;

    CoordType getValue( const Vector<CoordType>&, const int& = 0 );
    CoordType triangleInterpolation(const Triangle<CoordType>&, const Vector<CoordType>&, bool);
    CoordType verticesInterpolation(const Vector<CoordType>&, const Vector<CoordType>&);

    void eval3D(const Vertices<CoordType>&, Vector<CoordType>&, Vector<CoordType>&);
    void eval2D(const Vertices<CoordType>&, Vector<CoordType>&, Vector<CoordType>&);

};

#endif // SPATIALDESCRIPTORFUNCTIONM_H
