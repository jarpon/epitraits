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

    void setExternalProperties( const TriMesh<CoordType>& );
    void setExternalProperties( const Curve<CoordType>& );
    void setExternalProperties( const Vertices<CoordType>& );
    void setExternalProperties( const Vector<CoordType>& );

    void eval(const Vertices<CoordType>&, Vector<CoordType>&, Vector<CoordType>&);

private:

    const Curve<CoordType>* _curve;
    TriMeshQuery<CoordType> _triMeshQuery;

    int _externalProperties;
    const TriMesh<CoordType>*  _propertiesTriMesh;
    const Curve<CoordType>*    _propertiesCurve;
    const Vertices<CoordType>* _propertiesVertices;
    const Vector<CoordType>*   _propertiesVector;

    CoordType get3DValue( const Vector<CoordType>&, const int& = 0 );
//    CoordType get2DValue( const Vector<CoordType>&, const int& = 0 );

    Vector<CoordType> useTriMesh(const int&);
    Vector<CoordType> useCurve(const int&);
    Vector<CoordType> useVertices(const int&);
    Vector<CoordType> useVector(const int&);


    void eval3D(const Vertices<CoordType>&, Vector<CoordType>&, Vector<CoordType>&);
//    void eval2D(const Vertices<CoordType>&, Vector<CoordType>&, Vector<CoordType>&);

};

#endif // SPATIALDESCRIPTORFUNCTIONM_H
