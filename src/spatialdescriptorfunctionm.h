#ifndef SPATIALDESCRIPTORFUNCTIONM_H
#define SPATIALDESCRIPTORFUNCTIONM_H


#include "spatialdescriptor.h"
#include <trimesh.h>


template<class CoordType>
class SpatialDescriptorFunctionM : public SpatialDescriptor<CoordType>
{
  public:

    SpatialDescriptorFunctionM();

    void setTriMesh(const TriMesh<CoordType>&);
    void setPropertiesTriMesh( const Shape<CoordType>& );
    //void setPropertiesTriMesh(const TriMesh<CoordType>&);

    //void setPointPatterns(const Vertices<CoordType>&);

    void eval(const Vertices<CoordType>&,
              Vector<CoordType>&,
              Vector<CoordType>&);

//    CoordType getProperties(//TriMesh<CoordType>&,
//                             int&,
//                             Vector<CoordType>& );

    CoordType triangleInterpolation(const Triangle<CoordType>&,
                                    const Vector<CoordType>&,
                                    const Vector<CoordType>&);
    //TriMesh<CoordType> getPropertiesMap();


private:

    bool _externalProperties;
    const TriMesh<CoordType>* _triMesh;
    const Shape<CoordType>* _propertiesMap;

    CoordType getValue( const Vector<CoordType>&,
                             const int& = 0 );

//    const Vertices<CoordType>* _pointPatterns;

};



#endif // SPATIALDESCRIPTORFUNCTIONM_H
