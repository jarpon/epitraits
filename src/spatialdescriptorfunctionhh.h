#ifndef SPATIALDESCRIPTORFUNCTIONHH
#define SPATIALDESCRIPTORFUNCTIONHH

#include "spatialdescriptor.h"

template<class CoordType>
class SpatialDescriptorFunctionHH : public SpatialDescriptor<CoordType>
{
  public:

    SpatialDescriptorFunctionHH();

    void setVertices ( const Vertices<CoordType>& );
    void setOrder ( const Vector<CoordType>& );
//    void setOrder ( const Vector<string>& );

    void setVertices ( const Vertices<CoordType>&, const Vertices<CoordType>& );
    void setVerticesKind1 ( const Vertices<CoordType>& );
    void setVerticesKind2 ( const Vertices<CoordType>& );

    void eval(const Vertices<CoordType>&,Vector<CoordType>&,Vector<CoordType>&);

  private:

    Vertices<CoordType> _allVertices;
    Vector<CoordType> _allVerticesOrder;
    Vector<CoordType> _allVerticesLabelOrder;

    Vertices<CoordType> _verticesKind1;
    Vertices<CoordType> _verticesKind2;

};

#endif // SPATIALDESCRIPTORFUNCTIONHH

