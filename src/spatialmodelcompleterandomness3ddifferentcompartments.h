#ifndef SPATIALMODELCOMPLETERANDOMNESS3DDIFFERENTCOMPARTMENTS_H
#define SPATIALMODELCOMPLETERANDOMNESS3DDIFFERENTCOMPARTMENTS_H

#include "trimeshspatialmodeldifferentcompartments.h"

template<class CoordType>
class SpatialModelCompleteRandomness3DDifferentCompartments : public TriMeshSpatialModelDifferentCompartments<CoordType>
{
  public:

    SpatialModelCompleteRandomness3DDifferentCompartments();

    Vertices<CoordType> drawSample(const int);
};

#endif // SPATIALMODELCOMPLETERANDOMNESS3DDIFFERENTCOMPARTMENTS_H
