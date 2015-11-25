/*!
 * \class  SpatialModelCompleteRandomness3D
 * \author Javier Arpón (ja), INRA
 * \author Philippe Andrey (pa), INRA
 * \date    - creation (ja)
 * \date    - integration (pa)
 * \brief  Complete spatial randomness model in 3D to analyze different compartments
****************************************************************/

#include "spatialmodelcompleterandomness3ddifferentcompartments.h"

/*! Constructor.
****************************************************************/
template<class CoordType>
SpatialModelCompleteRandomness3DDifferentCompartments<CoordType>::SpatialModelCompleteRandomness3DDifferentCompartments() : TriMeshSpatialModelDifferentCompartments<CoordType>()
{
}


/*! Generates \c numVertices vertices distributed uniformly at random
 * within the specified 3D domain.
****************************************************************/
template<class CoordType>
Vertices<CoordType> SpatialModelCompleteRandomness3DDifferentCompartments<CoordType>::drawSample(const int numVertices)
{
  Vertices<CoordType> vertices( 3, numVertices );

  for (int v = 0; v < numVertices; ++v)
    this->drawPosition( vertices[v] );

  return vertices;
}

template class SpatialModelCompleteRandomness3DDifferentCompartments<float>;
template class SpatialModelCompleteRandomness3DDifferentCompartments<double>;
template class SpatialModelCompleteRandomness3DDifferentCompartments<long double>;
