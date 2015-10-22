/*!
 * \class  RegionAnalysis2D
 * \author Philippe Andrey (pa), INRA
 * \date   2015.04.09 - creation (pa)
 * \brief  Region analysis in two dimensions
****************************************************************/

#include "regionanalysis2d.h"
//#include "neighbourhood.h"
#include <neighbourhood.h>

#include <maths.h>
#include <programerror.h>

#include <cmath>

#define TRACE
#include <trace.h>

/*! Constructor.
****************************************************************/
template<class T>
RegionAnalysis2D<T>::RegionAnalysis2D() : RegionAnalysisBase<T>()
{
  _labelMatrix = 0;
  _valueMatrix = 0;
  _outputMatrix = 0;
}

/*! Destructor.
****************************************************************/
template<class T>
RegionAnalysis2D<T>::~RegionAnalysis2D()
{
}

/*! Returns \c true if a matrix of labels has been set.
****************************************************************/
template<class T>
bool RegionAnalysis2D<T>::hasLabelMatrix() const
{
  return _labelMatrix != 0;
}

/*! Sets the label matrix.
 *
 * Note that a reference to \c labelMatrix is taken (no hard copy).
 * Consequently, it is expected that \c labelMatrix is not destroyed
 * by the caller before this object.
****************************************************************/
template<class T>
void RegionAnalysis2D<T>::setLabelMatrix(PixelMatrix<T>& labelMatrix)
{
  _labelMatrix = &labelMatrix;
}

/*! Returns the currently set label matrix.
****************************************************************/
template<class T>
const PixelMatrix<T>& RegionAnalysis2D<T>::getLabelMatrix() const
{
  return *_labelMatrix;
}

/*! Returns \c true if a matrix of pixel values has been set.
****************************************************************/
template<class T>
bool RegionAnalysis2D<T>::hasValueMatrix() const
{
  return _valueMatrix != 0;
}

/*! Sets the matrix of pixel values.
 *
 * Note that a reference to \c valueMatrix is taken (no hard copy).
 * Consequently, it is expected that \c valueMatrix is not destroyed
 * by the caller before this object.
****************************************************************/
template<class T>
void RegionAnalysis2D<T>::setValueMatrix(const PixelMatrix<T>& valueMatrix)
{
  _valueMatrix = &valueMatrix;
}

/*! Returns the currently associated matrix of pixel values.
****************************************************************/
template<class T>
const PixelMatrix<T>& RegionAnalysis2D<T>::getValueMatrix() const
{
  return *_valueMatrix;
}

/*! Returns \c true if an output matrix has been set.
****************************************************************/
template<class T>
bool RegionAnalysis2D<T>::hasOutputMatrix() const
{
  return _outputMatrix != 0;
}

/*! Sets the output matrix.
 *
 * Note that a reference to \c outputMatrix is taken (no hard copy).
 * Consequently, it is expected that \c outputMatrix is not destroyed
 * by the caller before this object.
****************************************************************/
template<class T>
void RegionAnalysis2D<T>::setOutputMatrix(PixelMatrix<T>& outputMatrix)
{
  _outputMatrix = &outputMatrix;
}

/*! Returns the output matrix.
****************************************************************/
template<class T>
PixelMatrix<T>& RegionAnalysis2D<T>::getOutputMatrix() const
{
  return *_outputMatrix;
}

/*!
****************************************************************/
template<class T>
void RegionAnalysis2D<T>::run()
{
  const PixelMatrix<T>& labelMatrix = getLabelMatrix();
  Vector<Region>& regions = RegionAnalysisBase<T>::_regions;
  const int size1 = labelMatrix.getSize1();
  const int size2 = labelMatrix.getSize2();
  const int numRegions = static_cast<int>( labelMatrix.max().max() );
  int i, j, u, v, r, s, n, currentVertex;
  Vector<int> currentVertices;

  regions.setSize( numRegions );
  currentVertices.setZeros( numRegions );
  Vector<unsigned int> histogram = labelMatrix.histogram( 0.0, 1.0, numRegions+1 );

  for (r = 0; r < numRegions; ++r)
  {
    regions[r]._label = r+1;
    regions[r]._numFacets = 0;
    regions[r]._vertices = new Vertices<int>( 2, histogram[r+1] );
    regions[r]._edgeFlags.setSize( histogram[r+1] );
  }

  const Matrix<int> shifts = Neighbourhood::neighbourhood4();

  for (i = 0; i < size1; ++i)
    for (j = 0; j < size2; ++j)
    {
      r = static_cast<int>( labelMatrix(i,j) );
      if ( isRegionLabel(r) )
      {
        Region& region = regions[r-1];
        currentVertex = currentVertices[r-1]++;
        (*region._vertices)[currentVertex][X] = i;
        (*region._vertices)[currentVertex][Y] = j;
        region._edgeFlags[currentVertex] = false;

        for (n = 0; n < 4; ++n)
        {
          u = i + shifts(n,X);
          v = j + shifts(n,Y);
          if ( u >= 0 && u < size1
            && v >= 0 && v < size2 )
          {
            s = static_cast<int>( labelMatrix(u,v) );
            if ( s != r )
            {
              ++region._numFacets;
              region._edgeFlags[currentVertex] = true;
            }
          }
        }
      }
    }

  computeRAG();
}

/*! Computes the region adjacency graph, using the implicit boundary mode.
****************************************************************/
template<class T>
void RegionAnalysis2D<T>::computeRAGImplicitBoundaryMode()
{
  const PixelMatrix<T>& labelMatrix = getLabelMatrix();
  const Matrix<int> shifts = Neighbourhood::neighbourhood4();
  SquareMatrix<int>& adjacencyGraph = RegionAnalysisBase<T>::_adjacencyGraph;
  Vector<int> neighbour( 2 );
  int regionLabel, otherLabel;


  EVAL( shifts );

  adjacencyGraph.setZeros( numRegions() );

  for (int r = 0; r < numRegions(); ++r)
  {
    const Vertices<int>& vertices = getRegion(r).getVertices();
    regionLabel = getRegion(r).getLabel();
    for (int v = 0; v < vertices.getNumVertices(); ++v)
      for (int n = 0; n < 4; ++n)
      {
        neighbour = vertices[v] + shifts[n];
        if ( labelMatrix.contains(neighbour) )
        {
          otherLabel = static_cast<int>( labelMatrix(neighbour) );
          if ( isRegionLabel(otherLabel) && otherLabel != regionLabel )
            ++adjacencyGraph(regionLabel-1,otherLabel-1);
        }
      }
  }

  EVAL( adjacencyGraph.max().max() );
  EVAL( adjacencyGraph.sum().sum() );

  if ( adjacencyGraph.isSymmetrical() == false )
  {
    ProgramError programError;
    programError.setWhere( "void RegionAnalysis<T>::run()" );
    programError.setWhat( "Asymmetrical region adjacency graph" );
    throw programError;
  }
}

/*! Computes the region adjacency graph, using the explicit boundary mode.
****************************************************************/
template<class T>
void RegionAnalysis2D<T>::computeRAGExplicitBoundaryMode()
{
  SquareMatrix<int>& adjacencyGraph = RegionAnalysisBase<T>::_adjacencyGraph;
  const Matrix<int> neighbourhood8 = Neighbourhood::neighbourhood8();
  const PixelMatrix<T>& labelMatrix = getLabelMatrix();
  const int size1 = labelMatrix.getSize1();
  const int size2 = labelMatrix.getSize2();
  int i, j, n, label, label1, label2;
  Vector<int> position( 2 );
  Vector<int> neighbour1( 2 );
  Vector<int> neighbour2( 2 );

  adjacencyGraph.setZeros( numRegions() );

  for (i = 0; i < size1; ++i)
    for (j = 0; j < size2; ++j)
      {
        label = static_cast<int>( labelMatrix(i,j) );
        if ( label == getBoundaryValue() )
        {
          position[0] = i;
          position[1] = j;
          for (n = 0; n < 8; n += 2)
          {
            neighbour1 = position + neighbourhood8[n];
            neighbour2 = position + neighbourhood8[n+1];
            if ( labelMatrix.contains(neighbour1) && labelMatrix.contains(neighbour2) )
            {
              label1 = static_cast<int>( labelMatrix(neighbour1) );
              label2 = static_cast<int>( labelMatrix(neighbour2) );
              if ( label1 != label2 && isRegionLabel(label1) && isRegionLabel(label2) )
              {
                ++adjacencyGraph(label1-1,label2-1);
                ++adjacencyGraph(label2-1,label1-1);
                break; // quit loop over neighbourhood
              }
            }
          }
        }
      }
}

/*! Returns the values of the pixels belonging to the specified region.
****************************************************************/
template<class T>
Vector<T> RegionAnalysis2D<T>::regionValues(const Region& region) const
{
  const Vertices<int>& vertices = region.getVertices();
  const PixelMatrix<T>& pixelMatrix = getValueMatrix();
  Vector<T> values( vertices.getSize() );
  for (int v = 0; v < vertices.getSize(); ++v)
    values[v] = pixelMatrix( vertices[v] );
  return values;
}

/*! In the output matrix, fills the specified \c region with the given \c value.
 *
 * Pixels located outside of the specified \c region are not affected.
****************************************************************/
template<class T>
void RegionAnalysis2D<T>::outputFillRegion(const Region& region, const T value)
{
  if ( region.getSize() > 0 )
  {
    PixelMatrix<T>& pixelMatrix = *_outputMatrix;
    const Vertices<int>& vertices = region.getVertices();
    for (int v = 0; v < vertices.getSize(); ++v)
      pixelMatrix( vertices[v] ) = value;
  }
}

/*! Returns the compactness (circularity) of the specified region.
 *
 * Note that a rough estimation (over-estimation) of the perimeter is used.
****************************************************************/
template<class T>
float RegionAnalysis2D<T>::regionCompactness(const Region& region) const
{
  const float perimeter = region._numFacets;
  const float area = region.getSize();
  return 4.0 * M_PI * area / Maths::sqr(perimeter);
}

template<class T>
float RegionAnalysis2D<T>::regionBoundaryMeasure(const Region& region) const
{
  const Vertices<int>& vertices = region.getVertices();
  const PixelMatrix<T>& labelMatrix = getLabelMatrix();
  const int size1 = labelMatrix.getSize1();
  const int size2 = labelMatrix.getSize2();
  const float xPixelSize = spatialCalibration()[0];
  const float yPixelSize = spatialCalibration()[1];
  const int r = region.getLabel();
  float boundaryMeasure = 0.0;
  int i, j, ii, jj, n, s;

  Matrix<int> shifts( 4, 2 );
  shifts.setZeros();
  shifts(0,X) = -1;
  shifts(1,X) = 1;
  shifts(2,Y) = -1;
  shifts(3,Y) = 1;

  for (int v = 0; v < vertices.getNumVertices(); ++v)
  {
    i = vertices[v][X];
    j = vertices[v][Y];

    for (n = 0; n < 4; ++n)
    {
      ii = i + shifts(n,X);
      jj = j + shifts(n,Y);

      if ( ii >= 0 && ii < size1
        && jj >= 0 && jj < size2 )
      {
        s = static_cast<int>( labelMatrix(ii,jj) );
        if ( s != r )
        {
          if ( ii != i ) boundaryMeasure += yPixelSize;
          if ( jj != j ) boundaryMeasure += xPixelSize;
        }
      }
    }
  }

  return boundaryMeasure;
}

/*! Returns the spatial calibration of the associated label matrix.
****************************************************************/
template<class T>
Vector<float> RegionAnalysis2D<T>::spatialCalibration() const
{
  return getLabelMatrix().getPixelCalibration().getPixelSize();
}

template class RegionAnalysis2D<float>;
template class RegionAnalysis2D<double>;
template class RegionAnalysis2D<long double>;
