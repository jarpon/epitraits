/*!
 * \class  RegionAnalysis
 * \author Philippe Andrey (pa), INRA
 * \date    - création (pa)
 * \brief  Analyse de régions
****************************************************************/

#include "regionanalysis3d.h"
#include "componentlabelling.h"
#include "neighbourhood.h"
#include "sumfilter.h"

#include <maths.h>
#include <matherror.h>
#include <boundingbox.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <programerror.h>

#define TRACE
#include <trace.h>

/*! Constructor.
****************************************************************/
template<class T>
RegionAnalysis3D<T>::RegionAnalysis3D() : RegionAnalysisBase<T>()
{
  _labelMatrix = 0;
  _valueMatrix = 0;
  _outputMatrix = 0;
}

/*! Destructor.
****************************************************************/
template<class T>
RegionAnalysis3D<T>::~RegionAnalysis3D()
{
}

/*! Returns \c true if a matrix of labels has been set.
****************************************************************/
template<class T>
bool RegionAnalysis3D<T>::hasLabelMatrix() const
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
void RegionAnalysis3D<T>::setLabelMatrix(VoxelMatrix<T>& labelMatrix)
{
  _labelMatrix = &labelMatrix;
}

/*! Returns the currently set label matrix.
****************************************************************/
template<class T>
const VoxelMatrix<T>& RegionAnalysis3D<T>::getLabelMatrix() const
{
  return *_labelMatrix;
}

/*! Returns the currently set label matrix.
****************************************************************/
template<class T>
VoxelMatrix<T>& RegionAnalysis3D<T>::getLabelMatrix()
{
  return *_labelMatrix;
}

/*! Returns \c true if a matrix of voxel values has been set.
****************************************************************/
template<class T>
bool RegionAnalysis3D<T>::hasValueMatrix() const
{
  return _valueMatrix != 0;
}

/*! Sets the matrix of voxel values.
 *
 * Note that a reference to \c valueMatrix is taken (no hard copy).
 * Consequently, it is expected that \c valueMatrix is not destroyed
 * by the caller before this object.
****************************************************************/
template<class T>
void RegionAnalysis3D<T>::setValueMatrix(const VoxelMatrix<T>& valueMatrix)
{
  _valueMatrix = &valueMatrix;
}

/*! Returns the currently associated matrix of voxel values.
****************************************************************/
template<class T>
const VoxelMatrix<T>& RegionAnalysis3D<T>::getValueMatrix() const
{
  return *_valueMatrix;
}

/*! Returns \c true if an output matrix has been set.
****************************************************************/
template<class T>
bool RegionAnalysis3D<T>::hasOutputMatrix() const
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
void RegionAnalysis3D<T>::setOutputMatrix(VoxelMatrix<T>& outputMatrix)
{
  _outputMatrix = &outputMatrix;
}

/*! Returns the output matrix.
****************************************************************/
template<class T>
VoxelMatrix<T>& RegionAnalysis3D<T>::getOutputMatrix() const
{
  return *_outputMatrix;
}

/*! Applique une table de conversion des valeurs des labels.
 *
 * Chaque voxel est remplacé par la valeur qui lui correspond
 * dans la table.
 *
 * La taille de \c labelMap doit être égale au label maximum + 1
 * (+1 pour le fond qui correspond à la valeur 0).
 \todo à placer dans une classe LabelFilter
****************************************************************/
template<class T>
void RegionAnalysis3D<T>::applyLabelMap(const Vector<int>& labelMap)
{
  VoxelMatrix<T>& regionMatrix = *_labelMatrix;
  const int size1 = regionMatrix.getSize1();
  const int size2 = regionMatrix.getSize2();
  const int size3 = regionMatrix.getSize3();
  int i, j, k, label;

  for (i = 0; i < size1; i++)
  {
    for (j = 0; j < size2; j++)
    {
      for (k = 0; k < size3; k++)
      {
        label = int( regionMatrix(i,j,k) );
        regionMatrix(i,j,k) = labelMap[label];
      }
    }
  }
}

/*! Applique un vote majoritaire local sur l'image des labels de régions.
 *
 * On attribue à chaque voxel la valeur la plus représentée dans
 * son voisinage.
 * Le cas échéant, le voxel conserve sa valeur initiale si elle
 * fait partie des ex-aequos.
****************************************************************/
template<class T>
void RegionAnalysis3D<T>::applyMajorityFilter()
{
  ENTER( "void RegionAnalysis<T>::applyMajorityFilter()" );

  VoxelMatrix<T>& regionMatrix = *_labelMatrix;
  const VoxelMatrix<T> bufferMatrix = regionMatrix;
  const int size1 = regionMatrix.getSize1();
  const int size2 = regionMatrix.getSize2();
  const int size3 = regionMatrix.getSize3();
  Vectori values( 27 );
  Vectori counts( 27 );
  int numValues, p, pmax;
  int i, j, k, u, v, w, label;

  for (i = 0; i < size1; i++)
  {
    for (j = 0; j < size2; j++)
    {
      for (k = 0; k < size3; k++)
      {
        label = int( bufferMatrix(i,j,k) );
        if ( label > 0 )
        {
          values[0] = label;
          counts[0] = 0;
          numValues = 1;

          for (u = i-1; u <= i+1; u++)
          for (v = j-1; v <= j+1; v++)
          for (w = k-1; w <= k+1; w++)
          {
            if ( u >= 0 && u < size1
              && v >= 0 && v < size2
              && w >= 0 && w < size3 )
            {
              label = int( bufferMatrix(u,v,w) );
              if ( label > 0 )
              {
                values[numValues] = label; // sentinelle
                p = 0;
                while ( values[p] != label )
                {
                  p++;
                }
                if ( p == numValues )
                {
                  numValues++;
                  counts[p] = 1;
                }
                else
                {
                  counts[p]++;
                }
              }
            }
          }
          pmax = 0; // correspond à la valeur initiale du voxel
          for (p = 1; p < numValues; p++)
          {
            if ( counts[p] > counts[pmax] )
            {
              pmax = p;
            }
          }
          regionMatrix(i,j,k) = values[pmax];
          if ( regionMatrix(i,j,k) == 0 )
          {
            ProgramError programError;
            programError.setWhere( "void RegionAnalysis<T>::applyMajorityFilter()" );
            programError.setWhat( "Unexpected internal error: got null region label" );
            throw programError;
          }
        }
      }
    }
  }

  LEAVE();
}

/*! Recompute region labels to ensure a continuous labelling of
 * the regions, starting from 1.
 * For example, if there were initially three regions with labels
 * 2, 7 and 11, they will be relabelled 1, 2, and 3.
 *
 * The analysis (method #run) should be re-run after a call to this
 * function.
 *
 \todo à mettre dans une classe LabelFilter
****************************************************************/
template<class T>
int RegionAnalysis3D<T>::condenseRegionLabels()
{
  VoxelMatrix<T>& labelMatrix = getLabelMatrix();
  const int maxLabel = static_cast<int>( labelMatrix.max().max().max() );
  const Vector<unsigned int> histogram = labelMatrix.histogram( 0, 1, maxLabel+1 );
  Vector<int> labelMap( maxLabel+1 );
  int i, numLabels = 0;

  labelMap.setZeros();
  for (i = 1; i <= maxLabel; ++i)
    if ( histogram[i] > 0 )
      labelMap[i] = ++numLabels;
  applyLabelMap( labelMap );
  return numLabels;
}

template<class T>
void RegionAnalysis3D<T>::run()
{
  ENTER( "void RegionAnalysis<T>::run()" );

  const VoxelMatrix<T>& labelMatrix = getLabelMatrix();
  Vector<Region>& regions = RegionAnalysisBase<T>::_regions;
  const int size1 = labelMatrix.getSize1();
  const int size2 = labelMatrix.getSize2();
  const int size3 = labelMatrix.getSize3();
  const int numRegions = int( labelMatrix.max().max().max() );
  int i, j, k, u, v, w, r, s, n, currentVertex;
  Vector<int> currentVertices;

  regions.setSize( numRegions );
  currentVertices.setZeros( numRegions );
  Vector<unsigned int> histogram = labelMatrix.histogram( 0.0, 1.0, numRegions+1 );

  EVAL( histogram.sum() );

  PRINT("alloc");

  for (r = 0; r < numRegions; ++r)
  {
    regions[r]._label = r+1;
    regions[r]._numFacets = 0;
    regions[r]._vertices = new Vertices<int>( 3, histogram[r+1] );
    regions[r]._edgeFlags.setSize( histogram[r+1] );
  }

  PRINT("done");

  Matrix<int> shifts( 6, 3 );
  shifts.setZeros();
  shifts(0,X) = -1;
  shifts(1,X) = 1;
  shifts(2,Y) = -1;
  shifts(3,Y) = 1;
  shifts(4,Z) = -1;
  shifts(5,Z) = 1;

  for (i = 0; i < size1; ++i)
    for (j = 0; j < size2; ++j)
      for (k = 0; k < size3; ++k)
      {
        r = static_cast<int>( labelMatrix(i,j,k) );
        if ( isRegionLabel(r) && r > 0 )
        {
          Region& region = regions[r-1];
          currentVertex = currentVertices[r-1]++;
          (*region._vertices)[currentVertex][X] = i;
          (*region._vertices)[currentVertex][Y] = j;
          (*region._vertices)[currentVertex][Z] = k;
          region._edgeFlags[currentVertex] = false;

          // check utilité de cette portion
          for (n = 0; n < 6; ++n)
          {
            u = i + shifts(n,X);
            v = j + shifts(n,Y);
            w = k + shifts(n,Z);

            if ( u >= 0 && u < size1
              && v >= 0 && v < size2
              && w >= 0 && w < size3 )
            {
              s = static_cast<int>( labelMatrix(u,v,w) );
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

  LEAVE();
}

/*! Computes the region adjacency graph, using the implicit boundary mode.
****************************************************************/
template<class T>
void RegionAnalysis3D<T>::computeRAGImplicitBoundaryMode()
{
  const VoxelMatrix<T>& labelMatrix = getLabelMatrix();
  const Matrix<int> shifts = Neighbourhood::neighbourhood6();
  SquareMatrix<int>& adjacencyGraph = RegionAnalysisBase<T>::_adjacencyGraph;
  Vector<int> neighbour( 3 );
  int regionLabel, otherLabel;

  adjacencyGraph.setZeros( numRegions() );

  for (int r = 0; r < numRegions(); ++r)
  {
    const Vertices<int>& vertices = getRegion(r).getVertices();
    regionLabel = getRegion(r).getLabel();
    for (int v = 0; v < vertices.getNumVertices(); ++v)
      for (int n = 0; n < 6; ++n)
      {
        neighbour = vertices[v] + shifts[n];
        if ( labelMatrix.contains(neighbour) )
        {
          otherLabel = static_cast<int>( labelMatrix(neighbour) );
          if ( otherLabel > 0 && otherLabel != regionLabel )
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
void RegionAnalysis3D<T>::computeRAGExplicitBoundaryMode()
{
  SquareMatrix<int>& adjacencyGraph = RegionAnalysisBase<T>::_adjacencyGraph;
  const Matrix<int> neighbourhood26 = Neighbourhood::neighbourhood26();
  const VoxelMatrix<T>& labelMatrix = getLabelMatrix();
  const int size1 = labelMatrix.getSize1();
  const int size2 = labelMatrix.getSize2();
  const int size3 = labelMatrix.getSize3();
  int i, j, k, n, label, label1, label2;
  Vector<int> position( 3 );
  Vector<int> neighbour1( 3 );
  Vector<int> neighbour2( 3 );

  adjacencyGraph.setZeros( numRegions() );

  for (i = 0; i < size1; ++i)
    for (j = 0; j < size2; ++j)
      for (k = 0; k < size3; ++k)
      {
        label = static_cast<int>( labelMatrix(i,j,k) );
        if ( label == getBoundaryValue() )
        {
          position[0] = i;
          position[1] = j;
          position[2] = k;
          for (n = 0; n < 26; n += 2)
          {
            neighbour1 = position + neighbourhood26[n];
            neighbour2 = position + neighbourhood26[n+1];
            if ( labelMatrix.contains(neighbour1) && labelMatrix.contains(neighbour2) )
            {
              label1 = static_cast<int>( labelMatrix(neighbour1) );
              label2 = static_cast<int>( labelMatrix(neighbour2) );
              if ( isRegionLabel(label1) && isRegionLabel(label2) && label1 > 0 && label2 > 0 )
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
PixelMatrix<T> RegionAnalysis3D<T>::getLabel2DProjection( const T label) const
{
  VoxelMatrix<T>& regionMatrix = *_labelMatrix;
  const int size1 = regionMatrix.getSize1();
  const int size2 = regionMatrix.getSize2();
  const int size3 = regionMatrix.getSize3();

  PixelMatrix<T> labelProjection( size1, size2 );
  labelProjection.setZeros();

  for ( int k = 0; k < size3; ++k )
    for (int i = 0; i < size1; ++i )
      for ( int j = 0; j < size2; ++j )
        if ( regionMatrix(i,j,k) == label )
          labelProjection(i,j) = 1;

  PixelCalibration pixelCalibration;
  pixelCalibration.setPixelHeight( regionMatrix.getVoxelCalibration().getVoxelHeight() );
  pixelCalibration.setPixelWidth( regionMatrix.getVoxelCalibration().getVoxelWidth() );
  pixelCalibration.setLengthUnit( regionMatrix.getVoxelCalibration().getLengthUnit() );

  labelProjection.setPixelCalibration( pixelCalibration );

  return labelProjection;
}

/*! In the specified \c voxelMatrix, fills the specified \c region with the given \c value.
 *
 * Voxels located outside the specified \c region are not affected.
****************************************************************/
template<class T>
void RegionAnalysis3D<T>::fillRegion(
  const Region& region,
  const T value,
  VoxelMatrix<T>& voxelMatrix)
{
  if ( region.getSize() > 0 )
  {
    const Vertices<int>& vertices = *region._vertices;
    for (int v = 0; v < vertices.getSize(); ++v)
      voxelMatrix(vertices[v]) = value;
  }
}

/*! In the output matrix, fills the specified \c region with the given \c value.
 *
 * Voxels located outside the specified \c region are not affected.
****************************************************************/
template<class T>
void RegionAnalysis3D<T>::outputFillRegion(const Region& region, const T value)
{
  if ( region.getSize() > 0 )
  {
#if 0
    VoxelMatrix<T>& voxelMatrix = *_outputMatrix;
    const Vertices<int>& vertices = *region._vertices;
    int i, j, k;

    for (int v = 0; v < vertices.getSize(); ++v)
    {
      i = vertices[v][X];
      j = vertices[v][Y];
      k = vertices[v][Z];
      voxelMatrix(i,j,k) = value;
    }
#else
    fillRegion( region, value, *_outputMatrix );
#endif
  }
}

/*! Retourne les valeurs des voxels de \c voxelMatrix qui
 * appartiennent à la région \c region.
****************************************************************/
template<class T>
Vector<T> RegionAnalysis3D<T>::regionValues(const Region& region) const
{
  const VoxelMatrix<T>& voxelMatrix = getValueMatrix();
  const Vertices<int>& vertices = *region._vertices;
  Vector<T> values( vertices.getSize() );

  for (int v = 0; v < vertices.getSize(); ++v)
    values[v] = voxelMatrix( vertices[v] );

  return values;
}

template<class T>
Vector<T> RegionAnalysis3D<T>::randomOutsideValues(
  const Region& region,
  const VoxelMatrix<T>& voxelMatrix,
  RandomGenerator& randomGenerator) const
{
  const int numRegions = _regions.getSize();
  const int numValues = region.getSize();
  const int totalNumValues = this->computeRegionFeature( REGION_FEATURE_SIZE ).sum() - numValues;
  Vector<T> values( numValues );
  int v, rv, r, p, cumul;

  for (v = 0; v < numValues; v++)
  {
    rv = randomGenerator.uniformL( totalNumValues );
    cumul = 0;
    for (r = 0; r < numRegions; r++)
    {
      if ( &_regions[r] == &region ) continue;
      if ( _regions[r].getSize() + cumul > rv )
      {
        p = rv - cumul;
        const Vectori& vertex = (*_regions[r]._vertices)[p];
        values[v] = voxelMatrix( vertex[X], vertex[Y], vertex[Z] );
        break;
      }
      cumul += _regions[r].getSize();
    }
  }

  return values;
}


/*! Retourne le barycentre, coordonnées exprimées en voxel */
template<class T>
Vector<T> RegionAnalysis3D<T>::regionBarycenter(
  const Region& region,
  const VoxelMatrix<T>& voxelMatrix) const
{
  const Vertices<int>& vertices = *region._vertices;
  T svx = 0.0;
  T svy = 0.0;
  T svz = 0.0;
  T sv = 0.0;
  T value;
  int i, j, k;

  for (int v = 0; v < vertices.getSize(); ++v)
  {
    i = vertices[v][X];
    j = vertices[v][Y];
    k = vertices[v][Z];
    value = voxelMatrix(i,j,k);
    sv += value;
    svx += value * i;
    svy += value * j;
    svz += value * k;
  }

  Vector<T> barycenter( 3 );
  barycenter[X] = svx / sv;
  barycenter[Y] = svy / sv;
  barycenter[Z] = svz / sv;

  return barycenter;
}

template<class T>
T RegionAnalysis3D<T>::regionKurtosis(const Region& region,const VoxelMatrix<T>& voxelMatrix) const
{
  const Vector<T> barycenter = regionBarycenter( region, voxelMatrix );
  const int dimension = barycenter.getSize();
  SquareMatrix<T> covarianceMatrix( dimension );
  const Verticesi& vertices = *region._vertices;
  const int numVertices = vertices.getSize();
  T value, sumValues = 0.0, kurtosis = 0.0;
  Vector<T> buffer( 3 );
  int i, j, k, v;

  covarianceMatrix.setZeros();

  for (v = 0; v < numVertices; v++)
  {
    buffer[X] = i = vertices[v][X];
    buffer[Y] = j = vertices[v][Y];
    buffer[Z] = k = vertices[v][Z];
    value = voxelMatrix(i,j,k);
    buffer -= barycenter;
    covarianceMatrix +=  value * buffer.outerProduct(buffer);
    sumValues += value;
  }

  if ( sumValues == 0.0 )
  {
    return 0.0;
  }
  covarianceMatrix /= sumValues;

//  cout << "region size " << numVertices << " sumvalues = " << sumValues << endl;
//  barycenter.print("barycenter");
//  covarianceMatrix.print( "cov" );
  SquareMatrix<T> inverseCovarianceMatrix;
  try
  {
    inverseCovarianceMatrix = covarianceMatrix.inverse();
  }
  catch ( MathError /*unused*/ )
  {
    return 0.0;
  }

  for (v = 0; v < numVertices; v++)
  {
    buffer[X] = i = vertices[v][X];
    buffer[Y] = j = vertices[v][Y];
    buffer[Z] = k = vertices[v][Z];
    value = voxelMatrix(i,j,k);
    buffer -= barycenter;
    kurtosis += value * Maths::sqr( buffer.scalarProduct(inverseCovarianceMatrix*buffer) );
  }

  kurtosis /= sumValues;

  return kurtosis;
}

/*! This method clears regions with values (strictly) below a threshold.
 *
 * The values are given by \c regionValues. The vector given here
 * will generally be obtained by a previous call to #computeRegionFeature.
 * However, the method allows to provide values resulting from
 * arbitrary external computation.
 *
 * Regions below the threshold are set to 0 in the label matrix.
 * However, they remain in the list of regions (with a size and label of 0).
 * To erase them from the list, call again the #run method.
 *
 * \exception SizeError The length of \c regionValues differs from
 * the present number of regions.
****************************************************************/
template<class T>
void RegionAnalysis3D<T>::thresholdRegions(
  const Vector<T>& regionValues,
  const T threshold)
{
  if ( _regions.getSize() != regionValues.getSize() )
  {
    SizeError sizeError;
    sizeError.setWhere( "void RegionAnalysis3D<T>::thresholdRegions(...)" );
    sizeError.setWhat( "The number of values differs from the number of regions" );
    throw sizeError;
  }

  for (int r = 0; r < _regions.getSize(); ++r)
    if ( regionValues[r] < threshold )
    {
      fillRegion( _regions[r], 0, getLabelMatrix() );
      _regions[r].clear();
    }
}

/*! Returns the compactness (sphericity) of the specified region.
 *
 * Note that a rough estimation (over-estimation) of the surface
 * area is used.
****************************************************************/
template<class T>
float RegionAnalysis3D<T>::regionCompactness(const Region& region) const
{
  double c = static_cast<double>( region.getSize() ) / region._numFacets;
  c *= c;
  c /= region._numFacets;
  c *= 36.0 * M_PI;
  if ( c > 1.0 ) c = 0.0;

  return c;
}

template<class T>
float RegionAnalysis3D<T>::regionBoundaryMeasure(const Region& region) const
{
  const Vertices<int>& vertices = region.getVertices();
  const VoxelMatrix<T>& labelMatrix = getLabelMatrix();
  const int size1 = labelMatrix.getSize1();
  const int size2 = labelMatrix.getSize2();
  const int size3 = labelMatrix.getSize3();
  const float xVoxelSize = spatialCalibration()[0];
  const float yVoxelSize = spatialCalibration()[1];
  const float zVoxelSize = spatialCalibration()[2];
  const int r = region.getLabel();
  float boundaryMeasure = 0.0;
  int i, j, k, ii, jj, kk, n, s;

  Matrix<int> shifts( 6, 3 );
  shifts.setZeros();
  shifts(0,X) = -1;
  shifts(1,X) = 1;
  shifts(2,Y) = -1;
  shifts(3,Y) = 1;
  shifts(4,Z) = -1;
  shifts(5,Z) = 1;

  for (int v = 0; v < vertices.getNumVertices(); ++v)
  {
    i = vertices[v][X];
    j = vertices[v][Y];
    k = vertices[v][Z];

    for (n = 0; n < 6; n++)
    {
      ii = i + shifts(n,X);
      jj = j + shifts(n,Y);
      kk = k + shifts(n,Z);

      if ( ii >= 0 && ii < size1
        && jj >= 0 && jj < size2
        && kk >= 0 && kk < size3 )
      {
        s = static_cast<int>( labelMatrix(ii,jj,kk) );
        if ( s != r )
        {
          if ( ii != i ) boundaryMeasure += yVoxelSize * zVoxelSize;
          if ( jj != j ) boundaryMeasure += xVoxelSize * zVoxelSize;
          if ( kk != k ) boundaryMeasure += xVoxelSize * yVoxelSize;
        }
      }
    }
  }

  return boundaryMeasure;
}

#if 0
/*! Returns the volume (in physical units) of the specified region.
 *
 * Note that the spatial calibration of the label matrix is used.
****************************************************************/
template<class T>
float RegionAnalysis3D<T>::regionVolume(const Region& region) const
{
  return region.getSize() * _labelMatrix->getVoxelCalibration().getVoxelSize().prod();
}
#endif

#if 0














template<class T>
T RegionAnalysis<T>::barycenterValue(const Region& region) const
{
  Vectorf v = barycenter( region );
  const int i0 = int( v[X] );
  const int j0 = int( v[Y] );
  const int k0 = int( v[Z] );
  T result;

  if ( (*_labelMatrix)(i0,j0,k0) != region.label )
  {
    _barycenterErrorCount++;

    result = coreValue( region );
  if ( isnan(result) )
  {
    PRINT("A");
    EVAL(result);
    EVAL( v );
  }

  }
  else
  {
    int i, j, k, count = 0;
    result = 0.0;
    for (i = i0-2; i <= i0+2; i++)
    for (j = j0-2; j <= j0+2; j++)
    for (k = k0-2; k <= k0+2; k++)
    {
      if ( (*_labelMatrix)(i,j,k) == region.label )
      {
        result += (*_voxelMatrix)(i,j,k);
        count++;
      }
    }
    result /= count;

  if ( isnan(result) )
  {
    PRINT("B");
    EVAL(result);
    EVAL( v );
  }

  }


  return result;
}


#include <distancetransform.h>
#include "regionalextremafilter.h"

template<class T>
T RegionAnalysis<T>::coreValue(const Region& region) const
{
  const VoxelMatrix<T>& voxelMatrix = *_voxelMatrix;
  const VoxelMatrix<T>& labelMatrix = *_labelMatrix;
  const Verticesi& vertices = *region.vertices;
  const BoundingBoxi boundingBox = vertices.boundingBox();
  const Vectori& v1 = boundingBox.getVertex1();
  const Vectori& v2 = boundingBox.getVertex2();
  VoxelMatrix<T> boundingMatrix;

  boundingMatrix = labelMatrix.copy(
    v1[X], v2[X],
    v1[Y], v2[Y],
    v1[Z], v2[Z] );

  int i, j, k, count = 0;
#if 0
  for (i = 0; i < boundingMatrix.getSize1(); i++)
  for (j = 0; j < boundingMatrix.getSize2(); j++)
  for (k = 0; k < boundingMatrix.getSize3(); k++)
  {
    if ( boundingMatrix(i,j,k) == region.label )
    {
      count++;
    }
  }

  if ( count != region.getSize() )
  {
    ProgramError programError;
    programError.setWhere( "T RegionAnalysis<T>::coreValue(const Region&) const" );
    programError.setWhat( "Unexpected region size error" );
    EVAL( count );
    EVAL( region.getSize() );
    EVAL( boundingMatrix.getSize().prod() );
    throw programError;
  }
#endif

//   EVAL( count );


//   EVAL( boundingMatrix.range() );

  DistanceTransform<T> distanceTransform;
  distanceTransform.setForeground( region.label );
  distanceTransform.apply( boundingMatrix );

  VoxelMatrix<T> distanceMatrix = boundingMatrix;
//   EVAL( boundingMatrix.range() );
  RegionalExtremaFilter<T> regionalExtremaFilter;
  regionalExtremaFilter.apply( boundingMatrix );
  VoxelMatrix<T> maximaMatrix = boundingMatrix;
//   EVAL( boundingMatrix.range() );
  T sumValues = 0.0;

count = 0;
  for (i = 0; i < boundingMatrix.getSize1(); i++)
  for (j = 0; j < boundingMatrix.getSize2(); j++)
  for (k = 0; k < boundingMatrix.getSize3(); k++)
  {
    if ( boundingMatrix(i,j,k) )
    {
      sumValues += voxelMatrix(i+v1[X],j+v1[Y],k+v1[Z]);
      count++;
    }
  }

  if ( count == 0 )
  {
    EVAL( region.getSize() );
    EVAL( boundingMatrix.getSize().prod() );
    EVAL( v1 );
    EVAL( v2 );
    EVAL( distanceMatrix.range() );
    EVAL( maximaMatrix.range() );

    return 0.0;
  }

  return sumValues / count;
}

template<class T>
T RegionAnalysis<T>::distanceWeightedValue(const Region& region) const
{
  const VoxelMatrix<T>& voxelMatrix = *_voxelMatrix;
  const VoxelMatrix<T>& labelMatrix = *_labelMatrix;
  const Verticesi& vertices = *region.vertices;
  const BoundingBoxi boundingBox = vertices.boundingBox();
  const Vectori& v1 = boundingBox.getVertex1();
  const Vectori& v2 = boundingBox.getVertex2();
  VoxelMatrix<T> boundingMatrix;

  boundingMatrix = labelMatrix.copy(
    v1[X], v2[X],
    v1[Y], v2[Y],
    v1[Z], v2[Z] );

  int i, j, k, count = 0;
#if 0
  for (i = 0; i < boundingMatrix.getSize1(); i++)
  for (j = 0; j < boundingMatrix.getSize2(); j++)
  for (k = 0; k < boundingMatrix.getSize3(); k++)
  {
    if ( boundingMatrix(i,j,k) == region.label )
    {
      count++;
    }
  }

  if ( count != region.getSize() )
  {
    ProgramError programError;
    programError.setWhere( "T RegionAnalysis<T>::coreValue(const Region&) const" );
    programError.setWhat( "Unexpected region size error" );
    EVAL( count );
    EVAL( region.getSize() );
    EVAL( boundingMatrix.getSize().prod() );
    throw programError;
  }
#endif

//   EVAL( count );


//   EVAL( boundingMatrix.range() );

  DistanceTransform<T> distanceTransform;
  distanceTransform.setForeground( region.label );
  distanceTransform.apply( boundingMatrix );

  VoxelMatrix<T> distanceMatrix = boundingMatrix;
//   EVAL( boundingMatrix.range() );
/*  RegionalExtremaFilter<T> regionalExtremaFilter;
  regionalExtremaFilter.apply( boundingMatrix );*/
  VoxelMatrix<T> maximaMatrix = boundingMatrix;
//   EVAL( boundingMatrix.range() );
  T sumValues = 0.0;
  T sumWeights = 0.0;

count = 0;
  for (i = 0; i < boundingMatrix.getSize1(); i++)
  for (j = 0; j < boundingMatrix.getSize2(); j++)
  for (k = 0; k < boundingMatrix.getSize3(); k++)
  {
    if ( boundingMatrix(i,j,k) )
    {
      sumValues += boundingMatrix(i,j,k) * voxelMatrix(i+v1[X],j+v1[Y],k+v1[Z]);
      sumWeights += boundingMatrix(i,j,k);
      count++;
    }
  }

  if ( count == 0 )
  {
    EVAL( region.getSize() );
    EVAL( boundingMatrix.getSize().prod() );
    EVAL( v1 );
    EVAL( v2 );
    EVAL( distanceMatrix.range() );
    EVAL( maximaMatrix.range() );

    return 0.0;
  }

  return sumValues / sumWeights;
}


template<class T>
bool RegionAnalysis<T>::canMerge(const int r1, const int r2) const
{
  ENTER( "bool RegionAnalysis<T>::canMerge(const int, const int) const" );

  const Vectorf barycenter1 = barycenter( _regions[r1] );
  const Vectorf barycenter2 = barycenter( _regions[r2] );
  const int i1 = int( barycenter1[X] );
  const int j1 = int( barycenter1[Y] );
  const int k1 = int( barycenter1[Z] );
  const int i2 = int( barycenter2[X] );
  const int j2 = int( barycenter2[Y] );
  const int k2 = int( barycenter2[Z] );
  const int i = ( i1 + i2 ) / 2;
  const int j = ( j1 + j2 ) / 2;
  const int k = ( k1 + k2 ) / 2;

  LEAVE();

  return
    (*_voxelMatrix)(i,j,k) >= (*_voxelMatrix)(i1,j1,k1) ||
    (*_voxelMatrix)(i,j,k) >= (*_voxelMatrix)(i2,j2,k2);
}

class MergeRecord
{
  public:
    int region1;
    int region2;
    double score;

    bool operator<(const MergeRecord& mergeRecord) const
    {
      return score < mergeRecord.score;
    }
};

template<class T>
void RegionAnalysis<T>::mergeRegions(const unsigned int r, const unsigned int s)
{
  fillRegion( *_labelMatrix, _regions[s], _regions[r].label );
  _regions[r].vertices->append( *_regions[s].vertices );
  _regions[r].surface += _regions[s].surface - 2*_adjacencyGraph(r,s);
  _regions[s].clear();

  for (int t = 0; t < _regions.getSize(); t++)
  {
    if ( t != r && t != s )
    {
      _adjacencyGraph(r,t) += _adjacencyGraph(s,t);
      _adjacencyGraph(t,r) += _adjacencyGraph(t,s);
    }
    _adjacencyGraph(s,t) = 0;
    _adjacencyGraph(t,s) = 0;
  }
}

/*! Fusionne les régions voisines pour lesquelles la différence sur
 * \c regionFeature est inférieure ou égale au seuil \c threshold.
 *
 * La fusion est réalisée en commençant par fusionner les régions qui
 * diffèrent le moins.
 * Le processus est répété jusqu'à ce qu'il ne reste plus de régions
 * qui puissent être fusionnées.
*/
template<class T>
void RegionAnalysis<T>::mergeRegions2(
  RegionFeature regionFeature,
  const double threshold)
{
  ENTER( "void RegionAnalysis::mergeRegions2(...)" );

  const int numRegions = _regions.getSize();
  vector<MergeRecord> mergeRecordList;
  MergeRecord mergeRecord;
  int numMerges = 0;
  int r, s, t;
  double score;

//   Vectori labelMap( numRegions+1 );
//   for (r = 0; r <= numRegions; r++)
//   {
//     labelMap[r] = r;
//   }

  Vectord regionValues = computeRegionFeature( regionFeature );

  /*!
   * Initialisation de la liste des fusions
   */
  for (r = 0; r < numRegions-1; r++)
  {
    for (s = r+1; s < numRegions; s++)
    {
      if ( _adjacencyGraph(r,s) > 0 && canMergeByCompactness(r,s) )
      {
        score = fabs( regionValues[r] - regionValues[s] );
        if ( score <= threshold )
        {
          mergeRecord.region1 = r;
          mergeRecord.region2 = s;
          mergeRecord.score = score;
          mergeRecordList.push_back( mergeRecord );
        }
      }
    }
  }

  sort( mergeRecordList.begin(), mergeRecordList.end() );
  reverse( mergeRecordList.begin(), mergeRecordList.end() );

  /*!
   * Traitement de la liste des fusions
   */
  while ( mergeRecordList.empty() == false )
  {
    mergeRecord = mergeRecordList.back();
    r = mergeRecord.region1;
    s = mergeRecord.region2;
//     labelMap[s+1] = labelMap[r+1];
    numMerges++;

    /*!
     * Suppression de la région s au profit de la région r
     */
    fillRegion( *_labelMatrix, _regions[s], _regions[r].label );
    _regions[r].vertices->append( *_regions[s].vertices );
    _regions[r].surface += _regions[s].surface - 2*_adjacencyGraph(r,s);
    _regions[s].clear();

    for (t = 0; t < numRegions; t++)
    {
      if ( t != r && t != s )
      {
        _adjacencyGraph(r,t) += _adjacencyGraph(s,t);
        _adjacencyGraph(t,r) += _adjacencyGraph(t,s);
      }
      _adjacencyGraph(s,t) = 0;
      _adjacencyGraph(t,s) = 0;
    }

    regionValues[r] = computeRegionFeature( _regions[r], regionFeature );
    regionValues[s] = 0.0;

    /*!
     * Suppression de la liste de tous les items concernant
     * les deux régions (r et s) qui viennent d'être fusionnées
     */
    mergeRecordList.pop_back();
    for (int i = mergeRecordList.size()-1; i >= 0; i--)
    {
      mergeRecord = mergeRecordList[i];
      if ( mergeRecord.region1 == r || mergeRecord.region1 == s
        || mergeRecord.region2 == r || mergeRecord.region2 == s )
      {
        mergeRecordList.erase( mergeRecordList.begin() + i );
      }
    }

    for (s = 0; s < numRegions; s++)
    {
      if ( _adjacencyGraph(r,s) > 0 && canMergeByCompactness(r,s) )
      {
        score = fabs( regionValues[r] - regionValues[s] );
        if ( score <= threshold )
        {
          mergeRecord.region1 = s > r? r: s; //! region1: indice le plus petit
          mergeRecord.region2 = s > r? s: r; //! region2: indice le plus grand
          mergeRecord.score = score;
          mergeRecordList.push_back( mergeRecord );
        }
      }
    }

    sort( mergeRecordList.begin(), mergeRecordList.end() );
    reverse( mergeRecordList.begin(), mergeRecordList.end() );
  }

  /*!
   * Résolution des équivalences et
   * condensation de la table des labels
   */
//   int numEquivs = 0;
//   for (r = 1; r <= numRegions; r++)
//   {
//     if ( labelMap[r] == r )
//     {
//       labelMap[r] -= numEquivs;
//     }
//     else
//     {
//       labelMap[r] = labelMap[labelMap[r]];
//       numEquivs++;
//     }
// //     _regions[r-1].label = labelMap[r]; ATTN
//   }

//   applyLabelMap( labelMap ); ATTN

//   Vectorui histogram = _labelMatrix->histogram( 0, 1, numRegions+1 );
//   EVAL( _labelMatrix->range() );
//   for (r = 0; r < numRegions; r++)
//   {
//     EVAL( r );
//     EVAL( histogram[_regions[r].label] );
//     EVAL( _regions[r].getSize() );
//   }

  EVAL( numMerges );
  LEAVE();
}

template<class T>
void RegionAnalysis<T>::mergeRegions(RegionFeature regionFeature)
{
  ENTER( "void RegionAnalysis::mergeRegions(RegionFeature)" );

  const int numRegions = _regions.getSize();
  vector<MergeRecord> mergeRecordList;
  MergeRecord mergeRecord;
  int r, s, t;
  T diff;

  Vectord regionValues = computeRegionFeature( regionFeature );

  /*!
   * Initialisation de la liste des fusions
   */
  for (r = 0; r < numRegions-1; r++)
  {
    for (s = r+1; s < numRegions; s++)
    {
      if ( _adjacencyGraph(r,s) > 0 && canMergeByCompactness(r,s) )
      {
        diff = fabs( regionValues[r] - regionValues[s] );
        mergeRecord.region1 = r;
        mergeRecord.region2 = s;
        mergeRecord.score = diff;
        mergeRecordList.push_back( mergeRecord );
      }
    }
  }

  sort( mergeRecordList.begin(), mergeRecordList.end() );
  reverse( mergeRecordList.begin(), mergeRecordList.end() );

  /*!
   * Traitement de la liste des fusions
   */
  while ( mergeRecordList.empty() == false )
  {
    mergeRecord = mergeRecordList.back();
    r = mergeRecord.region1;
    s = mergeRecord.region2;

    /*!
     * Suppression de la région s au profit de la région r
     */
    mergeRegions( r, s );
    regionValues[r] = computeRegionFeature( _regions[r], regionFeature );
    regionValues[s] = 0.0;

    /*!
     * Suppression de la liste de tous les items concernant
     * les deux régions (r et s) qui viennent d'être fusionnées
     */
    mergeRecordList.pop_back();
    for (int i = mergeRecordList.size()-1; i >= 0; i--)
    {
      mergeRecord = mergeRecordList[i];
      if ( mergeRecord.region1 == r || mergeRecord.region1 == s
        || mergeRecord.region2 == r || mergeRecord.region2 == s )
      {
        mergeRecordList.erase( mergeRecordList.begin() + i );
      }
    }

    for (s = 0; s < numRegions; s++)
    {
      if ( _adjacencyGraph(r,s) > 0 && canMergeByCompactness(r,s) )
      {
        diff = fabs( regionValues[r] - regionValues[s] );
        mergeRecord.region1 = s > r? r: s;
        mergeRecord.region2 = s > r? s: r;
        mergeRecord.score = diff;
        mergeRecordList.push_back( mergeRecord );
      }
    }

    sort( mergeRecordList.begin(), mergeRecordList.end() );
    reverse( mergeRecordList.begin(), mergeRecordList.end() );
  }

  LEAVE();
}


template<class T>
void RegionAnalysis<T>::reconstructionClosing(
  Vectord& regionValues,
  const unsigned int size)
{
  ENTER( "void RegionAnalysis<T>::reconstructionClosing(...)" );

  const int numRegions = _regions.getSize();
  Vectord reconstruction( numRegions );
  Vectord dilation( numRegions );
  Vectord erosion( numRegions );
  int numChanges;
  int r, s;

  for (r = 0; r < numRegions; r++)
  {
    reconstruction[r] = regionValues[r];
  }

  // dilatation (résultat stocké dans reconstruction)
  for (unsigned int k = 0; k < size; k++)
  {
    for (r = 0; r < numRegions; r++)
    {
      dilation[r] = reconstruction[r];
      for (s = 0; s < numRegions; s++)
      {
        if ( _adjacencyGraph(r,s) > 0 )
        {
          if ( reconstruction[s] > dilation[r] )
          {
            dilation[r] = reconstruction[s];
          }
        }
      }
    }
    reconstruction = dilation;
  }

  // reconstruction par érosion
  numChanges = 1;

  while ( numChanges > 0 )
  {
    // érosion
    for (r = 0; r < numRegions; r++)
    {
      erosion[r] = reconstruction[r];
      for (s = 0; s < numRegions; s++)
      {
        if ( _adjacencyGraph(r,s) > 0 )
        {
          if ( reconstruction[s] < erosion[r] )
          {
            erosion[r] = reconstruction[s];
          }
        }
      }
    }

    numChanges = 0;

    // application du masque (maximum point à point)
    for (r = 0; r < numRegions; r++)
    {
      if ( erosion[r] < regionValues[r] )
      {
        if ( reconstruction[r] != regionValues[r] )
        {
          reconstruction[r] = regionValues[r];
          numChanges++;
        }
      }
      else if ( erosion[r] > regionValues[r] )
      {
        if ( reconstruction[r] != erosion[r] )
        {
          reconstruction[r] = erosion[r];
          numChanges++;
        }
      }
    }
  }

  for (r = 0; r < numRegions; r++)
  {
    regionValues[r] = reconstruction[r];
  }

  LEAVE();
}

template<class T>
Vectord RegionAnalysis<T>::computeContrasts(
  const Vectord& regionValues,
  const bool linear) const
{
  ENTER( "Vectord RegionAnalysis<T>::computeContrasts(...)" );

  const int numRegions = _regions.getSize();
  Vectord contrasts( numRegions );
  double contrast;
  double positiveContrast;
  double negativeContrast;
  int weight, sumWeights;
  int sumPositiveWeights;
  int sumNegativeWeights;
  int r, s;
  T diff;

  double minValue = regionValues.min();
  EVAL( minValue );

  if ( linear )
  {
    for (r = 0; r < numRegions; r++)
    {
      if ( _regions[r].getSize() == 0 )
      {
        contrasts[r] = 0.0;
        continue;
      }

      contrast = 0.0;
      sumWeights = 0;
//       int n = 0;
//
//       for (s = 0; s < numRegions; s++)
//       {
//         if ( _adjacencyGraph(r,s) > 0 )
//         {
//           n++;
//         }
//       }
//
//       Vectord values( n );
//       n = 0;
//       for (s = 0; s < numRegions; s++)
//       {
//         if ( _adjacencyGraph(r,s) > 0 )
//         {
//           values[n++] = regionValues[s];
//         }
//       }
//
//       double median = values.median();

      for (s = 0; s < numRegions; s++)
      {
        if ( _adjacencyGraph(r,s) > 0 && regionValues[s] > minValue )
        {
          int radius = pow( _regions[s].getSize(), 1.0/3.0 );
          weight = radius * _adjacencyGraph(r,s);
          diff = regionValues[r] - regionValues[s];
          contrast += weight * diff;
          sumWeights += weight;
        }
      }

      contrasts[r] = sumWeights > 0? contrast / sumWeights: regionValues[r]-minValue;

#if 0
      for (s = 0; s < numRegions; s++)
      {
        if ( _adjacencyGraph(r,s) > 0 )
        {
          weight = _regions[s].getSize();
          diff = regionValues[r] - regionValues[s];
          contrast += weight * diff;
          sumWeights += weight;
        }
      }
      contrasts[r] = sumWeights > 0? contrast / sumWeights: 0.0;
#endif
    }
  }
  else
  {
#if 0
    for (r = 0; r < numRegions; r++)
    {
      positiveContrast = 0.0;
      negativeContrast = 0.0;
      sumPositiveWeights = 0;
      sumNegativeWeights = 0;

      for (s = 0; s < numRegions; s++)
      {
        weight = _adjacencyGraph(r,s);
        if ( weight > 0 )
        {
          weight = 1;
          diff = _regions[r].feature( regionFeature ) - _regions[s].feature( regionFeature );
          if ( diff < 0 )
          {
            negativeContrast += weight * diff * diff;
            sumNegativeWeights += weight;
          }
          else
          {
            positiveContrast += weight * diff *diff;
            sumPositiveWeights += weight;
          }
        }
      }
      positiveContrast *= sumPositiveWeights;
      negativeContrast *= sumNegativeWeights;
      sumWeights = sumPositiveWeights + sumNegativeWeights;
      if ( sumWeights > 0 )
      {
        contrast = sqrt( positiveContrast ) - sqrt( negativeContrast );
        contrast /= sumWeights;
      }
      else
      {
        contrast = 0.0;
      }
      _regions[r].relativeContrast = contrast;
    }
#else
/*    for (r = 0; r < numRegions; r++)
    {
      contrast = 0.0;
      sumWeights = 0;

      for (s = 0; s < numRegions; s++)
      {
        weight = _adjacencyGraph(r,s);
        if ( weight > 0 )
        {
          weight = 1;
          diff =
            _regions[r].feature( regionFeature ) * _regions[r].feature( regionFeature ) -
            _regions[s].feature( regionFeature ) * _regions[s].feature( regionFeature );
          contrast += weight * diff;
          sumWeights += weight;
        }
      }
      if ( sumWeights > 0 )
      {
        contrast /= sumWeights;
      }
      _regions[r].relativeContrast = contrast;
    }*/
#endif
  }

  LEAVE();

  return contrasts;
}


template<class T>
void RegionAnalysis<T>::remapRegionLabels()
{
  int numRegions = 0;
  
  _labelMatrix->fill( 0.0 );
  
  for (unsigned int r = 0; r < _regions.getSize(); r++)
  {
    if ( _regions[r].getSize() > 0 )
    {
      _regions[r].label = ++numRegions;
      fillRegion( *_labelMatrix, _regions[r], _regions[r].label );
    }
  }

  EVAL( numRegions );
}













#if 0
template<class T>
void RegionAnalysis<T>::mapFeature(
  VoxelMatrix<T>& voxelMatrix,
  RegionFeature regionFeature) const
{
  const VoxelMatrix<T>& labels = *_labelMatrix;
  const int size1 = labels.getSize1();
  const int size2 = labels.getSize2();
  const int size3 = labels.getSize3();
  int i, j, k, r;

  Vectorf barycenterValues( _regions.getSize() );
  if ( regionFeature == REGION_FEATURE_BARYCENTER_VALUE )
  {
    VoxelMatrix<T> averageMatrix = *_voxelMatrix;
    SumFilter<T> sumFilter;
    sumFilter.apply( averageMatrix );
    averageMatrix /= 27;

    for (r = 0; r < _regions.getSize(); r++)
    {
      Vectorf b = barycenter( _regions[r] );
      barycenterValues[r] = averageMatrix( int(b[X]), int(b[Y]), int(b[Z]) );
    }
  }

  voxelMatrix.setSize( labels.getSize() );

  for (i = 0; i < size1; i++)
  for (j = 0; j < size2; j++)
  for (k = 0; k < size3; k++)
  {
    r = int( labels(i,j,k) );
    if ( r > 0 )
    {
      if ( regionFeature == REGION_FEATURE_BARYCENTER_VALUE )
      {
        voxelMatrix(i,j,k) = barycenterValues[r-1];
      }
      else
      {
        voxelMatrix(i,j,k) = _regions[r-1].feature( regionFeature );
      }
    }
    else
    {
      voxelMatrix(i,j,k) = 0.0;
    }
  }
}
#endif

template<class T>
float RegionAnalysis<T>::averageNumNeighbours() const
{
  const int numRegions = _regions.getSize();
  int sum = 0;
  int r, s;

  for (r = 0; r < numRegions; r++)
  {
    for (s = 0; s < numRegions; s++)
    {
      if ( _adjacencyGraph(r,s) > 0 )
      {
        sum++;
      }
    }
  }

  return float( sum ) / numRegions;
}


template<class T>
const Matrixi& RegionAnalysis<T>::getAdjacencyGraph() const
{
  return _adjacencyGraph;
}

template<class T>
Matrixd RegionAnalysis<T>::regionMatrix() const
{
  int n, r, numRegions = _regions.getSize();
  Matrixd matrix( numRegions, 7 );
  Vectorf center;

  for (r = 0; r < numRegions; r++)
  {
    const Region& region = _regions[r];
    n = 0;
    center = barycenter( region );
    matrix(r,n++) = r + 1;
    matrix(r,n++) = center[X];
    matrix(r,n++) = center[Y];
    matrix(r,n++) = center[Z];
    matrix(r,n++) = region.getSize();
    matrix(r,n++) = compactness( region );
    matrix(r,n++) = barycenterValue( region );
  }

  return matrix;
}

template<class T>
string RegionAnalysis<T>::regionMatrixHeader() const
{
  string header = string( "Label\tCenterX\tCenterY\tCenterZ\tSize\tCompactness"
    "\tAverage" );

  if ( _stackName.empty() == false )
  {
    header = "Stack\t" + header;
  }

  return header;
}

template<class T>
void RegionAnalysis<T>::saveRegionMatrix(const string& fileName) const
{
  saveRegionMatrix( regionMatrix(), fileName );
}

template<class T>
void RegionAnalysis<T>::saveRegionMatrix(const Matrixd& regionMatrix, const string& fileName) const
{
  ofstream ofs( fileName.c_str() );
  ofs << regionMatrixHeader() << endl;
  appendRegionMatrix( regionMatrix, ofs );
}

template<class T>
void RegionAnalysis<T>::appendRegionMatrix(ofstream& ofs) const
{
//   Matrixd matrix = regionMatrix();
  appendRegionMatrix( regionMatrix(), ofs );
/*  for (unsigned int r = 0; r < matrix.getSize1(); r++)
  {
    if ( _stackName.empty() == false )
    {
      ofs << _stackName << " ";
    }
    ofs << fixed << setprecision(2) << setw(8) << matrix[r] << endl;
  }*/
}

template<class T>
void RegionAnalysis<T>::appendRegionMatrix(const Matrixd& regionMatrix, ofstream& ofs) const
{
  for (unsigned int r = 0; r < regionMatrix.getSize1(); r++)
  {
    if ( _stackName.empty() == false )
    {
      ofs << _stackName << " ";
    }
    ofs << fixed << setprecision(2) << setw(8) << regionMatrix[r] << endl;
  }
}

template<class T>
void RegionAnalysis<T>::overlayContours(VoxelMatrix<T>& volumeMatrix) const
{
  int numRegions = _regions.getSize();
  int numVertices;
  int r, v, i, j, k;

  for (r = 0; r < numRegions; r++)
  {
    const Verticesi& vertices = *_regions[r].vertices;
    numVertices = vertices.getSize();
    for (v = 0; v < numVertices; v++)
    {
      if ( _regions[r].edgeFlags[v] )
      {
        i = vertices[v][X];
        j = vertices[v][Y];
        k = vertices[v][Z];
        volumeMatrix(i,j,k) = 255.0;
      }
    }
  }
}



template<class T>
Vector<T> RegionAnalysis<T>::regionValues(const Region& region) const
{
  Vector<T> values( region.size );
  const VoxelMatrix<T>& volume = *_voxelMatrix;
  const Verticesi& vertices = *region.vertices;
  int i, j, k, v;

  for (v = 0; v < region.size; v++)
  {
    i = vertices[v][X];
    j = vertices[v][Y];
    k = vertices[v][Z];
    values[v] = volume(i,j,k);
  }

  return values;
}

template<class T>
Vectorui RegionAnalysis<T>::computeRegionHistogram(
  const Region& region,
  const T binWidth) const
{
  return regionValues( region ).histogram( 0, binWidth, 256 );
}

template<class T>
void RegionAnalysis<T>::tophat()
{
  ENTER( "void RegionAnalysis<T>::tophat()" );

  const int numRegions = _regions.getSize();
  Vectord minima( numRegions );
  double maximum;
  int r, s;

  for (r = 0; r < numRegions; r++)
  {
    minima[r] = _regions[r].averageValue;
    for (s = 0; s < numRegions; s++)
    {
      if ( _adjacencyGraph(r,s) > 0 )
      {
        if ( _regions[s].averageValue < minima[r] )
        {
          minima[r] = _regions[s].averageValue;
        }
      }
    }
  }

  for (r = 0; r < numRegions; r++)
  {
    maximum = minima[r];
    for (s = 0; s < numRegions; s++)
    {
      if ( _adjacencyGraph(r,s) > 0 )
      {
        if ( minima[s] > maximum )
        {
          maximum = minima[s];
        }
      }
    }
    _regions[r].averageValue -= maximum;

    if ( _regions[r].averageValue < 0 )
    {
      PRINT( "Error: negative average" );
    }
  }

  LEAVE();
}

#if 0
template<class T>
void RegionAnalysis<T>::reconstructionClosing(const unsigned int size)
{
  ENTER( "void RegionAnalysis<T>::reconstructionClosing()" );

  const int numRegions = _regions.getSize();
  Vectord reconstruction( numRegions );
  Vectord dilation( numRegions );
  Vectord erosion( numRegions );
  int numChanges;
  int r, s;

  for (r = 0; r < numRegions; r++)
  {
    reconstruction[r] = _regions[r].averageValue;
  }

  // dilatation (résultat stocké dans reconstruction)
  for (unsigned int k = 0; k < size; k++)
  {
    for (r = 0; r < numRegions; r++)
    {
      dilation[r] = reconstruction[r];
      for (s = 0; s < numRegions; s++)
      {
        if ( _adjacencyGraph(r,s) > 0 )
        {
          if ( reconstruction[s] > dilation[r] )
          {
            dilation[r] = reconstruction[s];
          }
        }
      }
    }
    reconstruction = dilation;
  }

  // reconstruction par érosion
  numChanges = 1;

  while ( numChanges > 0 )
  {
    // érosion
    for (r = 0; r < numRegions; r++)
    {
      erosion[r] = reconstruction[r];
      for (s = 0; s < numRegions; s++)
      {
        if ( _adjacencyGraph(r,s) > 0 )
        {
          if ( reconstruction[s] < erosion[r] )
          {
            erosion[r] = reconstruction[s];
          }
        }
      }
    }

    numChanges = 0;

    // application du masque (maximum point à point)
    for (r = 0; r < numRegions; r++)
    {
      if ( erosion[r] < _regions[r].averageValue )
      {
        if ( reconstruction[r] != _regions[r].averageValue )
        {
          reconstruction[r] = _regions[r].averageValue;
          numChanges++;
        }
      }
      else if ( erosion[r] > _regions[r].averageValue )
      {
        if ( reconstruction[r] != erosion[r] )
        {
          reconstruction[r] = erosion[r];
          numChanges++;
        }
      }
    }
  }

  for (r = 0; r < numRegions; r++)
  {
    _regions[r].averageValue = reconstruction[r];
  }


  LEAVE();
}
#endif

template<class T>
void RegionAnalysis<T>::enhanceContrast(const double sigma)
{
  const int numRegions = _regions.getSize();
  Vectord variations( numRegions );
  double diff, weight;
  int r, s;

  variations.setZeros();

  for (r = 0; r < numRegions; r++)
  {
    for (s = r+1; s < numRegions; s++)
    {
      if ( _adjacencyGraph(r,s) > 0 )
      {
        diff = _regions[r].averageValue - _regions[s].averageValue;
        weight = diff / sigma;
        weight *= weight;
        weight = ( 1.0 - weight ) * exp( -0.5 * weight ) / (sigma*sigma);
        variations[r] -= weight * diff;
        variations[s] += weight * diff;
      }
    }
  }

  for (r = 0; r < numRegions; r++)
  {
    _regions[r].averageValue += variations[r];
  }
}

template<class T>
void RegionAnalysis<T>::bilateralFiltering(const double sigma)
{
  const int numRegions = _regions.getSize();
  Vectord numerators( numRegions );
  Vectord denominators( numRegions );
  double weight, delta;
  int r, s;

  for (r = 0; r < numRegions; r++)
  {
    numerators[r] = _regions[r].averageValue;
    denominators[r] = 1.0;
  }

  for (r = 0; r < numRegions; r++)
  {
    for (s = r+1; s < numRegions; s++)
    {
      if ( _adjacencyGraph(r,s) > 0 )
      {
        delta = ( _regions[r].averageValue - _regions[s].averageValue ) / sigma;
        weight = exp( -0.5 * delta * delta );
        numerators[r] += weight * _regions[s].averageValue;
        numerators[s] += weight * _regions[r].averageValue;
        denominators[r] += weight;
        denominators[s] += weight;
      }
    }
  }

  for (r = 0; r < numRegions; r++)
  {
    _regions[r].averageValue = numerators[r] / denominators[r];
  }
}

template<class T>
void RegionAnalysis<T>::boxBilateralFiltering(const double sigma)
{
  const int numRegions = _regions.getSize();
  Vectord numerators( numRegions );
  Vectord denominators( numRegions );
  bool stability = false;
  double delta;
  int r, s;

  while ( stability == false )
  {
    for (r = 0; r < numRegions; r++)
    {
      numerators[r] = _regions[r].averageValue;
      denominators[r] = 1.0;
    }

    for (r = 0; r < numRegions; r++)
    {
      for (s = r+1; s < numRegions; s++)
      {
        if ( _adjacencyGraph(r,s) > 0 )
        {
          delta = fabs( _regions[r].averageValue - _regions[s].averageValue );
          if ( delta < sigma )
          {
            numerators[r] += _regions[s].averageValue;
            numerators[s] += _regions[r].averageValue;
            denominators[r]++;
            denominators[s]++;
          }
        }
      }
    }

    stability = true;
    for (r = 0; r < numRegions; r++)
    {
      double average = numerators[r] / denominators[r];
      if ( fabs(average-_regions[r].averageValue) > 1.0 )
      {
        stability = false;
      }
      _regions[r].averageValue = average;
    }
  }

//   EVAL( niters );

}

template<class T>
void RegionAnalysis<T>::computeRelativeContrasts(
  RegionFeature regionFeature,
  const bool linear)
{
  const int numRegions = _regions.getSize();
  double contrast;
  double positiveContrast;
  double negativeContrast;
  int weight, sumWeights;
  int sumPositiveWeights;
  int sumNegativeWeights;
  int r, s;
  T diff;

  if ( linear )
  {
    for (r = 0; r < numRegions; r++)
    {
      contrast = 0.0;
      sumWeights = 0;

      for (s = 0; s < numRegions; s++)
      {
        weight = _adjacencyGraph(r,s);
        if ( weight > 0 )
        {
          weight = _regions[s].size;
          diff = _regions[r].feature( regionFeature ) - _regions[s].feature( regionFeature );
          contrast += weight * diff;
          sumWeights += weight;
        }
      }
      _regions[r].relativeContrast = sumWeights > 0? contrast / sumWeights: 0.0;
    }
  }
  else
  {
#if 0
    for (r = 0; r < numRegions; r++)
    {
      positiveContrast = 0.0;
      negativeContrast = 0.0;
      sumPositiveWeights = 0;
      sumNegativeWeights = 0;

      for (s = 0; s < numRegions; s++)
      {
        weight = _adjacencyGraph(r,s);
        if ( weight > 0 )
        {
          weight = 1;
          diff = _regions[r].feature( regionFeature ) - _regions[s].feature( regionFeature );
          if ( diff < 0 )
          {
            negativeContrast += weight * diff * diff;
            sumNegativeWeights += weight;
          }
          else
          {
            positiveContrast += weight * diff *diff;
            sumPositiveWeights += weight;
          }
        }
      }
      positiveContrast *= sumPositiveWeights;
      negativeContrast *= sumNegativeWeights;
      sumWeights = sumPositiveWeights + sumNegativeWeights;
      if ( sumWeights > 0 )
      {
        contrast = sqrt( positiveContrast ) - sqrt( negativeContrast );
        contrast /= sumWeights;
      }
      else
      {
        contrast = 0.0;
      }
      _regions[r].relativeContrast = contrast;
    }
#else
    for (r = 0; r < numRegions; r++)
    {
      contrast = 0.0;
      sumWeights = 0;

      for (s = 0; s < numRegions; s++)
      {
        weight = _adjacencyGraph(r,s);
        if ( weight > 0 )
        {
          weight = 1;
          diff =
            _regions[r].feature( regionFeature ) * _regions[r].feature( regionFeature ) -
            _regions[s].feature( regionFeature ) * _regions[s].feature( regionFeature );
          contrast += weight * diff;
          sumWeights += weight;
        }
      }
      if ( sumWeights > 0 )
      {
        contrast /= sumWeights;
      }
      _regions[r].relativeContrast = contrast;
    }
#endif
  }
}

template<class T>
int RegionAnalysis<T>::removeNegativeContrastRegions()
{
  ENTER( "int RegionAnalysis<T>::removeNegativeContrastRegions()" );

  const int numRegions = _regions.getSize();
  Vectori labelMap( numRegions+1 );
  int r, s, numRemoved, numRetained = 0;

  labelMap[0] = 0;

  for (r = 0; r < numRegions; r++)
  {
    if ( _regions[r].relativeContrast > 0 )
    {
      labelMap[r+1] = ++numRetained;
    }
    else
    {
      labelMap[r+1] = 0;
    }
  }

  numRemoved = numRegions - numRetained;

  if ( numRemoved > 0 )
  {
    Vector<Region> regions( numRetained );
    for (r = s = 0; r < numRegions; r++)
    {
      if ( labelMap[r+1] )
      {
        regions[s++] = _regions[r];
      }
    }
    _regions.swap( regions );

    applyLabelMap( labelMap );
  }

  LEAVE();

  return numRemoved;
}

template<class T>
int RegionAnalysis<T>::thresholdRegions(
  RegionFeature regionFeature,
  const double featureThreshold)
{
  ENTER( "int RegionAnalysis<T>::thresholdRegions(...)" );

  const int numRegions = _regions.getSize();
  Vectori labelMap( numRegions+1 );
  int r, s, numRemoved, numRetained = 0;

  labelMap[0] = 0;

  for (r = 0; r < numRegions; r++)
  {
    if ( _regions[r].feature(regionFeature) < featureThreshold )
    {
      labelMap[r+1] = 0;
    }
    else
    {
      labelMap[r+1] = ++numRetained;
    }
  }

  numRemoved = numRegions - numRetained;

  if ( numRemoved > 0 )
  {
    Vector<Region> regions( numRetained );
    for (r = s = 0; r < numRegions; r++)
    {
      if ( labelMap[r+1] )
      {
        regions[s++] = _regions[r];
      }
    }
    _regions.swap( regions );

    applyLabelMap( labelMap );
  }

  EVAL( numRemoved );
  LEAVE();

  return numRemoved;
}

template<class T>
bool RegionAnalysis<T>::canMergeByCompactness(const int r1, const int r2) const
{
  double compactness, compactness1, compactness2;

  return canMergeByCompactness( r1, r2, compactness1, compactness2, compactness );
}

template<class T>
bool RegionAnalysis<T>::canMergeByCompactness(
  const int r1, const int r2,
  double& compactness1, double& compactness2,
  double& compactness) const
{
  double S1 = _regions[r1].surface;
  double S2 = _regions[r2].surface;
  double V1 = _regions[r1].size;
  double V2 = _regions[r2].size;
  double V = V1 + V2;
  double S = S1 + S2 - 2.0 * _adjacencyGraph(r1,r2);
  compactness1 = (V1/S1) * (V1/S1) / S1;
  compactness2 = (V2/S2) * (V2/S2) / S2;
  compactness = (V/S) * (V/S) / S;

  return V*compactness >= V1*compactness1+V2*compactness2;
}

template<class T>
int RegionAnalysis<T>::mergeRegions()
{
  ENTER( "void RegionAnalysis::mergeRegions()" );

  const int numRegions = _regions.getSize();
  Vectori matchTable( numRegions );
  int r, s, t, bestMatch;
  double sph1, sph2;
  double sph, bestSph;
  double diff, minDiff;
//   double tolerance = 3.0;
  bool firstMatch;

  int numMatches, totalNumMerges = 0;

  Vectori labelMap( numRegions+1 );
  for (r = 0; r <= numRegions; r++)
  {
    labelMap[r] = r;
  }

  do
  {
    /*!
     * Détermination des appariements candidats
     */
    for (r = 0; r < numRegions; r++)
    {
      bestMatch = r;
      firstMatch = true;
      minDiff = 0.0;
      for (s = 0; s < numRegions; s++)
      {
        if ( _adjacencyGraph(r,s) > 0 )
        {
          if ( r == s )
          {
            PRINT( "-------> Adjacency Graph Error" );
          }
          if ( _regions[r].size <= 0 || _regions[s].size <= 0 )
          {
            PRINT( "-------> Region Size Error" );
          }
          diff = fabs( _regions[r].averageValue - _regions[s].averageValue );
//           if ( diff <= tolerance )// && ( firstMatch || diff < minDiff ) )
          {
            bool ok = canMergeByCompactness( r, s, sph1, sph2, sph );
            if ( ok && ( firstMatch || sph > bestSph ) )
            {
              bestMatch = s;
              firstMatch = false;
              minDiff = diff;
              bestSph = sph;
            }
          }
        }
      }

      matchTable[r] = bestMatch;
    }

    /*!
     * Suppression des appariements non-bijectifs,
     * détermination du nombre d'appariements effectifs
     * et construction de la table de mise à jour des
     * labels.
     */
    numMatches = 0;

    for (r = 0; r < numRegions; r++)
    {
      s = matchTable[r];

      if ( _regions[r].size <= 0 && r != matchTable[r] )
      {
        PRINT( "-------> Match table error" );
      }
      if ( s < r && matchTable[s] == r )
      {
        numMatches++;
        labelMap[r+1] = labelMap[s+1];
        if ( labelMap[s+1] != s+1 )
        {
          PRINT( "-------> Match error: already matched region" );
        }

        // la région r disparait au profit de s
        _regions[s].averageValue *= _regions[s].size;
        _regions[s].averageValue += _regions[r].size * _regions[r].averageValue;
        _regions[s].size += _regions[r].size;
        _regions[s].averageValue /= _regions[s].size;
        _regions[s].vertices->append( *_regions[r].vertices );
        _regions[s].surface += _regions[r].surface - 2*_adjacencyGraph(r,s);
        _regions[r].size = 0;
        _regions[r].surface = 0;
        for (t = 0; t < numRegions; t++)
        {
          if ( t != s && t != r )
          {
            _adjacencyGraph(s,t) += _adjacencyGraph(r,t);
            _adjacencyGraph(t,s) += _adjacencyGraph(t,r);
          }
          _adjacencyGraph(r,t) = 0;
          _adjacencyGraph(t,r) = 0;
        }
      }
    }

    EVAL( numMatches );
    totalNumMerges += numMatches;
  } while ( numMatches > 0 );

  /*!
   * Résolution des équivalences et
   * condensation de la table des labels
   */
  int numEquivs = 0;
  for (r = 1; r <= numRegions; r++)
  {
    if ( labelMap[r] == r )
    {
      labelMap[r] -= numEquivs;
    }
    else
    {
      labelMap[r] = labelMap[labelMap[r]];
      numEquivs++;
    }
  }

  applyLabelMap( labelMap );

  EVAL( totalNumMerges );
  LEAVE();

  return totalNumMerges;
}

template<class T>
int RegionAnalysis<T>::mergeRegionsByPriority()
{
  ENTER( "void RegionAnalysis::mergeRegionsByPriority()" );

  const int numRegions = _regions.getSize();
  int r, s, t, bestRegion1, bestRegion2;
  double compactness1, compactness2;
  double compactness; //, bestCompactness;
  bool matchFound;
  int numMerges = 0;

  T diff, bestDiff;

  Vectori labelMap( numRegions+1 );
  for (r = 0; r <= numRegions; r++)
  {
    labelMap[r] = r;
  }

  do
  {
    matchFound = false;
    for (r = 0; r < numRegions-1; r++)
    {
      for (s = r+1; s < numRegions; s++)
      {
        if ( _adjacencyGraph(r,s) > 0 && canMergeByCompactness(r,s,compactness1,compactness2,compactness) )
        {
#if 1
          diff = fabs( _regions[r].averageValue - _regions[s].averageValue );
          if ( !matchFound || diff < bestDiff )
          {
            bestRegion1 = r;
            bestRegion2 = s;
            bestDiff = diff;
            matchFound = true;
          }
#else
          if ( !matchFound || compactness > bestCompactness )
          {
            matchFound = true;
            bestCompactness = compactness;
            bestRegion1 = r;
            bestRegion2 = s;
          }
#endif
        }
      }
    }

    if ( matchFound )
    {
      r = bestRegion2;
      s = bestRegion1;
      labelMap[r+1] = labelMap[s+1];
      numMerges++;

      // la région r disparait au profit de s
      _regions[s].averageValue *= _regions[s].size;
      _regions[s].averageValue += _regions[r].size * _regions[r].averageValue;
      _regions[s].size += _regions[r].size;
      _regions[s].averageValue /= _regions[s].size;
      _regions[s].vertices->append( *_regions[r].vertices );
      _regions[s].surface += _regions[r].surface - 2*_adjacencyGraph(r,s);
      _regions[r].size = 0;
      _regions[r].surface = 0;
      for (t = 0; t < numRegions; t++)
      {
        if ( t != s && t != r )
        {
          _adjacencyGraph(s,t) += _adjacencyGraph(r,t);
          _adjacencyGraph(t,s) += _adjacencyGraph(t,r);
        }
        _adjacencyGraph(r,t) = 0;
        _adjacencyGraph(t,r) = 0;
      }
    }
  } while ( matchFound );

  /*!
   * Résolution des équivalences et
   * condensation de la table des labels
   */
  int numEquivs = 0;
  for (r = 1; r <= numRegions; r++)
  {
    if ( labelMap[r] == r )
    {
      labelMap[r] -= numEquivs;
    }
    else
    {
      labelMap[r] = labelMap[labelMap[r]];
      numEquivs++;
    }
  }

  applyLabelMap( labelMap );

  EVAL( numMerges );
  LEAVE();

  return numMerges;
}


template<class T>
int RegionAnalysis<T>::mergeRegionsByPriority2()
{
  ENTER( "void RegionAnalysis::mergeRegionsByPriority2()" );

  vector<MergeRecord> mergeRecordList;
  const int numRegions = _regions.getSize();
  int r, s, t;
  double compactness1, compactness2;
  double compactness;
  int numMerges = 0;

  MergeRecord mergeRecord;

  T diff;

  Vectori labelMap( numRegions+1 );
  for (r = 0; r <= numRegions; r++)
  {
    labelMap[r] = r;
  }

  /*!
   * Initialisation de la liste des fusions
   */
  for (r = 0; r < numRegions-1; r++)
  {
    for (s = r+1; s < numRegions; s++)
    {
      if ( _adjacencyGraph(r,s) > 0 && canMergeByCompactness(r,s,compactness1,compactness2,compactness) )
      {
        diff = fabs( _regions[r].averageValue - _regions[s].averageValue );
        mergeRecord.region1 = r;
        mergeRecord.region2 = s;
        mergeRecord.score = diff;
        mergeRecordList.push_back( mergeRecord );
      }
    }
  }

  sort( mergeRecordList.begin(), mergeRecordList.end() );
  reverse( mergeRecordList.begin(), mergeRecordList.end() );

  /*!
   * Traitement de la liste des fusions
   */
  while ( mergeRecordList.empty() == false )
  {
    mergeRecord = mergeRecordList.back();
    r = mergeRecord.region1;
    s = mergeRecord.region2;
/*    cout << mergeRecord.score << " " << r << " " << s << " ";
    cout << "(" << _regions[r].size << "," << _regions[s].size << ")" << endl;*/
    labelMap[s+1] = labelMap[r+1];
    numMerges++;

    /*!
     * Suppression de la région s au profit de la région r
     */
    _regions[r].averageValue *= _regions[r].size;
    _regions[r].averageValue += _regions[s].size * _regions[s].averageValue;
    _regions[r].size += _regions[s].size;
    _regions[r].averageValue /= _regions[r].size;
    _regions[r].vertices->append( *_regions[s].vertices );
    _regions[r].surface += _regions[s].surface - 2*_adjacencyGraph(r,s);
    _regions[s].size = 0;
    _regions[s].surface = 0;
    for (t = 0; t < numRegions; t++)
    {
      if ( t != r && t != s )
      {
        _adjacencyGraph(r,t) += _adjacencyGraph(s,t);
        _adjacencyGraph(t,r) += _adjacencyGraph(t,s);
      }
      _adjacencyGraph(s,t) = 0;
      _adjacencyGraph(t,s) = 0;
    }

    mergeRecordList.pop_back();
    for (int i = mergeRecordList.size()-1; i >= 0; i--)
    {
      mergeRecord = mergeRecordList[i];
      if ( mergeRecord.region1 == r || mergeRecord.region1 == s
        || mergeRecord.region2 == r || mergeRecord.region2 == s )
      {
        mergeRecordList.erase( mergeRecordList.begin() + i );
      }
    }

    for (s = 0; s < numRegions; s++)
    {
      if ( _adjacencyGraph(r,s) > 0 && canMergeByCompactness(r,s,compactness1,compactness2,compactness) )
      {
        diff = fabs( _regions[r].averageValue - _regions[s].averageValue );
        mergeRecord.region1 = s > r? r: s;
        mergeRecord.region2 = s > r? s: r;
        mergeRecord.score = diff;
        mergeRecordList.push_back( mergeRecord );
      }
    }

    sort( mergeRecordList.begin(), mergeRecordList.end() );
    reverse( mergeRecordList.begin(), mergeRecordList.end() );
  }

  /*!
   * Résolution des équivalences et
   * condensation de la table des labels
   */
  int numEquivs = 0;
  for (r = 1; r <= numRegions; r++)
  {
    if ( labelMap[r] == r )
    {
      labelMap[r] -= numEquivs;
    }
    else
    {
      labelMap[r] = labelMap[labelMap[r]];
      numEquivs++;
    }
  }

  applyLabelMap( labelMap );

  EVAL( numMerges );
  LEAVE();

  return numMerges;
}

/*! Fusionne les régions voisines pour lesquelles la différence sur
 * \c regionFeature est inférieure ou égale au seuil \c featureThreshold.
 *
 * La fusion est réalisée en commençant par fusionner les régions qui
 * diffèrent le moins.
 * Le processus est répété jusqu'à ce qu'il ne reste plus de régions
 * qui puissent être fusionnées.
*/
template<class T>
int RegionAnalysis<T>::mergeRegions(
  RegionFeature regionFeature,
  const double featureThreshold)
{
  ENTER( "void RegionAnalysis::mergeRegions(...)" );

  const int numRegions = _regions.getSize();
  vector<MergeRecord> mergeRecordList;
  MergeRecord mergeRecord;
  int numMerges = 0;
  int r, s, t;
  T score;

  Vectori labelMap( numRegions+1 );
  for (r = 0; r <= numRegions; r++)
  {
    labelMap[r] = r;
  }

  /*!
   * Initialisation de la liste des fusions
   */
  for (r = 0; r < numRegions-1; r++)
  {
    for (s = r+1; s < numRegions; s++)
    {
      if ( _adjacencyGraph(r,s) > 0 )
      {
        score = fabs( _regions[r].feature(regionFeature) - _regions[s].feature(regionFeature) );
        if ( score <= featureThreshold )
        {
          mergeRecord.region1 = r;
          mergeRecord.region2 = s;
          mergeRecord.score = score;
          mergeRecordList.push_back( mergeRecord );
        }
      }
    }
  }

  sort( mergeRecordList.begin(), mergeRecordList.end() );
  reverse( mergeRecordList.begin(), mergeRecordList.end() );

  /*!
   * Traitement de la liste des fusions
   */
  while ( mergeRecordList.empty() == false )
  {
    mergeRecord = mergeRecordList.back();
    r = mergeRecord.region1;
    s = mergeRecord.region2;
    labelMap[s+1] = labelMap[r+1];
    numMerges++;

    /*!
     * Suppression de la région s au profit de la région r
     */
    _regions[r].averageValue *= _regions[r].size;
    _regions[r].averageValue += _regions[s].size * _regions[s].averageValue;
    _regions[r].size += _regions[s].size;
    _regions[r].averageValue /= _regions[r].size;
    _regions[r].vertices->append( *_regions[s].vertices );
    _regions[r].surface += _regions[s].surface - 2*_adjacencyGraph(r,s);
    _regions[s].clear();
    for (t = 0; t < numRegions; t++)
    {
      if ( t != r && t != s )
      {
        _adjacencyGraph(r,t) += _adjacencyGraph(s,t);
        _adjacencyGraph(t,r) += _adjacencyGraph(t,s);
      }
      _adjacencyGraph(s,t) = 0;
      _adjacencyGraph(t,s) = 0;
    }

    /*!
     * Suppression de la liste de tous les items concernant
     * les deux régions (r et s) qui viennent d'être fusionnées
     */
    mergeRecordList.pop_back();
    for (int i = mergeRecordList.size()-1; i >= 0; i--)
    {
      mergeRecord = mergeRecordList[i];
      if ( mergeRecord.region1 == r || mergeRecord.region1 == s
        || mergeRecord.region2 == r || mergeRecord.region2 == s )
      {
        mergeRecordList.erase( mergeRecordList.begin() + i );
      }
    }

    for (s = 0; s < numRegions; s++)
    {
      if ( _adjacencyGraph(r,s) > 0 )
      {
        score = fabs( _regions[r].feature(regionFeature) - _regions[s].feature(regionFeature) );
        if ( score < featureThreshold )
        {
          mergeRecord.region1 = s > r? r: s; //! region1: indice le plus petit
          mergeRecord.region2 = s > r? s: r; //! region2: indice le plus grand
          mergeRecord.score = score;
          mergeRecordList.push_back( mergeRecord );
        }
      }
    }

    sort( mergeRecordList.begin(), mergeRecordList.end() );
    reverse( mergeRecordList.begin(), mergeRecordList.end() );
  }

  /*!
   * Résolution des équivalences et
   * condensation de la table des labels
   */
  int numEquivs = 0;
  for (r = 1; r <= numRegions; r++)
  {
    if ( labelMap[r] == r )
    {
      labelMap[r] -= numEquivs;
    }
    else
    {
      labelMap[r] = labelMap[labelMap[r]];
      numEquivs++;
    }
  }

  applyLabelMap( labelMap );

  EVAL( numMerges );
  LEAVE();

  return numMerges;
}

template<class T>
int RegionAnalysis<T>::mergeNeighbouringRegions()
{
  ENTER( "void RegionAnalysis::mergeNeighbouringRegions()" );

  Vector<Region>& regions = _regions;
  const int numRegions = regions.getSize();
  int r, s;

  Vectori labelMap( numRegions+1 );
  for (r = 0; r <= numRegions; r++)
  {
    labelMap[r] = r;
  }

  for (r = 0; r < numRegions; r++)
  {
    for (s = r+1; s < numRegions; s++)
    {
      if ( _adjacencyGraph(r,s) > 0 )
      {
        if ( labelMap[s+1] > labelMap[r+1] )
        {
          labelMap[s+1] = r+1;
        }
        else
        {
          labelMap[r+1] = labelMap[s+1];
        }
      }
    }
  }

  int numMerges = 0;
  for (r = 1; r <= numRegions; r++)
  {
    if ( labelMap[r] == r )
    {
      labelMap[r] -= numMerges;
    }
    else
    {
      labelMap[r] = labelMap[labelMap[r]];
      numMerges++;
    }
  }
  EVAL( numMerges );

  if ( numMerges > 0 )
  {
    applyLabelMap( labelMap );
  }

#if 0
  for (r = 0; r < numRegions; r++)
  {
    for (s = r+1; s < numRegions; s++)
    {
      if ( _adjacencyGraph(r,s) > 0 )
      {
        for (t = 0; t < numRegions; t++)
        {
          if ( r != t && _adjacencyGraph(s,t) )
          {
            _adjacencyGraph(r,t) = 1;
            _adjacencyGraph(t,r) = 1;
          }
        }
      }
    }
  }

  labelMap[0] = 0;
  for (r = 0; r < numRegions; r++)
  {
    labelMap[r+1] = r+1;
    for (s = 0; s < r; s++)
    {
      if ( _adjacencyGraph(r,s) )
      {
        labelMap[r+1] = labelMap[s+1];
        numMerges++;
        break;
      }
    }
    if ( labelMap[r+1] == r+1 )
    {
      labelMap[r+1] -= numMerges;
    }
  }

  EVAL( numMerges );
  EVAL( labelMap );

  if ( numMerges > 0 )
  {
    VoxelMatrix<T>& labels = *_labelMatrix;
    const int size1 = labels.getSize1();
    const int size2 = labels.getSize2();
    const int size3 = labels.getSize3();
    int i, j, k, label;

    for (i = 0; i < size1; i++)
    {
      for (j = 0; j < size2; j++)
      {
        for (k = 0; k < size3; k++)
        {
          label = int( labels(i,j,k) );
          labels(i,j,k) = labelMap[label];
        }
      }
    }
  }
#endif

  LEAVE();

  return numMerges;
}

template<class T>
int RegionAnalysis<T>::compareRegions(
  const int r,
  const int s,
  RegionFeature regionFeature,
  const bool strict) const
{
//   ENTER( "int RegionAnalysis<T>::compareRegions(...)" );

  if ( strict )
  {
    if ( _regions[r].feature(regionFeature) < _regions[s].feature(regionFeature) )
    {
      return -1;
    }
    else if ( _regions[r].feature(regionFeature) > _regions[s].feature(regionFeature) )
    {
      return 1;
    }
    else
    {
      return 0;
    }
  }
  else
  {
    Vector<T> values1, values2;
    double mean1, mean2;
    double var1, var2, var;
    double t;
    int n1, n2;

    values1 = regionValues( _regions[r] );
    values2 = regionValues( _regions[s] );
    values1.apply( sqrt );
    values2.apply( sqrt );
    mean1 = values1.mean();
    mean2 = values2.mean();
    var1 = values1.var();
    var2 = values2.var();
    n1 = values1.getSize();
    n2 = values2.getSize();
    var = ( (n1-1) * var1 + (n2-1) * var2 ) / ( n1+n2-2 );
    t = ( mean1 - mean2 ) / sqrt( var/n1 + var/n2 );

//     LEAVE();

    if ( t > 1.96 )
    {
      return 1;
    }
    else if ( t < -1.96 )
    {
      return -1;
    }
    else
    {
      return 0;
    }
  }
}

template<class T>
void RegionAnalysis<T>::removeNonLocalMaximaRegions(RegionFeature regionFeature)
{
  ENTER( "void RegionAnalysis<T>::removeNonLocalMaximaRegions(...)" );

  const int numRegions = _regions.getSize();
  Vectori labelMap( numRegions+1 );
  int r, s, numRemoved = 0;

  labelMap[0] = 0;
  for (r = 0; r < numRegions; r++)
  {
    labelMap[r+1] = r+1;
    for (s = 0; s < numRegions; s++)
    {
      if ( _adjacencyGraph(r,s) > 0 )
      {
        if ( compareRegions(r,s,regionFeature,true) < 0 )
        {
          labelMap[r+1] = 0;
          break;
        }
      }
    }
  }

  for (r = 0; r < numRegions; r++)
  {
    if ( labelMap[r+1] == r+1 )
    {
      labelMap[r+1] -= numRemoved;
    }
    else
    {
      numRemoved++;
    }
  }

  if ( numRemoved > 0 )
  {
    applyLabelMap( labelMap );
  }

  EVAL( numRemoved );
  LEAVE();
}

template<class T>
void RegionAnalysis<T>::removeUnsalientRegions(
  RegionFeature regionFeature,
  const double saliency)
{
  ENTER( "void RegionAnalysis<T>::removeUnsalientRegions(...)" );

  const int numRegions = _regions.getSize();
  Vectori labelMap( numRegions+1 );
  int r, s, numRemoved = 0;
  double diff;
  bool isSalient;

  labelMap[0] = 0;
  for (r = 0; r < numRegions; r++)
  {
    isSalient = false;
    for (s = 0; s < numRegions; s++)
    {
      if ( _adjacencyGraph(r,s) > 0 )
      {
        diff = _regions[r].feature( regionFeature ) - _regions[s].feature( regionFeature );
        if ( diff >= saliency )
        {
          isSalient = true;
          break;
        }
      }
    }

    if ( isSalient )
    {
      labelMap[r+1] = r + 1 - numRemoved;
    }
    else
    {
      labelMap[r+1] = 0;
      numRemoved++;
    }
  }

  if ( numRemoved > 0 )
  {
    applyLabelMap( labelMap );
  }

  EVAL( numRemoved );
  LEAVE();
}

#if 0
/*! Supprime les régions dont la taille est strictement inférieure
 * à \c minimumSize.
 *
 * Retourne le nombre de régions supprimées.
****************************************************************/
template<class T>
int RegionAnalysis<T>::removeSmallRegions(const int minimumSize)
{
  ENTER( "int RegionAnalysis<T>::removeSmallRegions(const int)" );

  const int numRegions = _regions.getSize();
  Vectori labelMap( numRegions+1 );
  int r, numRetained = 0;

  labelMap[0] = 0;

  for (r = 0; r < numRegions; r++)
  {
    if ( _regions[r].size < minimumSize )
    {
      labelMap[r+1] = 0;
    }
    else
    {
      labelMap[r+1] = ++numRetained;
    }
  }

  if ( numRetained != numRegions )
  {
    applyLabelMap( labelMap );
  }

  int numRemoved = numRegions - numRetained;
  EVAL( numRemoved );

  LEAVE();

  return numRegions - numRetained;
}
#endif

template<class T>
void RegionAnalysis<T>::removeDarkRegions(const double minAverage)
{
  int r, s, v, i, j, k;
  Vectori flags( _regions.getSize() );

  flags.fill( 1 );

  for (r = 0; r < int(_regions.getSize()); r++)
  {
    if ( _regions[r].averageValue < minAverage )
    {
      const Verticesi& vertices = *_regions[r].vertices;
      for (v = 0; v < int(vertices.getSize()); v++)
      {
        i = vertices[v][X];
        j = vertices[v][Y];
        k = vertices[v][Z];
        (*_labelMatrix)(i,j,k) = 0;
      }

      flags[r] = 0;
    }
  }

  Vector<Region> remainingRegions( flags.sum() );
  int numRemoved = 0;
  for (r = s = 0; r < int(_regions.getSize()); r++)
  {
    if ( flags[r] )
    {
      remainingRegions[s++] = _regions[r];
      flags[r] = _regions[r].label - numRemoved;
    }
    else
    {
      numRemoved++;
    }
  }

  EVAL( numRemoved );
  EVAL( flags );

  _regions = remainingRegions;

  VoxelMatrix<T>& labels = *_labelMatrix;
  const int size1 = labels.getSize1();
  const int size2 = labels.getSize2();
  const int size3 = labels.getSize3();
  int label;

  for (i = 0; i < size1; i++)
  {
    for (j = 0; j < size2; j++)
    {
      for (k = 0; k < size3; k++)
      {
        label = int( labels(i,j,k) );
        if ( label > 0 )
        {
          labels(i,j,k) = flags[label-1];
        }
      }
    }
  }
}
#endif

/*! Returns the spatial calibration of the associated label matrix.
****************************************************************/
template<class T>
Vector<float> RegionAnalysis3D<T>::spatialCalibration() const
{
  return _labelMatrix->getVoxelCalibration().getVoxelSize();
}

template class RegionAnalysis3D<float>;
template class RegionAnalysis3D<double>;
template class RegionAnalysis3D<long double>;
