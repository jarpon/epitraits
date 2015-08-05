#ifndef LITE
/*!
 * \class  RegionAnalysis
 * \author Philippe Andrey (pa), INRA
 * \date    - création (pa)
 * \brief  Analyse de régions
****************************************************************/

#include "regionanalysis.h"
#include <componentlabelling.h>
#include <sumfilter.h>

#include <maths.h>
#include <matherror.h>
#include <boundingbox.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <programerror.h>

//#define TRACE
#include <trace.h>

int Region::__numRegions = 0;

Region::Region()
{
   vertices = 0;
   realVertices = 0;
   __numRegions++;
}

//Region::Region(const Region& region)
//{
//  *this = region;
//   __numRegions++;
//}

Region::~Region()
{
   delete vertices;
   delete realVertices;
  __numRegions--;
//    EVAL( __numRegions );
}

void Region::clear()
{
  label = 0;
  //size = 0;
  surface = 0;
  delete vertices;
  vertices = 0;
  delete realVertices;
  realVertices = 0;
  edgeFlags.setSize( 0 );
}

unsigned int Region::getSize() const
{
  if ( vertices )
  {
    return vertices->getSize();
  }
  else
  {
    return 0;
  }
}

//Region& Region::operator=(const Region& region)
//{
//  if ( this != &region )
//  {
//    label = region.label;
//    size = region.size;
//    surface = region.surface;
//    delete vertices;
//    vertices = new Vertices( *region.vertices ); // detach???
//    edgeFlags = region.edgeFlags;
//  }
//
//  return *this;
//}


/*!*******************************************************
**********************************************************
**********************************************************/

template<class T>
RegionAnalysis<T>::RegionAnalysis()
{
  _regionMatrix = 0;
//  _voxelMatrix = 0;
}

template<class T>
RegionAnalysis<T>::~RegionAnalysis()
{
}

template<class T>
void RegionAnalysis<T>::setRegionMatrix(VoxelMatrix<T>& regionMatrix)
{
  _regionMatrix = &regionMatrix;
}

template<class T>
const VoxelMatrix<T>& RegionAnalysis<T>::getRegionMatrix() const
{
  return *_regionMatrix;
}

template<class T>
VoxelMatrix<T>& RegionAnalysis<T>::getRegionMatrix()
{
  return *_regionMatrix;
}

template<class T>
bool RegionAnalysis<T>::hasRegionMatrix() const
{
  return _regionMatrix != 0;
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
void RegionAnalysis<T>::applyLabelMap(const Vectori& labelMap)
{
  VoxelMatrix<T>& regionMatrix = *_regionMatrix;
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
void RegionAnalysis<T>::applyMajorityFilter()
{
  //ENTER( "void RegionAnalysis<T>::applyMajorityFilter()" );

  VoxelMatrix<T>& regionMatrix = *_regionMatrix;
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

  //LEAVE();
}

/*! Recalcule les valeurs de l'image des labels afin qu'il n'y
 * ait pas de trous dans la numérotation càd que tous les labels
 * entre 1 et la plus grande des valeurs de l'image des labels
 * soient représentés.
 \todo à mettre dans une classe LabelFilter
****************************************************************/
template<class T>
int RegionAnalysis<T>::condenseRegionLabels()
{
  VoxelMatrix<T>& regionMatrix = *_regionMatrix;
  const int maxLabel = int( regionMatrix.max().max().max() );
  const Vectorui histogram = regionMatrix.histogram( 0, 1, maxLabel+1 );
  Vectori labelMap( maxLabel+1 );
  int i, numLabels = 0;

  labelMap.setZeros();
  for (i = 1; i <= maxLabel; i++)
  {
    if ( histogram[i] > 0 )
    {
      labelMap[i] = ++numLabels;
    }
  }

  applyLabelMap( labelMap );

  return numLabels;
}


/* Run!!!
***********************************************/
template<class T>
void RegionAnalysis<T>::run()
{
  //ENTER( "void RegionAnalysis<T>::run()" );

  const VoxelMatrix<T>& regionMatrix = getRegionMatrix();
  Vector<Region>& regions = _regions;
  const int size1 = regionMatrix.getSize1();
  const int size2 = regionMatrix.getSize2();
  const int size3 = regionMatrix.getSize3();
  const int numRegions = int( regionMatrix.max().max().max() );
  int i, j, k, u, v, w, r, s, n, currentVertex;
  Vectori currentVertices;

  regions.setSize( numRegions );
  _adjacencyGraph.setZeros( numRegions );
  currentVertices.setZeros( numRegions );
  Vectorui histogram = regionMatrix.histogram( 0.0, 1.0, numRegions+1 );
  //EVAL(histogram);

  for (r = 0; r < numRegions; r++)
  {
    regions[r].label = r + 1;
    regions[r].surface = 0;
    regions[r].vertices = new Verticesi( 3, histogram[r+1] );
    regions[r].realVertices = new Vertices<double>( 3, histogram[r+1] );
    regions[r].edgeFlags.setSize( histogram[r+1] );
  }

  Matrixi shifts( 6, 3 );
  shifts.setZeros();
  shifts(0,X) = -1;
  shifts(1,X) = 1;
  shifts(2,Y) = -1;
  shifts(3,Y) = 1;
  shifts(4,Z) = -1;
  shifts(5,Z) = 1;

  const VoxelCalibration& calibration = regionMatrix.getVoxelCalibration();
  const float xVoxelSize = calibration.getVoxelSize()[0];
  const float yVoxelSize = calibration.getVoxelSize()[1];
  const float zVoxelSize = calibration.getVoxelSize()[2];

  for (i = 0; i < size1; i++)
  {
    for (j = 0; j < size2; j++)
    {
      for (k = 0; k < size3; k++)
      {
        r = int( regionMatrix(i,j,k) );
        if ( r > 0 )
        {
          Region& region = regions[r-1];
          currentVertex = currentVertices[r-1]++;
          //units in pixels
          (*region.vertices)[currentVertex][X] = i;
          (*region.vertices)[currentVertex][Y] = j;
          (*region.vertices)[currentVertex][Z] = k;
          //real size
          (*region.realVertices)[currentVertex][X] = i*xVoxelSize;
          (*region.realVertices)[currentVertex][Y] = j*yVoxelSize;
          (*region.realVertices)[currentVertex][Z] = k*zVoxelSize;
          region.edgeFlags[currentVertex] = false;

          for (n = 0; n < 6; n++)
          {
            u = i + shifts(n,X);
            v = j + shifts(n,Y);
            w = k + shifts(n,Z);

            if ( u >= 0 && u < size1
              && v >= 0 && v < size2
              && w >= 0 && w < size3 )
            {
              s = static_cast<int>( regionMatrix(u,v,w) );
              if ( s != r )
              {
                region.edgeFlags[currentVertex] = true;
                if ( s > 0 )
                {
                  _adjacencyGraph(r-1,s-1)++;
                }
                if ( u != i ) region.surface += yVoxelSize * zVoxelSize;
                if ( v != j ) region.surface += xVoxelSize * zVoxelSize;
                if ( w != k ) region.surface += xVoxelSize * yVoxelSize;
              }
            }
          }
        }
      }
    }
  }


  for (r = 0; r < numRegions; ++r)
  {
//    EVAL( _regions[r].vertices->getSize() );
//    EVAL( _regions[r].realVertices->getSize() );
//    EVAL( _regions[r].realVertices->cog() );
    _regions[r].eigenVectors = _regions[r].realVertices->principalAxes( _regions[r].eigenValues );
  }

  if ( _adjacencyGraph.isSymmetrical() == false )
  {
    ProgramError programError;
    programError.setWhere( "void RegionAnalysis<T>::run()" );
    programError.setWhat( "Asymmetrical region adjacency graph" );
    throw programError;
  }

  //EVAL( numRegions );

  //LEAVE();
}


template<class T>
int RegionAnalysis<T>::numRegions() const
{
  return _regions.getSize();
}

template<class T>
Vector<Region>& RegionAnalysis<T>::getRegions()
{
  return _regions;
}

template<class T>
const Vector<Region>& RegionAnalysis<T>::getRegions() const
{
  return _regions;
}


/*! Affecte la valeur \c value à tous les voxels de \c voxelMatrix
 * qui appartiennent à la région \c region.
****************************************************************/
template<class T>
void RegionAnalysis<T>::fillRegion(
  VoxelMatrix<T>& voxelMatrix,
  const Region& region,
  const T value) const
{
  if ( region.getSize() > 0 )
  {
    const Verticesi& vertices = *region.vertices;
    int i, j, k;

    for (unsigned v = 0; v < vertices.getSize(); v++)
    {
      i = vertices[v][X];
      j = vertices[v][Y];
      k = vertices[v][Z];
      voxelMatrix(i,j,k) = value;
    }
  }
}

/*! Génère une cartographie des valeurs \c values.
 *
 * Chaque région est remplie avec la valeur qui lui correspond
 * dans le tableau \c values.
 * Le résultat est stocké dans \c voxelMatrix, qui est redimensionnée
 * si nécessaire.
 * La valeur 0 est attribuée aux voxels qui n'appartiennent à
 * aucune région.
 *
 * \exception SizeError le nombre de valeurs passé diffère du
 * nombre de régions.
****************************************************************/
template<class T>
void RegionAnalysis<T>::fillRegions(
  VoxelMatrix<T>& voxelMatrix,
  const Vector<T>& values) const
{
  //ENTER( "void RegionAnalysis<T>::fillRegions(...)" );

  if ( _regions.getSize() != values.getSize() )
  {
    SizeError sizeError;
    sizeError.setWhere( "void RegionAnalysis<T>::fillRegions(...)" );
    sizeError.setWhat( "The number of values differs from the number of regions" );
    throw sizeError;
  }

  voxelMatrix.setSize( _regionMatrix->getSize() );
  voxelMatrix.setZeros();

  for (unsigned int r = 0; r < _regions.getSize(); r++)
  {
    fillRegion( voxelMatrix, _regions[r], values[r] );
  }

  //LEAVE();
}

/*! Retourne le vecteur des valeurs de l'ensemble des régions.
 *
 * Les valeurs sont prises dans \c voxelMatrix.
****************************************************************/
template<class T>
Vector<T> RegionAnalysis<T>::allRegionValues(const VoxelMatrix<T>& voxelMatrix) const
{
  Vector<T> allValues;

  for (int r = 0; r < numRegions(); r++)
  {
    allValues.append( regionValues(getRegions()[r],voxelMatrix) );
  }

  return allValues;
}

/*! Retourne les valeurs des voxels de \c voxelMatrix qui
 * appartiennent à la région \c region.
****************************************************************/
template<class T>
Vector<T> RegionAnalysis<T>::regionValues(
  const Region& region,
  const VoxelMatrix<T>& voxelMatrix) const
{
  const Verticesi& vertices = *region.vertices;
  Vector<T> values( vertices.getSize() );
  int i, j, k;

  for (unsigned int v = 0; v < vertices.getSize(); v++)
  {
    i = vertices[v][X];
    j = vertices[v][Y];
    k = vertices[v][Z];
    values[v] = voxelMatrix(i,j,k);
  }

  return values;
}

template<class T>
Vector<T> RegionAnalysis<T>::randomOutsideValues(
  const Region& region,
  const VoxelMatrix<T>& voxelMatrix,
  RandomGenerator& randomGenerator) const
{
  const int numRegions = _regions.getSize();
  const int numValues = region.getSize();
  const int totalNumValues = computeRegionFeature( REGION_FEATURE_SIZE, voxelMatrix ).sum() - numValues;
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
        const Vectori& vertex = (*_regions[r].vertices)[p];
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
Vector<T> RegionAnalysis<T>::regionCentroid(
  const Region& region,
  const VoxelMatrix<T>& voxelMatrix) const
{
  const Verticesi& vertices = *region.vertices;
  T svx = 0;
  T svy = 0;
  T svz = 0;
  T value;
  int i, j, k;

  for (unsigned int v = 0; v < vertices.getSize(); v++)
  {
    i = vertices[v][X];
    j = vertices[v][Y];
    k = vertices[v][Z];

    if ( voxelMatrix(i,j,k) != 0 )
    {
      svx += i;
      svy += j;
      svz += k;
      value += 1;
    }
  }

  Vector<T> centroid( 3 );
  centroid[X] = (svx / value) * voxelMatrix.getVoxelCalibration().getVoxelSize()[X];
  centroid[Y] = (svy / value) * voxelMatrix.getVoxelCalibration().getVoxelSize()[Y];
  centroid[Z] = (svz / value) * voxelMatrix.getVoxelCalibration().getVoxelSize()[Z];

  return centroid;
}

/*! Retourne le barycentre, coordonnées exprimées en voxel */
template<class T>
Vector<T> RegionAnalysis<T>::regionBarycenter(
  const Region& region,
  const VoxelMatrix<T>& voxelMatrix) const
{
  const Verticesi& vertices = *region.vertices;
  T svx = 0;
  T svy = 0;
  T svz = 0;
  T sv = 0;
  T value;
  int i, j, k;

  for (unsigned int v = 0; v < vertices.getSize(); v++)
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
  barycenter[X] = (svx / sv) * voxelMatrix.getVoxelCalibration().getVoxelSize()[X];
  barycenter[Y] = (svy / sv) * voxelMatrix.getVoxelCalibration().getVoxelSize()[Y];
  barycenter[Z] = (svz / sv) * voxelMatrix.getVoxelCalibration().getVoxelSize()[Z];

  return barycenter;
}

template<class T>
T RegionAnalysis<T>::regionKurtosis(const Region& region,const VoxelMatrix<T>& voxelMatrix) const
{
  const Vector<T> barycenter = regionBarycenter( region, voxelMatrix );
  const int dimension = barycenter.getSize();
  SquareMatrix<T> covarianceMatrix( dimension );
  const Verticesi& vertices = *region.vertices;
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


/*!intensity
****************************************************************/
template<class T>
T RegionAnalysis<T>::regionIntensity(const Region& region, const VoxelMatrix<T>& voxelMatrix) const
{
  //const Vector<T> intensities = regionIntensity( region, voxelMatrix );
  const Verticesi& vertices = *region.vertices;
  const int numVertices = vertices.getSize();
  T value, intensity = 0.0;
  Vector<T> buffer( 3 );
  int i, j, k, v;

  for (v = 0; v < numVertices; v++)
  {
    buffer[X] = i = vertices[v][X];
    buffer[Y] = j = vertices[v][Y];
    buffer[Z] = k = vertices[v][Z];
    value = voxelMatrix(i,j,k);
    intensity += value;
  }
  return intensity;
}


/*!integrated density
****************************************************************/
template<class T>
T RegionAnalysis<T>::regionIntegratedDensity(const Region& region, const VoxelMatrix<T>& voxelMatrix) const
{
  Vector <float> voxelSize = voxelMatrix.getVoxelCalibration().getVoxelSize();
  //EVAL( regionMatrix.max().max().max() );
  //EVAL (regionIntensity( region, voxelMatrix ) );
  //EVAL ( regionVolume( region ) );
  //EVAL ( voxelSize.prod() );
  return ( regionIntensity( region, voxelMatrix ) * voxelSize.prod() );
}


/*!volume
****************************************************************/
template<class T>
T RegionAnalysis<T>::regionVolume(const Region& region, const VoxelMatrix<T>& voxelMatrix) const
{
  Vector <float> voxelSize = voxelMatrix.getVoxelCalibration().getVoxelSize();
  return region.getSize()*voxelSize.prod();
}

/*!volume
****************************************************************/
template<class T>
T RegionAnalysis<T>::equivalentRadius(const Region& region, const VoxelMatrix<T>& voxelMatrix) const
{
  Vector <float> voxelSize = voxelMatrix.getVoxelCalibration().getVoxelSize();
  T volumeRegion = voxelSize.prod()*region.getSize();
  //return pow(3*volumeRegion/(4*M_PI),(1./3));
  return pow( (3.0*fabs(volumeRegion))/(4.0*M_PI), 1.0/3.0 );
}

/*!trying to get the flatness!
****************************************************************/
template<class T>
T RegionAnalysis<T>::regionFlatness(const Region& region) const
{
  //EVAL(region.realVertices->cog());
  //EVAL(region.realVertices->inertia());
  //EVAL(region.eigenValues);
  //EVAL(region.eigenVectors);
  return sqrt( region.eigenValues[1] / region.eigenValues[2] );
}

/*!trying to get the elongation!
****************************************************************/
template<class T>
T RegionAnalysis<T>::regionElongation(const Region& region) const
{
  return sqrt( region.eigenValues[0] / region.eigenValues[1] );
}

/*!trying to get the sphericity!
****************************************************************/
template<class T>
T RegionAnalysis<T>::regionSphericity(const Region& region, const VoxelMatrix<T>& voxelMatrix) const
{
  //VoxelMatrix<T>& regionMatrix = *_regionMatrix;
  Vector <float> voxelSize = voxelMatrix.getVoxelCalibration().getVoxelSize();
  T volumeRegion = voxelSize.prod()*region.getSize();
  //T surfaceArea = 4 * M_PI * pow( ( region.eigenValues[0] + region.eigenValues[1] + region.eigenValues[2]) / 6, 2);
  return ( ( 36 * M_PI * pow(volumeRegion, 2) ) / pow( region.surface, 3) );

  //return ( ( pow(M_PI,(2./3) )*pow( (6*volumeRegion),(2./3) ) ) / surfaceArea);
}

/*!trying to get the surface area!
****************************************************************
template<class T>
T RegionAnalysis<T>::regionSurfaceArea() const
{
  VoxelMatrix<T>& regionMatrix = *_regionMatrix;
  Vector <float> voxelSize = regionMatrix.getVoxelCalibration().getVoxelSize();
  return = voxelSize.prod()*regionMatrix.sum().sum().sum();
}
*/

template<class T>
Vector<T> RegionAnalysis<T>::computeBarycenter(
  const Region& region,
  RegionFeature regionFeature,
  const VoxelMatrix<T>& voxelMatrix) const
{
  Vector<T> featureValue(3);

  switch ( regionFeature )
  {
  case REGION_FEATURE_BARYCENTER_VALUE:
    featureValue = regionBarycenter( region, voxelMatrix );
  break;
  }
  return featureValue;
}

template<class T>
Vector<T> RegionAnalysis<T>::computeCentroid(
  const Region& region,
  RegionFeature regionFeature,
  const VoxelMatrix<T>& voxelMatrix) const
{
  Vector<T> featureValue(3);

  switch ( regionFeature )
  {
  case REGION_FEATURE_CENTROID:
    featureValue = regionCentroid( region, voxelMatrix );
  break;
  }
  return featureValue;
}

template<class T>
T RegionAnalysis<T>::computeRegionFeature(
  const Region& region,
  RegionFeature regionFeature,
  const VoxelMatrix<T>& voxelMatrix) const
{
  //ENTER( "T RegionAnalysis<T>::computeRegionFeature(...,...,...) const" );

  T featureValue;

  switch ( regionFeature )
  {
//    case REGION_FEATURE_COMPACTNESS:
//      featureValue = compactness( region );
//      break;

    case REGION_FEATURE_VOLUME:
      featureValue = regionVolume( region, voxelMatrix );
    break;

    case REGION_FEATURE_EQUIVALENT_RADIUS:
      featureValue = equivalentRadius( region, voxelMatrix );
    break;

    case REGION_FEATURE_INTENSITY:
      featureValue = regionIntensity( region, voxelMatrix );
    break;

    case REGION_FEATURE_INTEGRATED_DENSITY:
      featureValue = regionIntegratedDensity( region, voxelMatrix );
    break;

    case REGION_FEATURE_FLATNESS:
      featureValue = regionFlatness( region );
    break;

    case REGION_FEATURE_ELONGATION:
      featureValue = regionElongation( region );
    break;

    case REGION_FEATURE_SURFACE_AREA:
      featureValue = region.surface;
    break;

    case REGION_FEATURE_SPHERICITY:
      featureValue = regionSphericity( region, voxelMatrix );
    break;

    case REGION_FEATURE_SIZE:
      featureValue = region.getSize();
    break;

    case REGION_FEATURE_MEDIAN_VALUE:
      featureValue = regionValues(region,voxelMatrix).median();
    break;

    case REGION_FEATURE_AVERAGE_VALUE:
      featureValue = regionValues(region,voxelMatrix).mean();
    break;

    case REGION_FEATURE_MINIMUM_VALUE:
      featureValue = regionValues(region,voxelMatrix).min();
    break;

    case REGION_FEATURE_MAXIMUM_VALUE:
      featureValue = regionValues(region,voxelMatrix).max();
    break;

    case REGION_FEATURE_MAD_VALUE:
      featureValue = regionValues(region,voxelMatrix).mad();
    break;

    case REGION_FEATURE_KURTOSIS:
      featureValue = regionKurtosis( region, voxelMatrix );
    break;



//    case REGION_FEATURE_CORE_VALUE:
//      featureValue = coreValue( region );
//      break;

    default:
      ProgramError programError;
      throw programError;
  }

 // LEAVE();

  return featureValue;
}

/*! Retourne un vecteur donnant la valeur de l'attribut
 * \c regionFeature pour chacune des régions.
 *
 * La valeur est fixée à 0 pour les régions vides.
****************************************************************/
template<class T>
Vertices<T> RegionAnalysis<T>::computeRegionBarycenters(
  RegionFeature regionFeature,
  const VoxelMatrix<T>& voxelMatrix)
{
  //ENTER( "Vector<T> RegionAnalysis<T>::computeRegionFeature(...,...) const" );

  const int numRegions = _regions.getSize();
  Vector<T> currentBarycenter( 3 );
  Vertices<T> valuesB( 3, 0, 0, 0 );


  for (int r = 0; r < numRegions; r++)
  {
    if ( _regions[r].getSize() > 0 )
    {
      if ( regionFeature == REGION_FEATURE_BARYCENTER_VALUE )
      {
        currentBarycenter = computeBarycenter( _regions[r], regionFeature, voxelMatrix );
        valuesB.append(currentBarycenter);
      }
    }
    else
    {
      //values[r] = 0.0;
      return valuesB;
    }
  }

  //LEAVE();

  return valuesB;
}

/*! Retourne un vecteur donnant la valeur de l'attribut
 * \c regionFeature pour chacune des régions.
 *
 * La valeur est fixée à 0 pour les régions vides.
****************************************************************/
template<class T>
Vertices<T> RegionAnalysis<T>::computeRegionCentroids(
  RegionFeature regionFeature,
  const VoxelMatrix<T>& voxelMatrix)
{
  //ENTER( "Vector<T> RegionAnalysis<T>::computeRegionFeature(...,...) const" );

  const int numRegions = _regions.getSize();
  Vector<T> currentCentroid( 3 );
  Vertices<T> allCentroids( 3, 0, 0, 0 );


  for (int r = 0; r < numRegions; r++)
  {
    if ( _regions[r].getSize() > 0 )
    {
      if ( regionFeature == REGION_FEATURE_CENTROID )
      {
        currentCentroid = computeCentroid( _regions[r], regionFeature, voxelMatrix );
        allCentroids.append(currentCentroid);
      }
    }
    else
    {
      //values[r] = 0.0;
      return allCentroids;
    }
  }

  //LEAVE();

  return allCentroids;
}

/*! Retourne un vecteur donnant la valeur de l'attribut
 * \c regionFeature pour chacune des régions.
 *
 * La valeur est fixée à 0 pour les régions vides.
****************************************************************/
template<class T>
Vector<T> RegionAnalysis<T>::computeRegionFeature(
  RegionFeature regionFeature,
  const VoxelMatrix<T>& voxelMatrix) const
{
  //ENTER( "Vector<T> RegionAnalysis<T>::computeRegionFeature(...,...) const" );

  const int numRegions = _regions.getSize();
  Vector<T> values( numRegions );

  for (int r = 0; r < numRegions; r++)
  {
    if ( _regions[r].getSize() > 0 )
    {
      values[r] = computeRegionFeature( _regions[r], regionFeature, voxelMatrix );
    }
    else
    {
      values[r] = 0.0;
    }
  }

  //LEAVE();

  return values;
}

template<class T>
void RegionAnalysis<T>::mapRegionFeature(
  VoxelMatrix<T>& mapMatrix,
  RegionFeature regionFeature,
  const VoxelMatrix<T>& voxelMatrix) const
{
  fillRegions( mapMatrix, computeRegionFeature(regionFeature,voxelMatrix) );
}

template<class T>
void RegionAnalysis<T>::thresholdRegions(
  const Vector<T>& regionValues,
  const T threshold)
{
  if ( _regions.getSize() != regionValues.getSize() )
  {
    SizeError sizeError;
    sizeError.setWhere( "void RegionAnalysis<T>::thresholdRegions(...)" );
    sizeError.setWhat( "The number of values differs from the number of regions" );
    throw sizeError;
  }

  for (unsigned int r = 0; r < _regions.getSize(); r++)
  {
    if ( regionValues[r] < threshold )
    {
      fillRegion( getRegionMatrix(), _regions[r], 0 );
      _regions[r].clear();
    }
  }
}

template class RegionAnalysis<float>;
template class RegionAnalysis<double>;

#endif // LITE
