#ifndef LITE
#ifndef REGIONANALYSIS2_H
#define REGIONANALYSIS2_H

#include <voxelmatrix.h>
#include <alinesquarematrix.h>
#include <shape.h>
#include <vertices.h>

typedef enum
{
  REGION_FEATURE_SIZE,
  REGION_FEATURE_VOLUME,
  REGION_FEATURE_EQUIVALENT_RADIUS,
  REGION_FEATURE_CENTROID,
  REGION_FEATURE_INTENSITY,
  REGION_FEATURE_INTEGRATED_DENSITY,
  REGION_FEATURE_COMPACTNESS,
  REGION_FEATURE_FLATNESS,
  REGION_FEATURE_ELONGATION,
  REGION_FEATURE_SURFACE_AREA,
  REGION_FEATURE_SPHERICITY,
  REGION_FEATURE_AVERAGE_VALUE,
  REGION_FEATURE_MINIMUM_VALUE,
  REGION_FEATURE_MAXIMUM_VALUE,
  REGION_FEATURE_MEDIAN_VALUE,
  REGION_FEATURE_MAD_VALUE,
  REGION_FEATURE_KURTOSIS,
  REGION_FEATURE_BARYCENTER_VALUE,
  REGION_FEATURE_CORE_VALUE,
  REGION_FEATURE_CONTRAST,
  REGION_FEATURE_CONTRACTNESS,
  REGION_FEATURE_RELATIVE_CONTRAST
} RegionFeature;

class Region
{
  public:

    Region();
    Region(const Region&);
    Region& operator=(const Region&);
    ~Region();

    void clear();
    unsigned int getSize() const;

    void computeFeature(RegionFeature);
    double feature(RegionFeature) const;
    //SquareMatrix<float> feature(RegionFeature) const;

    int label;
   // int size;
    float surface;
    SquareMatrix<double> eigenVectors;
    Vector<double> eigenValues;

    Vertices<int>* vertices;
    Vertices<double>* realVertices;
    Vector<bool> edgeFlags;

    static int __numRegions;
};

inline ostream& operator<<(ostream& os, const Region& region)
{
  os << "Region" << endl;
  os << "  label = " << region.label << endl;
  //os << "  size  = " << region.size << endl;
  os << "  surface  = " << region.surface << endl;
  os << "  vertices = " << region.vertices << endl;
  os << "  realVertices = " << region.realVertices << endl;
  return os;
}

template<class T>
class RegionAnalysis
{
  public:

    RegionAnalysis();
    ~RegionAnalysis();

    void setRegionMatrix(VoxelMatrix<T>&);
    VoxelMatrix<T>& getRegionMatrix();
    const VoxelMatrix<T>& getRegionMatrix() const;
    bool hasRegionMatrix() const;

//    void setVoxelMatrix(const VoxelMatrix<T>&);
//    const VoxelMatrix<T>& getVoxelMatrix() const;
//    bool hasVoxelMatrix() const;

    void applyLabelMap(const Vectori&); // private?
    void applyMajorityFilter();
    int condenseRegionLabels();

    void run();
    int numRegions() const;
    Vector<Region>& getRegions();
    const Vector<Region>& getRegions() const;


    void fillRegion(VoxelMatrix<T>&,const Region&,const T) const;
    void fillRegions(VoxelMatrix<T>&,const Vector<T>&) const;

    Vector<T> allRegionValues(const VoxelMatrix<T>&) const;
    Vector<T> regionValues(const Region&,const VoxelMatrix<T>&) const;
    Vector<T> randomOutsideValues(const Region&, const VoxelMatrix<T>&, RandomGenerator&) const;
    Vector<T> regionBarycenter(const Region&,const VoxelMatrix<T>&) const;
    Vector<T> regionCentroid(const Region&,const VoxelMatrix<T>&) const;
    T regionKurtosis(const Region&,const VoxelMatrix<T>&) const;
    T regionVolume(const Region&, const VoxelMatrix<T>&) const;
    T regionIntensity(const Region&, const VoxelMatrix<T>&) const;
    T regionIntegratedDensity(const Region&, const VoxelMatrix<T>&) const;
    T regionFlatness(const Region&) const;
    T regionElongation(const Region&) const;
//    T regionSurfaceArea(const Region&) const;
    T regionSphericity(const Region&, const VoxelMatrix<T>&) const;
    T equivalentRadius(const Region&, const VoxelMatrix<T>&) const;
    T computeRegionFeature(const Region&, RegionFeature, const VoxelMatrix<T>&) const;
    Vector<T> computeBarycenter(const Region&, RegionFeature, const VoxelMatrix<T>&) const;
    Vector<T> computeCentroid(const Region&, RegionFeature, const VoxelMatrix<T>&) const;
    Vector<T> computeRegionFeature(RegionFeature,const VoxelMatrix<T>&) const;
    Vertices<T> computeRegionBarycenters(RegionFeature,const VoxelMatrix<T>&);
    Vertices<T> computeRegionCentroids(RegionFeature,const VoxelMatrix<T>&);
    void mapRegionFeature(VoxelMatrix<T>&,RegionFeature,const VoxelMatrix<T>&) const;

    void thresholdRegions(const Vector<T>&,const T);

#if 0
    double compactness(const Region&) const;
    //T barycenterValue(const Region&) const;
    T distanceWeightedValue(const Region&) const;
    T coreValue(const Region&) const;

    bool canMerge(const int,const int) const;
    void mergeRegions(const unsigned int,const unsigned int);
    void mergeRegions2(RegionFeature,const double);
    void mergeRegions(RegionFeature);
    void reconstructionClosing(Vectord&,const unsigned int);
    Vectord computeContrasts(const Vectord&,const bool) const;


    void remapRegionLabels();

    Vector<T> drawRandomValues(const Region&, const VoxelMatrix<T>&, Random&) const;

    void run();
    const Region& getBackground() const;
    const Matrixi& getAdjacencyGraph() const;

    Matrixd regionMatrix() const;
    void saveRegionMatrix(const string&) const;
    void saveRegionMatrix(const Matrixd&, const string&) const;
    void appendRegionMatrix(ofstream&) const;
    void appendRegionMatrix(const Matrixd&, ofstream&) const;
    string regionMatrixHeader() const;

    void overlayContours(VoxelMatrix<T>&) const;

    Vectori sizes() const;
    Vectord feature(RegionFeature) const;
//     void updateFeature(RegionFeature);
//     void mapFeature(VoxelMatrix<T>&,RegionFeature) const;

    float averageNumNeighbours() const;

    Vector<T> regionValues(const Region&) const;
    Vectorui computeRegionHistogram(const Region&, const T) const;

    void tophat();
//     void reconstructionClosing(const unsigned int);

    void enhanceContrast(const deeeeeeeouble);
    void bilateralFiltering(const double);
    void boxBilateralFiltering(const double);
    void computeRelativeContrasts(RegionFeature, const bool);
    int removeNegativeContrastRegions();
    int thresholdRegions(RegionFeature,const double);

    bool canMergeByCompactness(const int,const int) const;
    bool canMergeByCompactness(
      const int,const int,
      double&, double&,
      double&) const;
    int mergeRegions();
    int mergeRegionsByPriority();
    int mergeRegionsByPriority2();
    int mergeRegions(RegionFeature, const double);
    int mergeNeighbouringRegions();
    void removeNonLocalMaximaRegions(RegionFeature);
    void removeUnsalientRegions(RegionFeature, const double);
//     int removeSmallRegions(const int);
    void removeDarkRegions(const double);

    int compareRegions(const int, const int, RegionFeature, const bool) const;
#endif

  private:

    RegionAnalysis(const RegionAnalysis&);
    RegionAnalysis& operator=(const RegionAnalysis&);

    VoxelMatrix<T>* _regionMatrix;
//    const VoxelMatrix<T>* _voxelMatrix;
    Vector<Region> _regions;
    //Region* _region;
    SquareMatrixi _adjacencyGraph;
};

#endif
#endif // LITE
