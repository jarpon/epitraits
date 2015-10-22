#ifndef LITE
#ifndef REGIONANALYSIS3D_H
#define REGIONANALYSIS3D_H

//#include "regionanalysisbase.h"
//#include "voxelmatrix.h"
#include <regionanalysisbase.h>
#include <voxelmatrix.h>

template<class T>
class RegionAnalysis3D : public RegionAnalysisBase<T>
{
  public:

    RegionAnalysis3D();
    ~RegionAnalysis3D();

    int dimension() const { return 3; }

    bool hasLabelMatrix() const;
    void setLabelMatrix(VoxelMatrix<T>&);
    VoxelMatrix<T>& getLabelMatrix();
    const VoxelMatrix<T>& getLabelMatrix() const;

    bool hasValueMatrix() const;
    void setValueMatrix(const VoxelMatrix<T>&);
    const VoxelMatrix<T>& getValueMatrix() const;

    bool hasOutputMatrix() const;
    void setOutputMatrix(VoxelMatrix<T>&);
    VoxelMatrix<T>& getOutputMatrix() const;

    void applyLabelMap(const Vector<int>&); // private?
    void applyMajorityFilter();
    int condenseRegionLabels();

    void run();

    void outputFillRegion(const Region&,const T);

    Vector<T> allRegionValues(const VoxelMatrix<T>&) const;
    Vector<T> regionValues(const Region&) const;
    Vector<T> randomOutsideValues(const Region&, const VoxelMatrix<T>&, RandomGenerator&) const;
    Vector<T> regionBarycenter(const Region&,const VoxelMatrix<T>&) const;
    T regionKurtosis(const Region&,const VoxelMatrix<T>&) const;
    void mapRegionFeature(VoxelMatrix<T>&,RegionFeature) const;

    PixelMatrix<T> getLabel2DProjection(const T) const;

    void thresholdRegions(const Vector<T>&,const T);

    using RegionAnalysisBase<T>::numRegions;
    using RegionAnalysisBase<T>::getRegion;
    using RegionAnalysisBase<T>::getBoundaryValue;
    using RegionAnalysisBase<T>::isRegionLabel;

#if 0

    T barycenterValue(const Region&) const;
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

    const Region& getBackground() const;
    const Matrixi& getAdjacencyGraph() const;

    Matrixd regionMatrix() const;
    void saveRegionMatrix(const string&) const;
    void saveRegionMatrix(const Matrixd&, const string&) const;
    void appendRegionMatrix(ofstream&) const;
    void appendRegionMatrix(const Matrixd&, ofstream&) const;
    string regionMatrixHeader() const;

    void overlayContours(VoxelMatrix<T>&) const;

//     void mapFeature(VoxelMatrix<T>&,RegionFeature) const;

    float averageNumNeighbours() const;

    Vector<T> regionValues(const Region&) const;
    Vectorui computeRegionHistogram(const Region&, const T) const;

    void tophat();
//     void reconstructionClosing(const unsigned int);

    void enhanceContrast(const double);
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

    RegionAnalysis3D(const RegionAnalysis3D&);
    RegionAnalysis3D& operator=(const RegionAnalysis3D&);

    void computeRAGImplicitBoundaryMode();
    void computeRAGExplicitBoundaryMode();

    void fillRegion(const Region&,const T,VoxelMatrix<T>&);

    virtual Vector<float> spatialCalibration() const;
    float regionCompactness(const Region&) const;
    float regionBoundaryMeasure(const Region&) const;

    using RegionAnalysisBase<T>::computeRAG;
    using RegionAnalysisBase<T>::_regions;

    VoxelMatrix<T>* _labelMatrix;
    const VoxelMatrix<T>* _valueMatrix;
    VoxelMatrix<T>* _outputMatrix;
};

#endif
#endif // LITE
