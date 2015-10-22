#ifndef REGIONANALYSIS2D_H
#define REGIONANALYSIS2D_H

//#include "pixelmatrix.h"
//#include "regionanalysisbase.h"
#include <pixelmatrix.h>
#include <regionanalysisbase.h>

template<class T>
class RegionAnalysis2D : public RegionAnalysisBase<T>
{
  public:

    RegionAnalysis2D();
    ~RegionAnalysis2D();

    int dimension() const { return 2; }

    bool hasLabelMatrix() const;
    void setLabelMatrix(PixelMatrix<T>&);
    const PixelMatrix<T>& getLabelMatrix() const;

    bool hasValueMatrix() const;
    void setValueMatrix(const PixelMatrix<T>&);
    const PixelMatrix<T>& getValueMatrix() const;

    bool hasOutputMatrix() const;
    void setOutputMatrix(PixelMatrix<T>&);
    PixelMatrix<T>& getOutputMatrix() const;

    void run();

    Vector<T> regionValues(const Region&) const;

    void outputFillRegion(const Region&,const T);

    using RegionAnalysisBase<T>::numRegions;
    using RegionAnalysisBase<T>::getRegion;
    using RegionAnalysisBase<T>::getBoundaryValue;
    using RegionAnalysisBase<T>::isRegionLabel;

  private:

    RegionAnalysis2D(const RegionAnalysis2D&);
    RegionAnalysis2D& operator=(const RegionAnalysis2D&);

    void computeRAGImplicitBoundaryMode();
    void computeRAGExplicitBoundaryMode();
    using RegionAnalysisBase<T>::computeRAG;

    virtual Vector<float> spatialCalibration() const;
    float regionCompactness(const Region&) const;
    float regionBoundaryMeasure(const Region&) const;

    PixelMatrix<T> *_labelMatrix;
    const PixelMatrix<T>* _valueMatrix;
    PixelMatrix<T>* _outputMatrix;
};

#endif // REGIONANALYSIS2D_H
