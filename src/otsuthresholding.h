#ifndef OTSUTHRESHOLDING_H
#define OTSUTHRESHOLDING_H

#include "thresholding.h"

template<class T>
class OtsuThresholding : public Thresholding<T>
{
  public:

    OtsuThresholding();

    void setBins(const T&,const T&,const int);

    T computeThreshold(const PixelMatrix<T>&) const;
    T computeThreshold(const VoxelMatrix<T>&) const;
    int computeThreshold(const Vector<unsigned int>&) const;
    Vector<double> computeCriterion(const Vector<unsigned int>&) const;

    void apply(PixelMatrix<T>&);
    void apply(VoxelMatrix<T>&);
    //void applyAlternative(VoxelMatrix<T>&);

  private:

    T _firstBin;  //!< Borne inférieure de la première classe
    T _binWidth;  //!< Largeur des classes de l'histogramme
    int _numBins; //!< Nombre de classes de l'histogramme
};

#endif
