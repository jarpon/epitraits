#ifndef THRESHOLDING_H
#define THRESHOLDING_H

template<class T>
class Curve;

#include "pixelmatrixprocessing.h"
#include "voxelmatrixprocessing.h"

template<class T>
class Thresholding :
  public PixelMatrixProcessing<T>,
  public VoxelMatrixProcessing<T>
{
  public:

    Thresholding();

    void setThreshold(const T&);
    const T& getThreshold() const;

    void setForeground(const T&);
    void setBackground(const T&);
    const T& getForeground() const;
    const T& getBackground() const;

    void setForegroundIsAbove(const bool);
    bool foregroundIsAbove() const;

    void apply(PixelMatrix<T>&);
    void apply(VoxelMatrix<T>&);
    void levelSetMask(PixelMatrix<T>&,T);
    void levelSetMask(VoxelMatrix<T>&,T);

    void applyAlternative(VoxelMatrix<T>&);
    void applyInverted(VoxelMatrix<T>&);


    void contours(PixelMatrix<T>&) const;
    void extractCurve(Curve<int>&,
      const int,const int,
      const PixelMatrix<T>&,const int=0) const;

  protected:

    mutable T _threshold;   //!< Valeur de seuil

  private:

    T _foreground;  //!< Valeur d'avant-plan
    T _background;  //!< Valeur d'arrière-plan
    bool _foregroundIsAbove; //!< Sens du seuillage
};

#endif
