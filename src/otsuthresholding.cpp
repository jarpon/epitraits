/*!
 * \class  OtsuThresholding
 * \author Philippe Andrey (pa), INRA
 * \date   2005.04.24 - creation (pa)
 * \brief  Image binarization using Otsu's histogram thresholding
 * \details
 * Cette classe implémente la méthode d'Otsu (1979) de binarisation
 * automatique d'une image.
 * Cette méthode consiste à déterminer le seuil qui minimise la somme
 * des variances intra-classes.
 *
 * Exemple d'utilisation:
 *
 * \code
 * Image image( "inputImage.tif" );
 * PixelMatrixi pixelMatrix( image );
 * OtsuThresholdingi otsuThresholding;
 * otsuThresholding.apply( pixelMatrix );
 * pixelMatrix.convertToImage( image );
 * image.writeFile( "binarisedImage.tif" );
 * cout << "Image thresholded at " << otsuThresholding.getThreshold() << endl;
 * \endcode
 *
 * La méthode #computeThreshold fournit le seuil sans réaliser la
 * binarisation.
 *
 * Par défaut, les calculs sont fondés sur l'histogramme des 256
 * valeurs allant de 0 à 255.
 * La méthode #setBins permet de modifier les caractéristiques de
 * l'histogramme utilisé.
****************************************************************/

#include "otsuthresholding.h"

#include <maths.h>

/*! Constructor.
****************************************************************/
template<class T>
OtsuThresholding<T>::OtsuThresholding() : Thresholding<T>()
{
  _firstBin = -0.5;
  _binWidth = 1.0;
  _numBins = 256;
}

template<>
OtsuThresholding<int>::OtsuThresholding() : Thresholding<int>()
{
  _firstBin = 0;
  _binWidth = 1;
  _numBins = 256;
}

/*! Spécifie les caractéristiques de l'histogramme sur lequel la
 * méthode doit se fonder pour déterminer le seuil optimal.
****************************************************************/
template<class T>
void OtsuThresholding<T>::setBins(
  const T& firstBin,
  const T& binWidth,
  const int numBins)
{
  _firstBin = firstBin;
  _binWidth = binWidth;
  _numBins = numBins;
}

/*! Computes and applies Otsu's threshold to \c pixelMatrix.
****************************************************************/
template<class T>
void OtsuThresholding<T>::apply(PixelMatrix<T>& pixelMatrix)
{
  this->_threshold = computeThreshold( pixelMatrix );

  Thresholding<T>::apply( pixelMatrix );
}

/*! Computes and applies Otsu's threshold to \c voxelMatrix.
****************************************************************/
template<class T>
void OtsuThresholding<T>::apply(VoxelMatrix<T>& voxelMatrix)
{
  this->_threshold = computeThreshold( voxelMatrix );

  Thresholding<T>::apply( voxelMatrix );
}

/*! Computes and returns the vector of Otsu's criterion values
 * for all possible thresholds of the specified histogram.
****************************************************************/
template<class T>
Vector<double> OtsuThresholding<T>::computeCriterion(const Vector<unsigned int>& histogram) const
{
  const int numBins = histogram.getSize();
  Vector<double> criterion( numBins );
  int n0 = 0, n1 = 0; // class sizes
  int s0 = 0, s1 = 0; // sum of values
  double m0, m1; // mean values
  double p0, p1; // class proportions
  int g;

  for (g = 0; g < numBins; ++g)
  {
    criterion[g] = 0;
    n1 += histogram[g];
    s1 += g * histogram[g];
  }

  for (g = 1; g < numBins; ++g)
  {
    n0 += histogram[g-1];
    n1 -= histogram[g-1];
    s0 += (g-1) * histogram[g-1];
    s1 -= (g-1) * histogram[g-1];
    if ( n0 > 0 && n1 > 0 )
    {
      m0 = static_cast<double>(s0) / n0;
      m1 = static_cast<double>(s1) / n1;
      p0 = static_cast<double>(n0) / (n0+n1);
      p1 = static_cast<double>(n1) / (n0+n1);
      criterion[g] = p0 * p1 * Maths::sqr( m0-m1 );
    }
  }

  return criterion;
}

/*! Computes and returns the optimal threshold for the specified histogram.
****************************************************************/
template<class T>
int OtsuThresholding<T>::computeThreshold(const Vector<unsigned int>& histogram) const
{
  const Vector<double> criterion = computeCriterion( histogram );

  return criterion.maxPos();

#if 0
  double s0, ss0, m0, v0, s1, ss1, m1, v1, vmin, vgng;
  // initialisation des statistiques de la classe 0
  m0 = s0 = ss0 = v0 = 0;

  // initialisation des statistiques de la classe 1
  s1 = ss1 = 0;
  for (g = 0; g < histogram.getSize(); ++g)
  {
    n1 += histogram[g];
    s1 += g * histogram[g];
    ss1 += g * g * histogram[g];
  }
  m1 = s1 / n1;
  v1 = ss1 / n1 - m1 * m1;

  vmin = v1;
  gmin = 0;

  for (g = 0; g < histogram.getSize(); ++g)
  {
    n0 += histogram[g];
    n1 -= histogram[g];
    vgng = g * histogram[g];
    s0 += vgng;
    s1 -= vgng;
    ss0 += g * vgng;
    ss1 -= g * vgng;

    if ( n0 && n1 )
    {
      m0 = s0 / n0;
      m1 = s1 / n1;
      v0 = ss0 / n0 - m0 * m0;
      v1 = ss1 / n1 - m1 * m1;
      if ( v0 + v1 < vmin )
      {
        vmin = v0 + v1;
        gmin = g + 1;
      }
    }
  }

    return gmin;
#endif
}

/*! Computes and returns the optimal threshold for \c pixelMatrix.
****************************************************************/
template<class T>
T OtsuThresholding<T>::computeThreshold(const PixelMatrix<T>& pixelMatrix) const
{
  const Vector<unsigned int> histogram = pixelMatrix.histogram( _firstBin, _binWidth, _numBins );
  const int g = computeThreshold( histogram );

  return _firstBin + g * _binWidth;
}

/*! Computes and returns the optimal threshold for \c voxelMatrix.
****************************************************************/
template<class T>
T OtsuThresholding<T>::computeThreshold(const VoxelMatrix<T>& voxelMatrix) const
{
  const Vector<unsigned int> histogram = voxelMatrix.histogram( _firstBin, _binWidth, _numBins );
  const int g = computeThreshold( histogram );

  return _firstBin + g * _binWidth;
}

template class OtsuThresholding<int>;
template class OtsuThresholding<float>;
template class OtsuThresholding<double>;
template class OtsuThresholding<long double>;
