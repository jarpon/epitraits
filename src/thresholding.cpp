/*!
 * \class  Thresholding
 * \author Philippe Andrey (pa), INRA
 * \date   2005.04.24 - création (pa)
 * \brief  Binarisation par seuillage simple
 * \details
 * Cette classe réalise la binarisation d'une matrice-image autour
 * d'une valeur de seuil (0 par défaut, càd tous les pixels sont
 * sélectionnés).
 * Les valeurs d'avant-plan et d'arrière plan peuvent être spécifiées
 * (255 et 0, par défaut).
 * Par défaut, les pixels d'avant-plan sont ceux dont la valeur est
 * supérieure ou égale au seuil.
 * Utiliser #setForegroundIsAbove pour changer cette définition.
 *
 * Exemple d'utilisation:
 *
 * \code
 * Image image( "inputImage.tif" );
 * ImageMatrix<int> pixelMatrix( image );
 * Thresholding<int> thresholding;
 * thresholding.setForeground( 1 );
 * thresholding.setBackground( 0 );
 * thresholding.setThreshold( 100 );
 * thresholding.apply( pixelMatrix );
 * pixelMatrix.convertToImage( image );
 * image.writeFile( "binarisedImage.tif" );
 * \endcode
 *
****************************************************************/

#include "thresholding.h"

#include <curve.h>

#include <cstdlib>

//#define TRACE
#include <trace.h>

/*! Construit un opérateur de binarisation.
****************************************************************/
template<class T>
Thresholding<T>::Thresholding() :
  PixelMatrixProcessing<T>(),
  VoxelMatrixProcessing<T>()
{
  _threshold = 0;
  _background = 0;
  _foreground = 255;
  _foregroundIsAbove = true;
}

/*! Spécifie la valeur du seuil.
****************************************************************/
template<class T>
void Thresholding<T>::setThreshold(const T& threshold)
{
  _threshold = threshold;
}

/*! Retourne la valeur du seuil.
****************************************************************/
template<class T>
const T& Thresholding<T>::getThreshold() const
{
  return _threshold;
}

/*! Spécifie la valeur à attribuer aux pixels d'avant-plan.
****************************************************************/
template<class T>
void Thresholding<T>::setForeground(const T& foreground)
{
  _foreground = foreground;
}

/*! Spécifie la valeur à attribuer aux pixels d'arrière-plan.
****************************************************************/
template<class T>
void Thresholding<T>::setBackground(const T& background)
{
  _background = background;
}

/*! Retourne la valeur à attribuer aux pixels d'avant-plan.
****************************************************************/
template<class T>
const T& Thresholding<T>::getForeground() const
{
  return _foreground;
}

/*! Retourne la valeur à attribuer aux pixels d'arrière-plan.
****************************************************************/
template<class T>
const T& Thresholding<T>::getBackground() const
{
  return _background;
}

/*! Spécifie si les pixels dont la valeur est supérieure ou égale
 * au seuil constituent l'avant-plan ou l'arrière-plan.
****************************************************************/
template<class T>
void Thresholding<T>::setForegroundIsAbove(const bool foregroundIsAbove)
{
  _foregroundIsAbove = foregroundIsAbove;
}

/*! Retourne \c true si les pixels dont la valeur est supérieure
 * ou égale au seuil constituent l'avant-plan.
 * Retourne \c false s'ils constituent au contraire l'arrière-plan.
****************************************************************/
template<class T>
bool Thresholding<T>::foregroundIsAbove() const
{
  return _foregroundIsAbove;
}

/*! Binarise la matrice-image \c pixelMatrix au seuil qui a été spécifié.
 *
 * Par défaut:
 * - la valeur d'avant-plan est affectée aux pixels de valeur
 *   spérieure ou égale au seuil;
 * - la valeur d'arrière-plan est affectée aux autres pixels.
****************************************************************/
template<class T>
void Thresholding<T>::apply(PixelMatrix<T>& pixelMatrix)
{
  const int size1 = pixelMatrix.getSize1();
  const int size2 = pixelMatrix.getSize2();
  int i, j;
  T below, above;

  if ( foregroundIsAbove() )
  {
    above = getForeground();
    below = getBackground();
  }
  else
  {
    above = getBackground();
    below = getForeground();
  }

  this->notifyStarted();

  for (i = 0; i < size1; ++i)
  {
    for (j = 0; j < size2; ++j)
    {
      if ( pixelMatrix(i,j) < _threshold )
      {
        pixelMatrix(i,j) = below;
      }
      else
      {
        pixelMatrix(i,j) = above;
      }
    }
    this->notifyProgressed( i+1, size1 );
  }

  this->notifyFinished();
}

/*! Binarise l'image 3D \c voxelMatrix au seuil qui a été spécifié.
 *
 * Par défaut:
 * - la valeur d'avant-plan est affectée aux voxels de valeur
 *   supérieure ou égale au seuil;
 * - la valeur d'arrière-plan est affectée aux autres voxels.
****************************************************************/
template<class T>
void Thresholding<T>::apply(VoxelMatrix<T>& voxelMatrix)
{
  const int size1 = voxelMatrix.getSize1();
  const int size2 = voxelMatrix.getSize2();
  const int size3 = voxelMatrix.getSize3();
  int i, j, k;
  T below, above;

  if ( foregroundIsAbove() )
  {
    above = getForeground();
    below = getBackground();
  }
  else
  {
    above = getBackground();
    below = getForeground();
  }

  this->notifyStarted();

  for (i = 0; i < size1; ++i)
  {
    for (j = 0; j < size2; ++j)
    {
      for (k = 0; k < size3; ++k)
      {
        if ( voxelMatrix(i,j,k) < _threshold )
        {
          voxelMatrix(i,j,k) = below;
        }
        else
        {
          voxelMatrix(i,j,k) = above;
        }
      }
    }
    this->notifyProgressed( i+1, size1 );
  }

  this->notifyFinished();
}



/*! Binarise l'image 3D \c voxelMatrix au seuil qui a été spécifié.
 *
 * Par défaut:
 * - la valeur d'avant-plan est affectée aux voxels de valeur
 *   supérieure ou égale au seuil;
 * - la valeur d'arrière-plan est affectée aux autres voxels.
****************************************************************/
template<class T>
void Thresholding<T>::applyInverted(VoxelMatrix<T>& voxelMatrix)
{
  const int size1 = voxelMatrix.getSize1();
  const int size2 = voxelMatrix.getSize2();
  const int size3 = voxelMatrix.getSize3();
  int i, j, k;
  T below, above;

  if ( foregroundIsAbove() )
  {
    above = getForeground();
    below = getBackground();
  }
  else
  {
    above = getBackground();
    below = getForeground();
  }

  this->notifyStarted();

  for (i = 0; i < size1; ++i)
  {
    for (j = 0; j < size2; ++j)
    {
      for (k = 0; k < size3; ++k)
      {
        if ( voxelMatrix(i,j,k) > _threshold )
        {
          voxelMatrix(i,j,k) = below;
        }
        else
        {
          voxelMatrix(i,j,k) = above;
        }
      }
    }
    this->notifyProgressed( i+1, size1 );
  }

  this->notifyFinished();
}

/*! Binarise l'image 3D \c voxelMatrix au seuil qui a été spécifié.
 *
 * Par défaut:
 * - la valeur d'avant-plan est affectée aux voxels de valeur
 *   supérieure ou égale au seuil;
 * - la valeur d'arrière-plan est affectée aux autres voxels.
****************************************************************/
template<class T>
void Thresholding<T>::applyAlternative(VoxelMatrix<T>& voxelMatrix)
{
  const int size1 = voxelMatrix.getSize1();
  const int size2 = voxelMatrix.getSize2();
  const int size3 = voxelMatrix.getSize3();
  int i, j, k;

  this->notifyStarted();

  for (i = 0; i < size1; ++i)
  {
    for (j = 0; j < size2; ++j)
    {
      for (k = 0; k < size3; ++k)
      {
        if ( voxelMatrix(i,j,k) < _threshold )
        {
          voxelMatrix(i,j,k) = 0;
        }
      }
    }
    this->notifyProgressed( i+1, size1 );
  }

  this->notifyFinished();
}


/*! Generates the binary mask of the level set corresponding to \c value.
 *
 * All pixels in \c pixelMatrix equal to \c value are assigned the
 * foreground value.
 * The other pixels are assigned the background value.
****************************************************************/
template<class T>
void Thresholding<T>::levelSetMask(PixelMatrix<T>& pixelMatrix, const T value)
{
  const int size1 = pixelMatrix.getSize1();
  const int size2 = pixelMatrix.getSize2();
  int i, j;

  for (i = 0; i < size1; ++i)
  {
    for (j = 0; j < size2; ++j)
    {
      pixelMatrix(i,j) = pixelMatrix(i,j) == value? _foreground: _background;
    }
  }
}

/*! Generates the binary mask of the level set corresponding to \c value.
 *
 * All voxels in \c voxelMatrix equal to \c value are assigned the
 * foreground value.
 * The other voxels are assigned the background value.
****************************************************************/
template<class T>
void Thresholding<T>::levelSetMask(VoxelMatrix<T>& voxelMatrix, const T value)
{
  const int size3 = voxelMatrix.getSize3();
  for (int k = 0; k < size3; ++k)
    levelSetMask( voxelMatrix[k], value );
}

/*! Détermine les contours des régions dans la matrice-image \c im.
****************************************************************/
template<class T>
void Thresholding<T>::contours(PixelMatrix<T>& im) const
{
  Vector<T> prev_row( im.getSize2() ); // buffer ligne précédente
  Vector<T> curr_row( im.getSize2() ); // buffer ligne courante
  bool p0, p1, p2, p3, p4, p5, p6, p7, p8; // drapeaux voisinage == fg
  T fg = getForeground();
  T bg = getBackground();
  int i, j;

  curr_row = im[0];

  for (i = 1; i < im.getSize1()-1; ++i)
  {
    const Vector<T>& next_row = im[i+1];

    prev_row.swap( curr_row );
    curr_row = im[i];

    p2 = prev_row[0] == fg; p3 = prev_row[1] == fg;
    p0 = curr_row[0] == fg; p5 = curr_row[1] == fg;
    p7 = next_row[0] == fg; p8 = next_row[1] == fg;

    for (j = 1; j < im.getSize2()-1; ++j)
    {
      p1 = p2; p2 = p3; p3 = prev_row[j+1] == fg;
      p4 = p0; p0 = p5; p5 = curr_row[j+1] == fg;
      p6 = p7; p7 = p8; p8 = next_row[j+1] == fg;

      if ( p0 && p1 && p2 && p3 && p4 && p5 && p6 && p7 && p8 )
      {
        im(i,j) = bg;
      }
    }
  }
}

/*! Extrait le contour \c curve qui passe par le point \c i_start,
 * \c j_start dans la matrice-image \c pixelMatrix.
 *
 * Le point de départ doit être à la valeur d'avant-plan.
 *
 * Si \c decimationLevel > 0, la décimation est activée.
 * Le contour extrait contient d'autant moins de points que la
 * valeur de \c decimationLevel est élevée.
 *
 * Référence: T. Pavlidis, "Algorithms for Graphics and Image
 * Processing". Springer Verlag, Berlin-Heidelberg, 1982.
 * [Algorithme 7.1, page 143].
****************************************************************/
template<class T>
void Thresholding<T>::extractCurve(
  Curve<int>& curve,
  const int i_start,
  const int j_start,
  const PixelMatrix<T>& pixelMatrix,
  const int decimationLevel) const
{
  typedef enum { NORTH, SOUTH, EAST, WEST } Direction;
  int i, i_prev, di, di_prev = 100, ilim = pixelMatrix.getSize1()-1;
  int j, j_prev, dj, dj_prev = 100, jlim = pixelMatrix.getSize2()-1;
  Vectori vertex( curve.getSpace() );
  bool found, starting_point = true;
  bool decimationEnabled = decimationLevel > 0;
  Direction direction = SOUTH;
  T fg = _foreground;
  int loop_count;

  vertex[X] = i = i_start;
  vertex[Y] = j = j_start;
  curve.append( vertex );

  while ( starting_point || i != i_start || j != j_start )
  {
    i_prev = i;
    j_prev = j;
    found = false;
    loop_count = 0;
    while ( !found && loop_count < 3 )
    {
      switch( direction )
      {
        case NORTH:
          if ( i > 0 )
          {
            if ( j < jlim && pixelMatrix(i-1,j+1) == fg ) i--, j++, direction = EAST, found = true;
            else if ( pixelMatrix(i-1,j) == fg ) i--, found = true;
            else if ( j > 0 && pixelMatrix(i-1,j-1) == fg ) i--, j--, found = true;
            else direction = WEST;
          }
          else
          {
            direction = WEST;
          }
          break;
        case SOUTH:
          if ( i < ilim )
          {
            if ( j > 0 && pixelMatrix(i+1,j-1) == fg ) i++, j--, direction = WEST, found = true;
            else if ( pixelMatrix(i+1,j) == fg ) i++, found = true;
            else if ( j < jlim && pixelMatrix(i+1,j+1) == fg ) i++, j++, found = true;
            else direction = EAST;
          }
          else
          {
            direction = EAST;
          }
          break;
        case EAST:
          if ( j < jlim )
          {
            if ( i < ilim && pixelMatrix(i+1,j+1) == fg ) i++, j++, direction = SOUTH, found = true;
            else if ( pixelMatrix(i,j+1) == fg ) j++, found = true;
            else if ( i > 0 && pixelMatrix(i-1,j+1) == fg ) i--, j++, found = true;
            else direction = NORTH;
          }
          else
          {
            direction = NORTH;
          }
          break;
        case WEST:
          if ( j > 0 )
          {
            if ( i > 0 && pixelMatrix(i-1,j-1) == fg ) i--, j--, direction = NORTH, found = true;
            else if ( pixelMatrix(i,j-1) == fg ) j--, found = true;
            else if ( i < ilim && pixelMatrix(i+1,j-1) == fg ) i++, j--, found = true;
            else direction = SOUTH;
          }
          else
          {
            direction = SOUTH;
          }
          break;
      }
      loop_count++;
    }

    starting_point = false;

#if 0
    if ( i != int(i_start) || j != int(j_start) )
    {
      vertex[X] = i;
      vertex[Y] = j;
      curve.append( vertex );
    }
  }
#else
    di = i - i_prev;
    dj = j - j_prev;
    i_prev = i;
    j_prev = j;

    if ( di == di_prev && dj == dj_prev )
    {
      // la direction du contour n'a pas changé: il suffit
      // alors de déplacer le dernier point de la courbe
      curve[curve.getSize()-1][X] = i;
      curve[curve.getSize()-1][Y] = j;
    }
    else if ( (di != di_prev && dj != dj_prev) || !decimationEnabled )
    {
      // la direction du contour a changé: on ajoute un point
      di_prev = di;
      dj_prev = dj;
      if ( i != int(i_start) || j != int(j_start) )
      {
        vertex[X] = i;
        vertex[Y] = j;
        curve.append( vertex );
      }
    }
    else // di ou dj n'a pas changé, et la décimation est activée
    {
      di_prev = di;
      dj_prev = dj;

      Vectori& v1 = curve[curve.getSize()-1]; // dernier point rentré
      Vectori& v2 = curve[curve.getSize()-2]; // avant-dernier point
      int i1 = v1[X], j1 = v1[Y];
      int i2 = v2[X], j2 = v2[Y];
      int area; // double de l'aire du triangle

      area = j*i1 - i*j1 + j1*i2 - i1*j2 + j2*i - i2*j;

      if ( abs(area) > decimationLevel ) // ajout d'un nouveau point
      {
        if ( i != int(i_start) || j != int(j_start) )
        {
          vertex[X] = i;
          vertex[Y] = j;
          curve.append( vertex );
        }
      }
      else // on déplace le dernier point
      {
        v1[X] = i;
        v1[Y] = j;
      }
    }
  }
#endif
}

template class Thresholding<int>;
template class Thresholding<float>;
template class Thresholding<double>;
template class Thresholding<long double>;
