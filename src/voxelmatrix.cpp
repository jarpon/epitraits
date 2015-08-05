/*!
 * \class  VoxelMatrix
 * \author Philippe Andrey (pa), INRA
 * \date   2005.04.08 - création de VolumeMatrix (pa)
 * \date   2007.12.18 - nouveau nom VoxelMatrix (pa)
 * \brief  Matrice de voxels (image 3D)
 * \details
 * Cette classe permet la représentation et la manipulation d'images
 * tridimensionnelles sous la forme d'empilements de matrices de pixels
 * (voir classe PixelMatrix).
 * Comme illustré ci-dessous, l'axe d'empilement des matrices de pixels
 * est l'axe des Z:
 *
 * \image html volumematrix.png
 *
 * Les conventions utilisées sont les suivantes:
 * - la dimension 1 correspond à l'axe des X (indice i)
 * - la dimension 2 correspond à l'axe des Y (indice j)
 * - la dimension 3 correspond à l'axe des Z (indice k)

 * Avec ces notations, le repère d'axes XYZ est direct.
****************************************************************/

#include "voxelmatrix.h"

#include <filereaderror.h>
#include <filewriteerror.h>
#include <sizeerror.h>
#include <typeerror.h>
#include <valueerror.h>
#include <holesfilling.h>

#include <fstream>
#include <typeinfo>

//#define TRACE
#include <trace.h>

/*! Crée une matrice de taille nulle.
****************************************************************/
template<class T>
VoxelMatrix<T>::VoxelMatrix() : _n(0), _p(0), _h(0), _planes(0)
{
}

/*! Crée un volume de \c size3 plans dont chacun est de taille
 * \c size1 lignes par \c size2 colonnes.
 *
 * Les valeurs des voxels de sont pas initialisées.
****************************************************************/
template<class T>
VoxelMatrix<T>::VoxelMatrix(
  const int size1,
  const int size2,
  const int size3) : _n(0), _p(0), _h(0), _planes(0)
{
  setSize( size1, size2, size3 );
}

/*! Crée une copie de la matrice \c voxelMatrix.
****************************************************************/
template<class T>
VoxelMatrix<T>::VoxelMatrix(const VoxelMatrix& voxelMatrix) :
  _n(0), _p(0), _h(0), _planes(0)
{
  setSize(
    voxelMatrix.getSize1(),
    voxelMatrix.getSize2(),
    voxelMatrix.getSize3() );

  for (int k = 0; k < _h; ++k)
  {
    _planes[k] = voxelMatrix._planes[k];
  }

  _voxelCalibration = voxelMatrix._voxelCalibration;
}

/*! Crée une matrice à partir des données stockées dans le fichier
 * nommé \c fileName.
 *
 * \exception FileReadError (propagée) Le fichier n'est pas accessible
 * en lecture
****************************************************************/
template<class T>
VoxelMatrix<T>::VoxelMatrix(const string& fileName) :
  _n(0), _p(0), _h(0), _planes(0)
{
  try
  {
    load( fileName );
  }

  catch( Exception& error )
  {
    error.addCall( "VoxelMatrix<T>::VoxelMatrix(const string&)" );
    throw error;
  }
}

/*! Détruit ce volume.
****************************************************************/
template<class T>
VoxelMatrix<T>::~VoxelMatrix()
{
  delete[] _planes;
}

/*! Spécifie la taille de cette matrice de voxels.
 *
 * Le nombre de plans est donné par \c size3.
 * Les nombres de lignes et de colonnes de chacun des plans sont
 * donnés respectivement par \c size1 et \c size2.
****************************************************************/
template<class T>
void VoxelMatrix<T>::setSize(
  const int size1,
  const int size2,
  const int size3)
{
  if ( _n != size1 || _p != size2 || _h != size3 )
  {
    delete[] _planes;
    _n = size1;
    _p = size2;
    _h = size3;
    _planes = _h > 0? new PixelMatrix<T>[_h]: 0;
    for (int k = 0; k < _h; ++k)
    {
      _planes[k].setSize( _n, _p );
    }
  }
}

/*! Spécifie la taille de cette matrice de voxels.
****************************************************************/
template<class T>
void VoxelMatrix<T>::setSize(const Vector<int>& size)
{
  setSize( size[0], size[1], size[2] );
}

/*! Retourne la taille de cette matrice le long de la première
 * dimension (càd dimension X, ce qui correspond au nombre de
 * lignes des plans images).
****************************************************************/
template<class T>
int VoxelMatrix<T>::getSize1() const
{
  return _n;
}

/*! Retourne la taille de cette matrice le long de la deuxième
 * dimension (càd dimension Y, ce qui correspond au nombre de
 * colonnes des plans images).
****************************************************************/
template<class T>
int VoxelMatrix<T>::getSize2() const
{
  return _p;
}

/*! Retourne la taille de cette matrice le long de la troisième
 * dimension (càd dimension Z, ce qui correspond au nombre de
 * plans images empilés).
****************************************************************/
template<class T>
int VoxelMatrix<T>::getSize3() const
{
  return _h;
}

/*! Retourne la taille de cette matrice de voxels le long de la ième
 * dimension (\c i = 1, 2 ou 3).
 *
 * \sa getSize1, getSize2, getSize3
****************************************************************/
template<class T>
int VoxelMatrix<T>::getSize(const int i) const
{
  switch ( i )
  {
    case 1: return _n;
    case 2: return _p;
    case 3: return _h;
    default:
    {
      ValueError valueError;
      valueError.setWhere( "unsigned int VoxelMatrix<T>::getSize(...) const" );
      valueError.setWhat( "Invalid argument value (neither 1, 2 or 3)" );
      throw valueError;
    }
  }
}

/*! Retourne les trois dimensions de cette matrice dans un vecteur.
****************************************************************/
template<class T>
Vector<int> VoxelMatrix<T>::getSize() const
{
  Vector<int> v( 3 );
  v[0] = getSize1();
  v[1] = getSize2();
  v[2] = getSize3();
  return v;
}

/*! Retourne \c true si cette matrice et la matrice \c voxelMatrix
 * ont les mêmes dimensions.
****************************************************************/
template<class T>
inline bool VoxelMatrix<T>::sameSizeAs(const VoxelMatrix& voxelMatrix) const
{
  return
    _n == voxelMatrix._n &&
    _p == voxelMatrix._p &&
      _h == voxelMatrix._h;
}

/*! Sets the spatial calibration of the voxels in this image.
 *
 * \sa getVoxelCalibration
****************************************************************/
template<class T>
void VoxelMatrix<T>::setVoxelCalibration(const VoxelCalibration& voxelCalibration)
{
  _voxelCalibration = voxelCalibration;
}

/*! Returns the spatial calibration of the voxels in this image.
 *
 * \sa setVoxelCalibration
****************************************************************/
template<class T>
const VoxelCalibration& VoxelMatrix<T>::getVoxelCalibration() const
{
  return _voxelCalibration;
}

/*! Affecte à cette matrice de voxels les dimensions et les valeurs
 * des voxels de la matrice \c voxelMatrix.
****************************************************************/
template<class T>
VoxelMatrix<T>& VoxelMatrix<T>::operator=(const VoxelMatrix& voxelMatrix)
{
  if ( this != &voxelMatrix )
  {
    setSize( voxelMatrix.getSize() );

    for (int k = 0; k < _h; ++k)
    {
      _planes[k] = voxelMatrix._planes[k];
    }

    _voxelCalibration = voxelMatrix._voxelCalibration;
  }

  return *this;
}

/*! Echange les dimensions ainsi que les valeurs des voxels entre
 * cette matrice de voxels et la matrice \c voxelMatrix.
****************************************************************/
template<class T>
void VoxelMatrix<T>::swap(VoxelMatrix& voxelMatrix)
{
  const int n0 = _n;
  const int p0 = _p;
  const int h0 = _h;
  PixelMatrix<T>* planes0 = _planes;
  VoxelCalibration voxelCalibration = _voxelCalibration;

  _n = voxelMatrix._n;
  _p = voxelMatrix._p;
  _h = voxelMatrix._h;
  _planes = voxelMatrix._planes;
  _voxelCalibration = voxelMatrix._voxelCalibration;

  voxelMatrix._n = n0;
  voxelMatrix._p = p0;
  voxelMatrix._h = h0;
  voxelMatrix._planes = planes0;
  voxelMatrix._voxelCalibration = voxelCalibration;
}

/*! Retourne \c true si cette matrice de voxels est identique
 * à la matrice \c voxelMatrix.
****************************************************************/
template<class T>
bool VoxelMatrix<T>::operator==(const VoxelMatrix& voxelMatrix)
{
  if ( this != &voxelMatrix )
  {
    for (int d = 1; d <= 3; ++d)
    {
      if ( getSize(d) != voxelMatrix.getSize(d) )
      {
        return false;
      }
    }

    for (int k = 0; k < _h; ++k)
    {
      if ( _planes[k] != voxelMatrix._planes[k] )
      {
        return false;
      }
    }
  }

  return true;
}

/*! Retourne \c true si cette matrice de voxels diffère de la
 * matrice \c voxelMatrix.
****************************************************************/
template<class T>
bool VoxelMatrix<T>::operator!=(const VoxelMatrix& voxelMatrix)
{
  return !( *this == voxelMatrix );
}

/*! Met tous les voxels de cette matrice à 0.
 *
 * \sa setOnes, fill
****************************************************************/
template<class T>
void VoxelMatrix<T>::setZeros()
{
  for (int k = 0; k < _h; ++k)
  {
    _planes[k].setZeros();
  }
}

/*! Met tous les voxels de cette matrice à 1.
 *
 * \sa setZeros, fill
****************************************************************/
template<class T>
void VoxelMatrix<T>::setOnes()
{
  for (int k = 0; k < _h; ++k)
  {
    _planes[k].setOnes();
  }
}

/*! Attribue à tous les voxels de cette matrice la valeur \c value.
 *
 * \sa setZeros, setOnes
****************************************************************/
template<class T>
void VoxelMatrix<T>::fill(const T& value)
{
  for (int k = 0; k < _h; ++k)
  {
    _planes[k].fill( value );
  }
}

/*! Applique à tous les voxels de cette matrice la fonction \c f.
****************************************************************/
template<class T>
void VoxelMatrix<T>::apply(T(*f)(T))
{
  for (int k = 0; k < _h; ++k)
  {
    _planes[k].apply( f );
  }
}

/*! Inverse l'ordre des plans images de cette matrice de voxels.
****************************************************************/
template<class T>
void VoxelMatrix<T>::reverse()
{
  for (int kInf = 0, kSup = _h-1; kInf < kSup; kInf++, kSup--)
  {
    _planes[kInf].swap( _planes[kSup] );
  }
}

/*! Retourne une copie du bloc compris entre les lignes \c i1 à \c i2,
 * les colonnes \c j1 à \c j2 et les plans \c k1 à \c k2.
 *
 * Attention, la validité du domaine spécifié n'est pas testée.
 *
 * \sa paste
****************************************************************/
template<class T>
VoxelMatrix<T> VoxelMatrix<T>::copy(
  const int i1, const int i2,
  const int j1, const int j2,
  const int k1, const int k2) const
{
  VoxelMatrix voxelMatrix( i2-i1+1, j2-j1+1, k2-k1+1 );

  for (int k = k1; k <= k2; ++k)
  {
    voxelMatrix._planes[k-k1] = _planes[k].copy( i1, i2, j1, j2 );
  }

  voxelMatrix._voxelCalibration = _voxelCalibration;

  return voxelMatrix;
}

/*! Colle la matrice de voxels \c voxelMatrix à la position
 * \c i, \c j, \c k dans cette matrice.
 *
 * Attention, la validité des arguments n'est pas testée.
 *
 * \sa copy
****************************************************************/
template<class T>
void VoxelMatrix<T>::paste(
  const int i,
  const int j,
  const int k,
  const VoxelMatrix& voxelMatrix)
{
  const int pmax = k + voxelMatrix._h - 1;

  for (int p = k; p <= pmax; ++p)
  {
    _planes[p].paste( i, j, voxelMatrix._planes[p] );
  }
}

/*! Concatène les valeurs des voxels de cette matrice et retourne
 * le vecteur résultant.
 *
 * La concaténation est effectuée plan par plan et ligne à ligne
 * au sein de chaque plan.
****************************************************************/
template<class T>
Vector<T> VoxelMatrix<T>::cat() const
{
  Vector<T> v( getSize().prod() );

  for (int k = 0; k < _h; ++k)
  {
    v.paste( k*getSize1()*getSize2(), _planes[k].cat() );
  }

  return v;
}

/*! Ajoute à chacun des voxels de cette matrice le voxel correspondant
 * dans la matrice \c voxelMatrix.
 *
 * \exception SizeError les dimensions de cette matrice et de
 * \c volumeMatrix diffèrent.
****************************************************************/
template<class T>
void VoxelMatrix<T>::operator+=(const VoxelMatrix& voxelMatrix)
{
  if ( !sameSizeAs(voxelMatrix) )
  {
    SizeError sizeError;
    sizeError.setWhere( "void VoxelMatrix<T>::operator+=(const VoxelMatrix&)" );
    sizeError.setWhat( "This matrix and matrix argument differ in size" );
    throw sizeError;
  }

  for (int k = 0; k < _h; ++k)
  {
    _planes[k] += voxelMatrix._planes[k];
  }
}

/*! Soustrait à chacun des voxels de cette matrice le voxel correspondant
 * dans la matrice \c voxelMatrix.
 *
 * \exception SizeError les dimensions de cette matrice et de
 * \c volumeMatrix diffèrent.
****************************************************************/
template<class T>
void VoxelMatrix<T>::operator-=(const VoxelMatrix& voxelMatrix)
{
  if ( !sameSizeAs(voxelMatrix) )
  {
    SizeError sizeError;
    sizeError.setWhere( "void VoxelMatrix<T>::operator-=(const VoxelMatrix&)" );
    sizeError.setWhat( "This matrix and matrix argument differ in size" );
    throw sizeError;
  }

  for (int k = 0; k < _h; ++k)
  {
    _planes[k] -= voxelMatrix._planes[k];
  }
}

/*! Multiplie chacun des voxels de cette matrice par le voxel correspondant
 * dans la matrice \c voxelMatrix.
 *
 * \exception SizeError les dimensions de cette matrice et de
 * \c volumeMatrix diffèrent.
****************************************************************/
template<class T>
void VoxelMatrix<T>::operator*=(const VoxelMatrix& voxelMatrix)
{
  if ( getSize() != voxelMatrix.getSize() )
  {
    SizeError sizeError;
    sizeError.setWhere( "void VoxelMatrix<T>::operator*=(const VoxelMatrix&)" );
    sizeError.setWhat( "This matrix and matrix argument differ in size" );
    throw sizeError;
  }

  const int size1 = getSize1();
  const int size2 = getSize2();
  const int size3 = getSize3();
  int i, j, k;

  for (k = 0; k < size3; ++k)
  {
    for (i = 0; i < size1; ++i)
    {
      for (j = 0; j < size2; ++j)
      {
        _planes[k](i,j) *= voxelMatrix(i,j,k);
      }
    }
  }
}


/*! Ajoute \c value à chacun des voxels de cette matrice.
****************************************************************/
template<class T>
void VoxelMatrix<T>::operator+=(const T& value)
{
  for (int k = 0; k < _h; ++k)
  {
    _planes[k] += value;
  }
}

/*! Retranche \c value à chacun des voxels de cette matrice.
****************************************************************/
template<class T>
void VoxelMatrix<T>::operator-=(const T& value)
{
  for (int k = 0; k < _h; ++k)
  {
    _planes[k] -= value;
  }
}

/*! Multiplie par \c value chacun des voxels de cette matrice.
****************************************************************/
template<class T>
void VoxelMatrix<T>::operator*=(const T& value)
{
  for (int k = 0; k < _h; ++k)
  {
    _planes[k] *= value;
  }
}

/*! Divise par \c value chacun des voxels de cette matrice.
****************************************************************/
template<class T>
void VoxelMatrix<T>::operator/=(const T& value)
{
  for (int k = 0; k < _h; ++k)
  {
    _planes[k] /= value;
  }
}

/*! Détermine et retourne l'histogramme des valeurs des voxels
 * de cette matrice.
****************************************************************/
template<class T>
Vectorui VoxelMatrix<T>::histogram(
  const T firstBin, //!< Borne gauche de la première classe
  const T binWidth, //!< Largeur de chacune des classes
  const unsigned int numBins //!< Nombre de classes à utiliser
  ) const
{
  Vectorui h( numBins );

  h.setZeros();

  for (int k = 0; k < _h; ++k)
  {
    h += _planes[k].histogram( firstBin, binWidth, numBins );
  }

  return h;
}

/*! Détermine et retourne l'histogramme des valeurs des voxels
 * de cette matrice qui appartiennent au masque \c mask.
 *
 * \exception SizeError la taille du masque diffère de celle de
 * volume.
****************************************************************/
template<class T>
Vectorui VoxelMatrix<T>::histogram(
  const T firstBin, //!< Borne gauche de la première classe
  const T binWidth, //!< Largeur de chacune des classes
  const unsigned int numBins, //!< Nombre de classes à utiliser
  const VoxelMatrix<T>& mask  //!< Masque à utiliser
  ) const
{
  if ( getSize() != mask.getSize() )
  {
    SizeError sizeError;
    sizeError.setWhere( "Vectorui VoxelMatrix<T>::histogram(...)" );
    sizeError.setWhat( "Incompatible mask size" );
    throw sizeError;
  }

  Vectorui h( numBins );

  h.setZeros();

  for (int k = 0; k < _h; ++k)
  {
    h += _planes[k].histogram( firstBin, binWidth, numBins, mask[k] );
  }

  return h;
}

/*! Détermine et retourne une matrice de pixels contenant, en chaque
 * position \c i, \c j, la somme des valeurs des voxels au travers
 * de l'empilement des plans de cette matrice.
 *
 * \exception SizeError Cette matrice est de taille nulle
****************************************************************/
template<class T>
PixelMatrix<T> VoxelMatrix<T>::sum() const
{
  if ( !getSize1() || !getSize2() || !getSize3() )
  {
    SizeError sizeError;
    sizeError.setWhere( "PixelMatrix<T> VoxelMatrix<T>::sum() const" );
    sizeError.setWhat( "This matrix is empty" );
    throw sizeError;
  }

  PixelMatrix<T> m( getSize1(), getSize2() );
  int i, j, k;
  T s;

  for (i = 0; i < getSize1(); ++i)
  for (j = 0; j < getSize2(); ++j)
  {
    s = _planes[0][i][j];
    for (k = 1; k < getSize3(); ++k)
    {
      s += _planes[k][i][j];
    }
    m(i,j) = s;
  }

  return m;
}

/*! Détermine et retourne une matrice de pixels contenant, en chaque
 * position \c i, \c j, la plus petite des valeurs des voxels au travers
 * de l'empilement des plans de cette matrice.
 *
 * \exception SizeError Cette matrice est de taille nulle
 *
 * \sa max, range
****************************************************************/
template<class T>
PixelMatrix<T> VoxelMatrix<T>::min() const
{
  if ( !getSize1() || !getSize2() || !getSize3() )
  {
    SizeError sizeError;
    sizeError.setWhere( "PixelMatrix<T> VoxelMatrix<T>::min() const" );
    sizeError.setWhat( "This matrix is empty" );
    throw sizeError;
  }

  PixelMatrix<T> m( getSize1(), getSize2() );
  int i, j, k;
  T minValue;

  for (i = 0; i < getSize1(); ++i)
  for (j = 0; j < getSize2(); ++j)
  {
    minValue = _planes[0][i][j];
    for (k = 1; k < getSize3(); ++k)
    {
      if ( _planes[k][i][j] < minValue )
      {
        minValue = _planes[k][i][j];
      }
    }
    m(i,j) = minValue;
  }

  return m;
}

/*! Obtain the Z projection of the maximum intensities on the voxelmatrix
****************************************************************/
template<class T>
PixelMatrix<T> VoxelMatrix<T>::getZProjection()
{
//  if ( !getSize1() || !getSize2() || !getSize3() )
//  {
//    SizeError sizeError;
//    sizeError.setWhere( "PixelMatrix<T> VoxelMatrix<T>::min() const" );
//    sizeError.setWhat( "This matrix is empty" );
//    throw sizeError;
//  }

  PixelMatrix<T> zProj( getSize1(), getSize2() );
  zProj.setZeros();
  int i, j, k;

  for (i = 0; i < getSize1(); ++i)
  for (j = 0; j < getSize2(); ++j)
  {
    for (k = 1; k < getSize3(); ++k)
    {
      if ( _planes[k][i][j] > zProj(j,i) )
      {
        zProj(j,i) = _planes[k][i][j];
      }
    }
  }

  return zProj;
}

/*! Détermine et retourne une matrice de pixels contenant, en chaque
 * position \c i, \c j, la plus grande des valeurs des voxels au travers
 * de l'empilement des plans de cette matrice.
 *
 * \exception SizeError Cette matrice est de taille nulle
 *
 * \sa min, range
****************************************************************/
template<class T>
PixelMatrix<T> VoxelMatrix<T>::max() const
{
  if ( !getSize1() || !getSize2() || !getSize3() )
  {
    SizeError sizeError;
    sizeError.setWhere( "PixelMatrix<T> VoxelMatrix<T>::max() const" );
    sizeError.setWhat( "This matrix is empty" );
    throw sizeError;
  }

  PixelMatrix<T> m( getSize1(), getSize2() );
  int i, j, k;
  T maxValue;

  for (i = 0; i < getSize1(); ++i)
  for (j = 0; j < getSize2(); ++j)
  {
    maxValue = _planes[0][i][j];
    for (k = 1; k < getSize3(); ++k)
    {
      if ( _planes[k][i][j] > maxValue )
      {
        maxValue = _planes[k][i][j];
      }
    }
    m(i,j) = maxValue;
  }

  return m;
}

/*! Retourne un vecteur contenant les deux bornes du domaine dans
 * lequel sont comprises les valeurs des voxels de cette matrice.
 *
 * Le premier élément du vecteur retourné contient la borne inférieure.
 * La borne supérieure est donnée par le second élément.
 *
 * Le code
 * \code
 * r = voxelMatrix.range();
 * mini = r[0];
 * maxi = r[1];
 * \endcode
 * est équivalent à
 * \code
 * mini = voxelMatrix.min().min().min();
 * maxi = voxelMatrix.max().max().max();
 * \endcode
 * La première version devrait cependant être généralement plus rapide
 * que la seconde.
 *
 * \exception SizeError Cette matrice est de dimension nulle.
 *
 * \sa min, max
****************************************************************/
template<class T>
Vector<T> VoxelMatrix<T>::range() const
{
  if ( getSize().prod() == 0 )
  {
    SizeError sizeError;
    sizeError.setWhere( "Vector<T> VoxelMatrix<T>::range() const" );
    sizeError.setWhat( "This matrix is empty." );
    throw sizeError;
  }

  Vector<T> rtmp, r( _planes[0].range() );

  for (int k = 1; k < getSize3(); ++k)
  {
    rtmp = _planes[k].range();

    if ( rtmp[1] > r[1] )
    {
      r[1] = rtmp[1];
    }
    else if ( rtmp[0] < r[0] )
    {
      r[0] = rtmp[0];
    }
  }

  return r;
}

/*! Affiche cette matrice de voxels sur \c cout, éventuellement
 * précédée du message \c s s'il est fourni.
****************************************************************/
template<class T>
void VoxelMatrix<T>::print(const string& s) const
{
  if ( !s.empty() )
  {
    cout << s << endl;
  }

  for (int k = 0; k < _h; ++k)
  {
    cout << "-- Plane " << k << endl;
    _planes[k].print();
  }
}

/*! Enregistre cette matrice de voxels dans le fichier \c fileName.
 *
 * Si le fichier existe déjà, il est écrasé si \c overwrite est
 * égal à \c true (\c false par défaut).
 *
 * \exception FileWriteError Le fichier existe déjà et \c overwrite
 * vaut \c false.
 *
 * \sa load
****************************************************************/
template<class T>
void VoxelMatrix<T>::save(const string& fileName, const bool overwrite) const
{
  if ( !overwrite )
  {
    ifstream ifs( fileName.c_str() );
    if ( !ifs.fail() )
    {
      ifs.close();
      FileWriteError fileWriteError;
      fileWriteError.setWhere( "void VoxelMatrix<T>::save(...) const" );
      fileWriteError.setWhat( "File already exists" );
      fileWriteError.setFile( fileName );
      throw fileWriteError;
    }
  }

  ofstream ofs( fileName.c_str(), ios::out | ios::binary );
  if ( ofs.fail() )
  {
    ofs.close();
    FileWriteError fileWriteError;
    fileWriteError.setWhere( "void VoxelMatrix<T>::save(...) const" );
    fileWriteError.setWhat( "Cannot open file for writing" );
    fileWriteError.setFile( fileName );
    throw fileWriteError;
  }

  const int format = 1;
  switch ( format )
  {
    case 0: saveUnderFormat0( ofs ); break;
    case 1: saveUnderFormat1( ofs ); break;
  }

  ofs.close();
}

template<class T>
void VoxelMatrix<T>::saveUnderFormat0(ofstream& ofs) const
{
  const unsigned int size1 = getSize1();
  const unsigned int size2 = getSize2();
  const unsigned int size3 = getSize3();

  ofs.write( (char*)&size1, sizeof(unsigned int) );
  ofs.write( (char*)&size2, sizeof(unsigned int) );
  ofs.write( (char*)&size3, sizeof(unsigned int) );

  for (unsigned int k = 0; k < size3; ++k)
    for (unsigned int i = 0; i < size1; ++i)
    {
      ofs.write( (char*)_planes[k][i].data(), size2*sizeof(T) );
    }
}

template<class T>
void VoxelMatrix<T>::saveUnderFormat1(ofstream& ofs) const
{
  const int zero = 0;
  const int format = 1;
  const int dataType = dataTypeId();
  const int size1 = getSize1();
  const int size2 = getSize2();
  const int size3 = getSize3();
  const int unit = _voxelCalibration.getLengthUnit().getScale();
  const float xsize = _voxelCalibration.getVoxelSize()[0];
  const float ysize = _voxelCalibration.getVoxelSize()[1];
  const float zsize = _voxelCalibration.getVoxelSize()[2];

  ofs.write( (char*)&zero, sizeof(int) );
  ofs.write( (char*)&zero, sizeof(int) );
  ofs.write( (char*)&zero, sizeof(int) );
  ofs.write( (char*)&format, sizeof(int) );
  ofs.write( (char*)&dataType, sizeof(int) );
  ofs.write( (char*)&size1, sizeof(int) );
  ofs.write( (char*)&size2, sizeof(int) );
  ofs.write( (char*)&size3, sizeof(int) );
  ofs.write( (char*)&unit, sizeof(int) );
  ofs.write( (char*)&xsize, sizeof(float) );
  ofs.write( (char*)&ysize, sizeof(float) );
  ofs.write( (char*)&zsize, sizeof(float) );

  for (int k = 0; k < size3; ++k)
    for (int i = 0; i < size1; ++i)
      ofs.write( (char*)_planes[k][i].data(), size2*sizeof(T) );
}

/*! Initialise (ou réinitialise) cette matrice de voxels à partir
 * des données stockées dans le fichier nommé \c fileName.
 *
 * \exception FileReadError Le fichier n'est pas accessible en
 * lecture.
 *
 * \sa save
****************************************************************/
template<class T>
void VoxelMatrix<T>::load(const string& fileName)
{
  ENTER( "void VoxelMatrix<T>::load(const string&)" );

  ifstream ifs( fileName.c_str(), ios::in | ios::binary );
  if ( ifs.fail() )
  {
    FileReadError fileReadError;
    fileReadError.setWhere( "void VoxelMatrix<T>::load(const string&)" );
    fileReadError.setWhat( "Cannot open file for reading" );
    fileReadError.setFile( fileName );
    throw fileReadError;
  }

  unsigned int size1;
  unsigned int size2;
  unsigned int size3;

  ifs.read( (char*)&size1, sizeof(unsigned int) );
  ifs.read( (char*)&size2, sizeof(unsigned int) );
  ifs.read( (char*)&size3, sizeof(unsigned int) );

  if ( size1 == 0 && size2 == 0 && size3 == 0 ) // format >= 1
  {
    load( ifs );
  }
  else // format == 0
  {
    PRINT("Old format");
    EVAL( size1 );
    EVAL( size2 );
    EVAL( size3 );
    setSize( size1, size2, size3 );

    for (unsigned int k = 0; k < size3; ++k)
      for (unsigned int i = 0; i < size1; ++i)
        ifs.read( (char*)_planes[k][i].data(), size2*sizeof(T) );
  }

  ifs.close();

  LEAVE();
}

template<class T>
void VoxelMatrix<T>::loadHeader(const string& filename)
{
  ifstream ifs( filename.c_str(), ios::in | ios::binary );
  if ( ifs.fail() )
  {
    FileReadError fileReadError;
    fileReadError.setWhere( "void VoxelMatrix<T>::loadHeader(const string&)" );
    fileReadError.setWhat( "Cannot open file for reading" );
    fileReadError.setFile( filename );
    throw fileReadError;
  }

  unsigned int size1;
  unsigned int size2;
  unsigned int size3;

  ifs.read( (char*)&size1, sizeof(unsigned int) );
  ifs.read( (char*)&size2, sizeof(unsigned int) );
  ifs.read( (char*)&size3, sizeof(unsigned int) );

  if ( size1 == 0 && size2 == 0 && size3 == 0 ) // format >= 1
  {
    int format;
    int dataType;
    int size1, size2, size3;
    int unitScale;
    Vector<float> voxelSize( 3 );

    ifs.read( (char*)&format, sizeof(int) );
    ifs.read( (char*)&dataType, sizeof(int) );
    ifs.read( (char*)&size1, sizeof(int) );
    ifs.read( (char*)&size2, sizeof(int) );
    ifs.read( (char*)&size3, sizeof(int) );
    ifs.read( (char*)&unitScale, sizeof(int) );
    ifs.read( (char*)&voxelSize[0], sizeof(float) );
    ifs.read( (char*)&voxelSize[1], sizeof(float) );
    ifs.read( (char*)&voxelSize[2], sizeof(float) );

    _n = size1;
    _p = size2;
    _h = size3;

    if ( dataType != dataTypeId() )
    {
      TypeError typeError;
      typeError.setWhere( "void VoxelMatrix<T>::load(ifstream&)" );
      typeError.setWhat( "Type mismatch between template instantiation and file contents" );
      throw typeError;
    }

    LengthUnit lengthUnit;
    lengthUnit.setScale( unitScale );
    _voxelCalibration.setLengthUnit( lengthUnit );
    _voxelCalibration.setVoxelSize( voxelSize );
  }
  else
  {
    _n = size1;
    _p = size2;
    _h = size3;
  }

  ifs.close();
}

template<class T>
void VoxelMatrix<T>::load(ifstream& ifs)
{
  ENTER( "void VoxelMatrix<T>::load(ifstream&)" );

  int format;
  int dataType;
  int size1, size2, size3;
  int unitScale;
  Vector<float> voxelSize( 3 );

  ifs.read( (char*)&format, sizeof(int) );
  ifs.read( (char*)&dataType, sizeof(int) );
  ifs.read( (char*)&size1, sizeof(int) );
  ifs.read( (char*)&size2, sizeof(int) );
  ifs.read( (char*)&size3, sizeof(int) );
  ifs.read( (char*)&unitScale, sizeof(int) );
  ifs.read( (char*)&voxelSize[0], sizeof(float) );
  ifs.read( (char*)&voxelSize[1], sizeof(float) );
  ifs.read( (char*)&voxelSize[2], sizeof(float) );

  if ( dataType != dataTypeId() )
  {
    TypeError typeError;
    typeError.setWhere( "void VoxelMatrix<T>::load(ifstream&)" );
    typeError.setWhat( "Type mismatch between template instantiation and file contents" );
    throw typeError;
  }

  EVAL( format );
  EVAL( dataType );
  EVAL( size1 );
  EVAL( size2 );
  EVAL( size3 );
  EVAL( unitScale );
  EVAL( voxelSize[0] );
  EVAL( voxelSize[1] );
  EVAL( voxelSize[2] );

  setSize( size1, size2, size3 );

  LengthUnit lengthUnit;
  lengthUnit.setScale( unitScale );
  _voxelCalibration.setLengthUnit( lengthUnit );
  _voxelCalibration.setVoxelSize( voxelSize );

  for (int k = 0; k < size3; ++k)
    for (int i = 0; i < size1; ++i)
      ifs.read( (char*)_planes[k][i].data(), size2*sizeof(T) );

  LEAVE();
}

/*! Reslice voxelmatrix changing axis Y for depth
****************************************************************/
template<class T>
VoxelMatrix<T> VoxelMatrix<T>::resliceYZ()
{
  const int size1 = getSize1();
  const int size2 = getSize2();
  const int size3 = getSize3();
  int i, j, k;
  VoxelMatrix<T> temp(size1,size3,size2);

  for (k = 0; k < size3; ++k)
    for (i = 0; i < size1; ++i)
      for (j = 0; j < size2; ++j)
        temp[j][i][k] =_planes[k][i][j];

  temp.setVoxelCalibration(_voxelCalibration);

  return temp;
}

/*! Reslice voxelmatrix changing axis X for depth
****************************************************************/
template<class T>
VoxelMatrix<T> VoxelMatrix<T>::resliceXZ()
{
  const int size1 = getSize1();
  const int size2 = getSize2();
  const int size3 = getSize3();
  int i, j, k;
  VoxelMatrix<T> temp(size3,size2,size1);

  for (k = 0; k < size3; ++k)
    for (i = 0; i < size1; ++i)
      for (j = 0; j < size2; ++j)
        temp[i][k][j] =_planes[k][i][j];

  temp.setVoxelCalibration(_voxelCalibration);

  return temp;
}

/*! Rotate voxelmatrix changing axis X and Y
****************************************************************/
template<class T>
VoxelMatrix<T> VoxelMatrix<T>::rotateXY()
{
  const int size1 = getSize1();
  const int size2 = getSize2();
  const int size3 = getSize3();
  int i, j, k;
  VoxelMatrix<T> temp(size2,size1,size3);

  for (k = 0; k < size3; ++k)
  {
    for (i = 0; i < size1; ++i)
    {
      for (j = 0; j < size2; ++j)
      {
        temp[k][j][i] =_planes[k][i][j];
      }
    }
  }

  temp.setVoxelCalibration(_voxelCalibration);

  return temp;
}

/*! Fill holes line by line
****************************************************************/
template<class T>
void VoxelMatrix<T>::fillIt(VoxelMatrix<float>& voxelMatrix)
{
  const int size1 = voxelMatrix.getSize1();
  const int size2 = voxelMatrix.getSize2();
  const int size3 = voxelMatrix.getSize3();

  HolesFillingf holesFilling;
  for (int k = 0; k < size3; ++k) holesFilling.apply( voxelMatrix[k] );
  //voxelMatrix.fillLineGaps(voxelMatrix);

  VoxelMatrix<float> tempYZ = voxelMatrix.resliceYZ();
  for (int j = 0; j < size2; ++j)  holesFilling.apply( tempYZ[j] );
  //tempYZ.fillLineGaps(tempYZ);
  voxelMatrix = tempYZ.resliceYZ();
  VoxelMatrix<float> tempXZ = voxelMatrix.resliceXZ();
  //tempXZ.fillLineGaps(tempXZ);
  for (int i = 0; i < size1; ++i)  holesFilling.apply( tempXZ[i] );
  voxelMatrix = tempXZ.resliceXZ();

  for (int k = 0; k < size3; ++k) holesFilling.apply( voxelMatrix[k] );

  holesFilling.apply( voxelMatrix );

}

/*! Fill image in Y direction
****************************************************************/
template<class T>
void VoxelMatrix<T>::fillLineGaps(VoxelMatrix<T>& voxelMatrix)
{
  const int size1 = getSize1();
  const int size2 = getSize2();
  const int size3 = getSize3();
  int i, j, k;
  Vector<int> startPos(3), endPos(3);
  bool data = false, finished = false, toFill = false;

  for (k = 0; k < size3; ++k)
  {
    data = false; finished = false; toFill = false;
    startPos.setZeros();
    endPos.setZeros();
    for (i = 0; i < size1; ++i)
    {
      data = false; finished = false; toFill = false;
      startPos.setZeros();
      endPos.setZeros();
      for (j = 0; j < size2; ++j)
      {
        if ( ( voxelMatrix._planes[k](i,j) == 1 ) && ( data == false ) && ( finished == false ) )
        {
          data = true;
        }
        else if ( ( voxelMatrix._planes[k](i,j) == 0 ) && ( data == true ) && ( finished == false ) )
        {
          finished = true;
          startPos[0] = i;
          startPos[1] = j;
          startPos[2] = k;
        }
        else if ( ( voxelMatrix._planes[k](i,j) == 1 ) && ( data == true ) && ( finished == true ) )
        {
          toFill = true;
          endPos[0] = i;
          endPos[1] = j;
          endPos[2] = k;
        }
        if ( ( toFill == true ) )//&& ( (size2/(endPos[1]-startPos[1])) > 3.5 ) )
        {
          for (j = startPos[1]; j < endPos[1]; ++j)
            voxelMatrix._planes[k](i,j) = 1;
          data = false; finished = false; toFill = false;
          startPos.setZeros();
          endPos.setZeros();
        }
      }
    }
  }
}

/*! Performs type-conversion: returns this matrix as a matrix of
 * values of type \c NewType.
 *
 * Conversions are only instantiated at present between \c int,
 * \c float, \c double, and \c long \c double types.
****************************************************************/
template<class T>
template<typename NewType>
VoxelMatrix<NewType> VoxelMatrix<T>::toType() const
{
  const int size1 = getSize1();
  const int size2 = getSize2();
  const int size3 = getSize3();
  VoxelMatrix<NewType> voxelMatrix;
  voxelMatrix.setSize( getSize() );
  int i, j, k;

  for (i = 0; i < size1; ++i)
    for (j = 0; j < size2; ++j)
      for (k = 0; k < size3; ++k)
        voxelMatrix(i,j,k) = static_cast<NewType>( operator()(i,j,k) );

  return voxelMatrix;
}

template<class T>
int VoxelMatrix<T>::dataTypeId() const
{
  T _t;
  int _int;
  unsigned int _unsignedInt;
  float _float;
  double _double;
  long double _longDouble;

  if ( typeid(_t) == typeid(_int) ) return 2;
  if ( typeid(_t) == typeid(_unsignedInt) ) return 3;
  if ( typeid(_t) == typeid(_float) ) return 5;
  if ( typeid(_t) == typeid(_double) ) return 6;
  if ( typeid(_t) == typeid(_longDouble) ) return 7;

  return 0;
}

template class VoxelMatrix<char>;
template class VoxelMatrix<unsigned char>;
template class VoxelMatrix<short>;
template class VoxelMatrix<unsigned short>;
template class VoxelMatrix<int>;
template class VoxelMatrix<unsigned int>;
template class VoxelMatrix<float>;
template class VoxelMatrix<double>;
template class VoxelMatrix<long double>;

template VoxelMatrix<int> VoxelMatrix<unsigned char>::toType() const;
template VoxelMatrix<float> VoxelMatrix<unsigned char>::toType() const;
template VoxelMatrix<double> VoxelMatrix<unsigned char>::toType() const;
template VoxelMatrix<long double> VoxelMatrix<unsigned char>::toType() const;

template VoxelMatrix<unsigned char> VoxelMatrix<int>::toType() const;
template VoxelMatrix<float> VoxelMatrix<int>::toType() const;
template VoxelMatrix<double> VoxelMatrix<int>::toType() const;
template VoxelMatrix<long double> VoxelMatrix<int>::toType() const;

template VoxelMatrix<unsigned char> VoxelMatrix<float>::toType() const;
template VoxelMatrix<int> VoxelMatrix<float>::toType() const;
template VoxelMatrix<double> VoxelMatrix<float>::toType() const;
template VoxelMatrix<long double> VoxelMatrix<float>::toType() const;

template VoxelMatrix<unsigned char> VoxelMatrix<double>::toType() const;
template VoxelMatrix<int> VoxelMatrix<double>::toType() const;
template VoxelMatrix<float> VoxelMatrix<double>::toType() const;
template VoxelMatrix<long double> VoxelMatrix<double>::toType() const;

template VoxelMatrix<unsigned char> VoxelMatrix<long double>::toType() const;
template VoxelMatrix<int> VoxelMatrix<long double>::toType() const;
template VoxelMatrix<float> VoxelMatrix<long double>::toType() const;
template VoxelMatrix<double> VoxelMatrix<long double>::toType() const;
