#ifndef VOXELMATRIX_H
#define VOXELMATRIX_H

#include "pixelmatrix.h"
#include "voxelcalibration.h"

template<class T>
class VoxelMatrix
{
  public:

    VoxelMatrix();
    VoxelMatrix(const int,const int,const int);
    VoxelMatrix(const VoxelMatrix&);
    VoxelMatrix(const string&);
    ~VoxelMatrix();

    void setSize(const int,const int,const int);
    void setSize(const Vector<int>&);
    int getSize1() const;
    int getSize2() const;
    int getSize3() const;
    int getSize(const int) const;
    Vector<int> getSize() const;
    bool sameSizeAs(const VoxelMatrix&) const;

    void setVoxelCalibration(const VoxelCalibration&);
    const VoxelCalibration& getVoxelCalibration() const;

    VoxelMatrix& operator=(const VoxelMatrix&);
    void swap(VoxelMatrix&);
    bool operator==(const VoxelMatrix&);
    bool operator!=(const VoxelMatrix&);

    void setZeros();
    void setOnes();
    void fill(const T&);
    void apply(T(*)(T));

    void reverse();
    PixelMatrix<T>& operator[](const int);
    const PixelMatrix<T>& operator[](const int) const;
    T& operator()(const int,const int,const int);
    const T& operator()(const int,const int,const int) const;
    T& operator()(const Vector<int>&);
    const T& operator()(const Vector<int>&) const;
    bool contains(const int,const int,const int) const;
    bool contains(const Vector<int>&) const;

    VoxelMatrix copy(
      const int,const int,
      const int,const int,
      const int,const int) const;
    void paste(
      const int,
      const int,
      const int,
      const VoxelMatrix&);
    Vector<T> cat() const;

    void operator+=(const VoxelMatrix&);
    void operator-=(const VoxelMatrix&);
    void operator*=(const VoxelMatrix&);
    void operator+=(const T&);
    void operator-=(const T&);
    void operator*=(const T&);
    void operator/=(const T&);

    Vector<unsigned int> histogram(
      const T,
      const T,
      const unsigned int) const;
    Vector<unsigned int> histogram(
      const T,
      const T,
      const unsigned int,
      const VoxelMatrix&) const;
    PixelMatrix<T> sum() const;
    PixelMatrix<T> min() const;
    PixelMatrix<T> max() const;
    Vector<T> range() const;

    void print(const string& ="") const;
    void save(const string&,const bool=false) const;
    void load(const string&);
    void loadHeader(const string&);

    template<typename NewType>
    VoxelMatrix<NewType> toType() const;

    PixelMatrix<T> getZProjection();
    VoxelMatrix<T> resliceXZ();
    VoxelMatrix<T> resliceYZ();
    VoxelMatrix<T> rotateXY();
    void fillIt(VoxelMatrix<float>&);
    void fillLineGaps(VoxelMatrix&);

  private:

    int dataTypeId() const;
    void saveUnderFormat0(ofstream&) const;
    void saveUnderFormat1(ofstream&) const;
    void load(ifstream&);

    int _n; //!< Nombre de lignes de chacun des plans
    int _p; //!< Nombre de colonnes de chacun des plans
    int _h; //!< Nombre de plans images empilés
    PixelMatrix<T>* _planes; //!< Empilement des plans images
    VoxelCalibration _voxelCalibration;
};

typedef VoxelMatrix<short> VoxelMatrixs;
typedef VoxelMatrix<unsigned short> VoxelMatrixus;
typedef VoxelMatrix<int> VoxelMatrixi;
typedef VoxelMatrix<unsigned int> VoxelMatrixui;
typedef VoxelMatrix<double> VoxelMatrixd;
typedef VoxelMatrix<long double> VoxelMatrixld;


/*! Retourne la référence du plan d'indice \c k le long de la troisième
 * dimension (Z) de cette matrice.
****************************************************************/
template<class T>
inline PixelMatrix<T>& VoxelMatrix<T>::operator[](const int k)
{
  return _planes[k];
}

/*! Retourne la référence du plan d'indice \c k le long de la troisième
 * dimension (Z) de cette matrice.
****************************************************************/
template<class T>
inline const PixelMatrix<T>& VoxelMatrix<T>::operator[](const int k) const
{
  return _planes[k];
}

/*! Retourne la référence du voxel situé ligne \c i, colonne \c j
 * dans le plan \c k.
****************************************************************/
template<class T>
inline T& VoxelMatrix<T>::operator()(
  const int i,
  const int j,
  const int k)
{
  return _planes[k][i][j];
}

/*! Retourne la référence du voxel situé ligne \c i, colonne \c j
 * dans le plan \c k.
****************************************************************/
template<class T>
inline const T& VoxelMatrix<T>::operator()(
  const int i,
  const int j,
  const int k) const
{
  return _planes[k][i][j];
}

/*! This convenience function returns the reference
 * to the voxel at the specified position.
****************************************************************/
template<class T>
inline T& VoxelMatrix<T>::operator()(const Vector<int>& v)
{
  return _planes[v[2]][v[0]][v[1]];
}

/*! This convenience function returns the reference
 * to the voxel at the specified position.
****************************************************************/
template<class T>
inline const T& VoxelMatrix<T>::operator()(const Vector<int>& v) const
{
  return _planes[v[2]][v[0]][v[1]];
}

/*! Ecrit la matrice \c voxelMatrix sur le flux de sortie \c os.
 *
 * Le format d'affichage peut être contrôlé par un appel préalable
 * aux méthodes precision() et width() sur \c os.
 *
 * \relates VoxelMatrix
****************************************************************/
template<class T>
inline ostream& operator<<(ostream& os, const VoxelMatrix<T>& voxelMatrix)
{
  const streamsize p = os.precision();
  const streamsize w = os.width();

  for (int k = 0; k < voxelMatrix.getSize3(); ++k)
  {
    os.precision( p );
    os.width( w );
    os << voxelMatrix[k];
    if ( k < voxelMatrix.getSize3()-1 )
    {
      os << endl;
    }
  }

  return os;
}

#endif
