#ifndef CDFTOOLS2_H
#define CDFTOOLS2_H

#include <alinematrix.h>

#include <vector>

template<typename T>
class CDFTools
{
  public:

    CDFTools();
    ~CDFTools();

    Vector<T> cdf(const Vector<T>&) const;
    Vector<T> cdf(const Vector<T>&, const Vector<T>&) const;

    Vector<T> average(const vector<const Vector<T>*>&) const;
    Vector<T> average(const Matrix<T>&) const;
    Vector<T> average(const vector<const Vector<T>*>&, const Vector<T>&) const;
    Vector<T> average(const Matrix<T>&, const Vector<T>&) const;

    Vector<T> percentile(const vector<const Vector<T>*>&, const float) const;
    Vector<T> percentile(const Matrix<T>&, const float) const;
    Vector<T> percentile(const vector<const Vector<T>*>&, const Vector<T>&, const float) const;
    Vector<T> percentile(const Matrix<T>&, const Vector<T>&, const float) const;

    Matrix<T> percentiles(const vector<const Vector<T>*>&, const Vector<T>&, const Vectorf&) const;
    Matrix<T> percentiles(const Matrix<T>&, const Vector<T>&, const Vectorf&) const;

    void differences(
      const Vector<T>&,
      const Vector<T>&,
      const Vector<T>&,
      const Vector<T>&,
      Vector<T>&,
      Vector<T>&,
      Vector<T>&) const;
    void differences(
      const Vector<T>&,
      const Vector<T>&,
      Vector<T>&,
      Vector<T>&,
      Vector<T>&) const;
    void differences(
      const Vector<T>&,
      const Vector<T>&,
      Vector<T>&) const;

    Vector<T> areasDifference(
      const Vector<T>&,
      const Vector<T>&,
      const Vector<T>&,
      const Vector<T>&) const;

    int rank(const Vector<T>&, const Matrix<T>&, const Vector<T>&, const Vector<T>&) const;
    Vector<float> rankAndMaxDiff(const Vector<T>&, const Matrix<T>&, const Vector<T>&, const Vector<T>&) const;

};

#endif
