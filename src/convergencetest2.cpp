
// 2015.10.13

#include "convergencetest2.h"

#include <cmath>

template<class T>
ConvergenceTest2<T>::ConvergenceTest2()
{
  //_range = 30;
  _range = 1000;
  _data = 0;
}

template<class T>
bool ConvergenceTest2<T>::hasData() const
{
  return _data != 0;
}

template<class T>
void ConvergenceTest2<T>::setData(const Vector<T>& data)
{
  _data = &data;
}

template<class T>
const Vector<T>& ConvergenceTest2<T>::getData() const
{
  return *_data;
}

template<class T>
void ConvergenceTest2<T>::setRange(const int range)
{
  _range = range;
}

template<class T>
int ConvergenceTest2<T>::getRange() const
{
  return _range;
}

template<class T>
bool ConvergenceTest2<T>::isPositive(const int atStep) const
{
  const Vector<T>& data = getData();
  bool result = false;

  if ( atStep >= 2*_range )
  {
    Vector<T> sub2 = data.copy( atStep-2*_range, atStep-_range-1 );
    Vector<T> sub1 = data.copy( atStep-_range, atStep-1 );
    sub1 -= sub2;
    const float z = fabs( sub1.mean() ) / sub1.sem( false );
    //z /= sqrt( sub1.var()/sub1.getSize() );
    //result = z < 1.96;
    result = z < 1.56;
  }

  return result;
}

template<class T>
bool ConvergenceTest2<T>::isNegative(const int atStep) const
{
  return !isPositive( atStep );
}

template class ConvergenceTest2<float>;
template class ConvergenceTest2<double>;
template class ConvergenceTest2<long double>;

