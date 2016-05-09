#ifndef CONVERGENCETEST2_H
#define CONVERGENCETEST2_H

#include <alinevector.h>

template<class T>
class ConvergenceTest2
{
  public:

    ConvergenceTest2();

    bool hasData() const;
    void setData(const Vector<T>&);
    const Vector<T>& getData() const;

    void setRange(const int);
    int getRange() const;

    bool isPositive(const int) const;
    bool isNegative(const int) const;

  private:

    const Vector<T>* _data;
    int _range;
};

#endif // CONVERGENCETEST2_H
