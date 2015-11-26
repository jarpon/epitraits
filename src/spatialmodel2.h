#ifndef SPATIALMODEL2_H
#define SPATIALMODEL2_H

#include <pixelmatrix.h>
#include <pixelsampler.h>

#include <shapeset.h>

template<class CoordType,class PixelType>
class SpatialModel2
{
  public:

    SpatialModel2();
    virtual ~SpatialModel2();

    virtual void setRandomGenerator(RandomGenerator&);
    void setMask(const PixelMatrix<PixelType>&);
    void setPixelSize(const Vector<CoordType>&);
    RandomGenerator& getRandomGenerator();
    const PixelMatrix<PixelType>& getMask() const;
    const Vector<CoordType>& getPixelSize() const;
    virtual void initialize() = 0;

//    virtual Vertices<CoordType> drawSample(const int) = 0;
//    virtual Vertices<CoordType> drawSample(const int, const int) = 0;

//    ShapeSet<CoordType> drawSamples(const int,const int);
//    ShapeSet<CoordType> drawSamples(const int,const int, const int);

  private:

    RandomGenerator* _randomGenerator;
    const PixelMatrix<PixelType>* _mask;
    Vector<CoordType> _pixelSize;
};

#endif // SPATIALMODEL2_H
