#ifndef SPATIALMODELEVALUATOR_H
#define SPATIALMODELEVALUATOR_H

#include <spatialdescriptor.h>
#include <spatialmodel.h>

#include <dataset.h>

template<class CoordType,class PixelType>
class SpatialModelEvaluator
{
  public:

    SpatialModelEvaluator();

    void setModel(SpatialModel<CoordType,PixelType>&);
    void setDescriptor(SpatialDescriptor<CoordType>&);
    void addDescriptor(SpatialDescriptor<CoordType>&);
    void setNumMonteCarloSamples(const int);
    void setPrecision(const CoordType&);

    float eval(const Vertices<CoordType>&,DataSet* = 0);
    Vector<float> evalSDIandMaxDiff(const Vertices<CoordType>&,DataSet* = 0);
    void eval(const Vertices<CoordType>&,vector<float>&,vector<int>&,DataSet* =0);
    void evalSDIandMaxDiff(const Vertices<CoordType>&,vector<float>&,vector<int>&,vector<float>&,DataSet* =0);
    //void evalSDIandMaxDiff(const Vertices<CoordType>&,vector<float>&,vector<int>&,vector<float>&,DataSet* =0);

  private:

    void eval(
      const Vertices<CoordType>&,
      SpatialDescriptor<CoordType>&,
      const ShapeSet<CoordType>&,
      const ShapeSet<CoordType>&,
      float&,int&,DataSet* =0);

    void eval(
      const Vertices<CoordType>&,
      SpatialDescriptor<CoordType>&,
      const ShapeSet<CoordType>&,
      const ShapeSet<CoordType>&,
      float&,int&,float&,DataSet* =0);

    SpatialModel<CoordType,PixelType>* _model;
    vector<SpatialDescriptor<CoordType>*> _descriptors;
    int _numMonteCarloSamples;
    CoordType _precision;
};

#endif // SPATIALMODELEVALUATOR_H
