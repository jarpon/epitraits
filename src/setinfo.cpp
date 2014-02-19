#include <iostream>
#include <dataset.h>
#include <voxelmatrix.h>

void setData(DataSet& dataset, const string& infoColumn, int& numImage, float value )
{
  dataset.setValue( infoColumn, numImage, value );
  return;
}

void setData( DataSet& dataset, const string& infoColumn, int& numImage, const string value )
{
  dataset.setValue( infoColumn, numImage, value );
  return;
}
