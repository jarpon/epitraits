#include "cdftools2.h"
#include <dataset.h>

#define TRACE
#include <trace.h>

using namespace std;

void test2Distributions(const string&nameFile1, const string& descriptorName, const string& nameFile2 )
{
  EVAL( "test2Distributions..." );
  EVAL( nameFile1 );
  const DataSet file1Dataset( nameFile1 );

//  Vector<float> descriptorSDI1 = file1Dataset.getValues<float>( descriptorName );
//  Vector<float> descriptorSDI2 = file1Dataset.getValues<float>( descriptorName );

//  descriptorSDI1.sort();
//  descriptorSDI2.sort();


//  const DataSet file1Dataset( "/home/jarpon/data/projects/testStatisticalTests/newOne_15july/toTest.csv" );
//  const DataSet file1Dataset( "/home/jarpon/data/projects/testStatisticalTests/newOne_15july/toTest40.csv" );



  Vector<float> x1 = file1Dataset.getValues<float>( "F1-SDI" );
  Vector<float> x2 = file1Dataset.getValues<float>( "F2-SDI" );

  x1.sort();
  x2.sort();

  Vector<float> y1, y2;

  CDFTools<float> cdfTools;

  y1 = cdfTools.cdf( x1 );
  y2 = cdfTools.cdf( x2 );

  Vector<float> differenceBetweenCurves;

  differenceBetweenCurves = cdfTools.areasDifference( y1, x1, y2, x2 );
  EVAL(differenceBetweenCurves);

  //repetition
  Vector<float> x1G = file1Dataset.getValues<float>( "G1-SDI" );
  Vector<float> x2G = file1Dataset.getValues<float>( "G2-SDI" );

  x1G.sort();
  x2G.sort();

  y1.setZeros();
  y2.setZeros();
  y1 = cdfTools.cdf( x1G );
  y2 = cdfTools.cdf( x2G );

  differenceBetweenCurves = cdfTools.areasDifference( y1, x1G, y2, x2G );

  EVAL(differenceBetweenCurves);


  Vector<float> x1H = file1Dataset.getValues<float>( "H1-SDI" );
  Vector<float> x2H = file1Dataset.getValues<float>( "H2-SDI" );

  x1H.sort();
  x2H.sort();

  y1.setZeros();
  y2.setZeros();
  y1 = cdfTools.cdf( x1H );
  y2 = cdfTools.cdf( x2H );

  differenceBetweenCurves = cdfTools.areasDifference( y1, x1H, y2, x2H );

  EVAL(differenceBetweenCurves);
}
