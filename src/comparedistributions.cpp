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

  Vector<float> descriptorSDI1 = file1Dataset.getValues<float>( "F1-SDI" );
  Vector<float> descriptorSDI2 = file1Dataset.getValues<float>( "F2-SDI" );

  descriptorSDI1.sort();
  descriptorSDI2.sort();

  Vector<float> y1, y2;
//  y1.setSize( ecdf1.getSize() );
//  y2.setSize( ecdf2.getSize() );

  CDFTools<float> cdfTools;

  y1 = cdfTools.cdf( descriptorSDI1 );
  y2 = cdfTools.cdf( descriptorSDI2 );

  Vector<float> differenceBetweenCurves;
  Vector<float> ecdf1;
  ecdf1 = cdfTools.cdf( descriptorSDI1, y1 );
  Vector<float> ecdf2;
  ecdf2 = cdfTools.cdf( descriptorSDI2, y2 );
  EVAL(ecdf1);
  EVAL(ecdf2);
  EVAL(y1);
  EVAL(y2);

  differenceBetweenCurves = cdfTools.areasDifference( descriptorSDI1, descriptorSDI2, y1, y2 );

  EVAL(differenceBetweenCurves);

  //repetition
  Vector<float> descriptorSDI11 = file1Dataset.getValues<float>( "G1-SDI" );
  Vector<float> descriptorSDI22 = file1Dataset.getValues<float>( "G2-SDI" );

  descriptorSDI11.sort();
  descriptorSDI22.sort();

  y1.setZeros();
  y2.setZeros();
  y1 = cdfTools.cdf( descriptorSDI11 );
  y2 = cdfTools.cdf( descriptorSDI22 );

  ecdf1.setZeros();
  ecdf2.setZeros();
  ecdf1 = cdfTools.cdf( descriptorSDI11, y1 );
  ecdf2 = cdfTools.cdf( descriptorSDI22, y2 );
  EVAL(ecdf2);

  differenceBetweenCurves = cdfTools.areasDifference( descriptorSDI11, descriptorSDI22, y1, y2 );

  EVAL(differenceBetweenCurves);


  Vector<float> descriptorSDI1H = file1Dataset.getValues<float>( "H1-SDI" );
  Vector<float> descriptorSDI2H = file1Dataset.getValues<float>( "H2-SDI" );

  descriptorSDI1H.sort();
  descriptorSDI2H.sort();

  y1.setZeros();
  y2.setZeros();
  y1 = cdfTools.cdf( descriptorSDI1H );
  y2 = cdfTools.cdf( descriptorSDI2H );
  ecdf1.setZeros();
  ecdf2.setZeros();
  ecdf1 = cdfTools.cdf( descriptorSDI1H, y1 );
  ecdf2 = cdfTools.cdf( descriptorSDI2H, y2 );
  EVAL(ecdf2);

  differenceBetweenCurves = cdfTools.areasDifference( descriptorSDI1H, descriptorSDI2H, y1, y2 );

  EVAL(differenceBetweenCurves);
}
