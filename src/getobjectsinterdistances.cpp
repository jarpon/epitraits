//#include <componentlabelling.h>
#include <dataset.h>
//#include "regionanalysis2.h"
//#include <regionanalysis.h>
#include <cmath>
//#include <marchingcubes.h>
//#include <thresholding.h>
#include <trimesh.h>
#include <sstream>

#define TRACE
#include <trace.h>

void chromocentersInterdistances(const string& filename,
                           const string& analysisDir, DataSet& individualChromocentersDataset)
{
  const DataSet datasetNucleus( analysisDir + filename + "_chromocenters.csv" );
  const int numCompartments = datasetNucleus.size()[0];

  Vertices<float> vertices ( 3, numCompartments, 0, 0 );
  for ( int i = 0; i < numCompartments; ++i )
  {
    vertices[i][0] = datasetNucleus.getValue<float>( "centroidCoordX", i );
    vertices[i][1] = datasetNucleus.getValue<float>( "centroidCoordY", i );
    vertices[i][2] = datasetNucleus.getValue<float>( "centroidCoordZ", i );
    EVAL(vertices[i]);
  }

  float tempDistance;
  Vector<float> temp1, temp2;

  for (int i = 0; i < numCompartments; ++i)
  {
    temp1 = vertices[i];

    for (int j = 0; j < numCompartments; ++j)
    {
      if ( i == j )
        tempDistance = 0;
      else
      {
        temp2 = vertices[j];
        tempDistance = temp1.distance(temp2);
      }
      ostringstream iss; //we suppose as much 99 labels
      iss << (j+1);
      individualChromocentersDataset.setValue ( iss.str(), i, tempDistance );
    }
  }
}
