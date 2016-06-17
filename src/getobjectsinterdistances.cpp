//#include <componentlabelling.h>
#include <dataset.h>
//#include "regionanalysis2.h"
//#include "regionanalysis.h"
#include <cmath>
//#include <marchingcubes.h>
//#include <thresholding.h>
#include <trimesh.h>
#include <sstream>

#define TRACE
#include <trace.h>

void chromocentersInterdistances(const string& filename,
                           const string& parentDir,
                                 int& totalNumCCs,
                                 DataSet& completeDataset,
                                 DataSet& individualChromocentersDataset)
{
  //const DataSet datasetNucleus( analysisDir + filename + "_chromocenters.csv" );

  const DataSet ccsInfo( parentDir + "/analysis/ccs.data" );

  Vector<string> tempFileNames;
  tempFileNames = ccsInfo.getValues<string>( ccsInfo.variableNames()[0] );

  int lastPos, numCCS = 0;

  for ( int j = 0; j < tempFileNames.getSize(); ++j )
    if ( tempFileNames[j] == filename )
    {
      lastPos = j;
      ++ numCCS;
    }

  if ( numCCS == 0 )
  {
    EVAL("Nucleus not found");
    return;
  }


  Vertices<float> vertices( 3, numCCS, 0, 0 );
  int k = 0;
  for ( int j = lastPos - numCCS + 1 ; j < lastPos + 1; ++j, ++k )
  {
    vertices[k][0] = ccsInfo.getValue<float>( "centroidCoordX", j );
    vertices[k][1] = ccsInfo.getValue<float>( "centroidCoordY", j );
    vertices[k][2] = ccsInfo.getValue<float>( "centroidCoordZ", j );
    EVAL(vertices[k]);
  }


  //const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/nuclei/" + filename + ".tm" );
  const TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/territories/" + filename + "-01.tm" );
  float distanceToBorder;
  Vector<float> vertexTriMesh(3);

  float tempDistance;
  Vector<float> temp1, temp2;

  for (int i = 0; i < numCCS; ++i)
  {
    temp1 = vertices[i];
    nucleusTriMesh.closestPoint( temp1, vertexTriMesh );
    distanceToBorder = temp1.distance( vertexTriMesh );

    for (int j = 0; j < numCCS; ++j)
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
    individualChromocentersDataset.setValue ( "distanceToTheBorder", i, distanceToBorder );
    completeDataset.setValue ( "name", i+totalNumCCs, filename );
    completeDataset.setValue ( "distanceToTheBorder", i+totalNumCCs, distanceToBorder );

  }
  totalNumCCs += numCCS;
}
