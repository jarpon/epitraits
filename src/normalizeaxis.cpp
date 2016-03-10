#include <trimesh.h>

void normalizeAxisTriMesh( const string& filename, const string& parentDir )
{
  TriMesh<float> shapeCopy;
  shapeCopy.load( parentDir + "/shapes/nuclei/" + filename + ".tm" );
  shapeCopy.applyPAT();

  shapeCopy.save( parentDir + "/shapes/nuclei/" + filename + "-normalized.tm" );
}
