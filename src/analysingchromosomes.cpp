#include <componentlabelling.h>
#include <dataset.h>
//#include "regionanalysis2.h"
#include <regionanalysis3d.h>
#include <cmath>
#include <marchingcubes.h>
#include <thresholding.h>
#include <trimesh.h>

#define TRACE
#include <trace.h>

void chromosomesAnalysis(VoxelMatrix<float>& chromosomeMask, const string& filename, const string& parentDir,
                           const int& numNucleus, int& totalNumChromosomes,
                           DataSet& nucleiDataset, DataSet& chromosomesDataset, DataSet& individualChromosomesDataset)
{
  VoxelMatrix<float> originalVoxelMatrix( parentDir + "/originals_vm/" + filename + ".vm" );
  VoxelMatrix<float> nucleusMask( parentDir + "/segmented_nuclei/" + filename + ".vm" );

//  RegionAnalysis3D<float> regionAnalysis;
//  regionAnalysis.setLabelMatrix( nucleusMask );
//  regionAnalysis.run();

////  ComponentLabelling<float> componentLabelling;
////  componentLabelling.apply( chromosomeMask );
////  EVAL( componentLabelling.getNumLabels() );
////  EVAL( chromosomeMask.max().max().max() );

//  RegionAnalysis3D<float> regionAnalysisChromosomes;
//  regionAnalysisChromosomes.setLabelMatrix( chromosomeMask );
//  regionAnalysisChromosomes.run();

  //uncomment



//  //TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + "_nucleus.tm" );
//  TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/nuclei/" + filename + ".tm" );
//  //TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + ".tm" );

//  Vector<float> centroid(3);
//  Vector<float> vertexTriMesh(3);

//  Vertices<float> centroids;
//  centroids = regionAnalysisChromosomes.computeRegionCentroids(REGION_FEATURE_CENTROID, originalVoxelMatrix);
//  nucleiDataset.setValue ( "name", numNucleus, filename );//filename
//  //nucleiDataset.setValue ( "chromosomesNumber", numNucleus, componentLabelling.getNumLabels() );//number of ccs obtained in the nucleus
//  nucleiDataset.setValue ( "chromosomesNumber", numNucleus, chromosomeMask.max().max().max() );//number of ccs obtained in the nucleus
//  nucleiDataset.setValue ( "chromosomesVolume_vm", numNucleus, regionAnalysisChromosomes.computeRegionFeature(REGION_FEATURE_VOLUME,originalVoxelMatrix).sum() );//total ccs volume
//  //nucleiDataset.setValue ( "equivalentRadius_tm", numNucleus, abs(triMesh.equivalentRadius()) );
//  //nucleiDataset.setValue ( "nucleusVolume_tm", numNucleus, abs(triMesh.volume()) );//nucleus volume got from the trimesh
//  //nucleiDataset.setValue ( "surfaceArea_tm", numNucleus, abs(triMesh.area() ) );//nucleus volume got from the trimesh
//  nucleiDataset.setValue ( "chromosomesRelativeVolume", numNucleus, ( regionAnalysisChromosomes.computeRegionFeature(REGION_FEATURE_VOLUME,originalVoxelMatrix)[0] / regionAnalysis.computeRegionFeature(REGION_FEATURE_VOLUME,originalVoxelMatrix)[0] ) );//relative volume of ccs regarding the complete nucleus
//  nucleiDataset.setValue ( "chromosomesIntensity", numNucleus, regionAnalysisChromosomes.computeRegionFeature(REGION_FEATURE_INTENSITY,originalVoxelMatrix).sum() );//total absolute intensity of ccs
//  nucleiDataset.setValue ( "chromosomesIntegratedDensity", numNucleus, regionAnalysisChromosomes.computeRegionFeature(REGION_FEATURE_INTEGRATED_DENSITY,originalVoxelMatrix).sum() );//integrated density of ccs taking into account real volume
//  nucleiDataset.setValue ( "chromosomes-NucleusIntensityRate", numNucleus, regionAnalysisChromosomes.computeRegionFeature(REGION_FEATURE_INTENSITY,originalVoxelMatrix).sum() / regionAnalysis.computeRegionFeature(REGION_FEATURE_INTENSITY,originalVoxelMatrix)[0] );//intensity RHF: rate of heterochromatin (int.Density of ccs/ int.Density of nuclei)
//  nucleiDataset.setValue ( "chromosomes-NucleusVolumeRate", numNucleus, regionAnalysisChromosomes.computeRegionFeature(REGION_FEATURE_VOLUME,originalVoxelMatrix).sum() / regionAnalysis.computeRegionFeature(REGION_FEATURE_VOLUME,originalVoxelMatrix)[0] );//volume RHF: rate of heterochromatin (int.Density of ccs/ int.Density of nuclei)


///**///chromocenters individual information
//  for (int numChromosome = 0; numChromosome < chromosomeMask.max().max().max(); numChromosome++ )
//  {
//    chromosomesDataset.setValue ( "name", numChromosome+totalNumChromosomes, filename );
//    chromosomesDataset.setValue ( "idChromosomes", numChromosome+totalNumChromosomes, numChromosome+1 );

//    centroid = centroids[numChromosome];
//    chromosomesDataset.setValue ( "centroidCoordX", numChromosome+totalNumChromosomes, centroid[X] );
//    chromosomesDataset.setValue ( "centroidCoordY", numChromosome+totalNumChromosomes, centroid[Y] );
//    chromosomesDataset.setValue ( "centroidCoordZ", numChromosome+totalNumChromosomes, centroid[Z] );

//    MarchingCubes<float> marchingCubes;
//    TriMesh<float> triMesh;
//    VoxelMatrix<float> currentLabeledVM = chromosomeMask;
//    Thresholding<float> thresholding;
//    thresholding.setForeground( 1.0 );
//    thresholding.setBackground( 0.0 );
//    thresholding.levelSetMask( currentLabeledVM, numChromosome+1 );
//    triMesh = marchingCubes.buildMesh( currentLabeledVM, 0.5, true );
//    triMesh.scale( originalVoxelMatrix.getVoxelCalibration().getVoxelSize() );

//    nucleusTriMesh.closestPoint( centroid, vertexTriMesh );
//    float distanceToBorder = centroid.distance( vertexTriMesh );
//    float chromosomeVolume_tm = fabs(triMesh.volume());
//    float eqRadius_tm = triMesh.equivalentRadius();

//    float chromosomeVolume = regionAnalysisChromosomes.computeRegionFeature(REGION_FEATURE_VOLUME,originalVoxelMatrix)[numChromosome];
//    float eqRadius = regionAnalysisChromosomes.computeRegionFeature(REGION_FEATURE_EQUIVALENT_RADIUS,originalVoxelMatrix)[numChromosome];
//    float chromosomeRelativeVolume = regionAnalysisChromosomes.computeRegionFeature(REGION_FEATURE_VOLUME,originalVoxelMatrix)[numChromosome] / regionAnalysisChromosomes.computeRegionFeature(REGION_FEATURE_VOLUME,originalVoxelMatrix).sum();
//    float flatness = regionAnalysisChromosomes.computeRegionFeature(REGION_FEATURE_FLATNESS,originalVoxelMatrix)[numChromosome];
//    float elongation = regionAnalysisChromosomes.computeRegionFeature(REGION_FEATURE_ELONGATION,originalVoxelMatrix)[numChromosome];
//    float sphericity_vm = regionAnalysisChromosomes.computeRegionFeature(REGION_FEATURE_SPHERICITY,originalVoxelMatrix)[numChromosome];
//    float sphericity_tm = ( 36 * M_PI * pow(abs(triMesh.volume()) , 2) ) / pow( abs(triMesh.area() ) , 3);
//    float surfaceArea_vm = regionAnalysisChromosomes.computeRegionFeature(REGION_FEATURE_SURFACE_AREA,originalVoxelMatrix)[numChromosome];
//    float surfaceArea_tm = abs(triMesh.area() ) ;
//    float chromosomesIntensity = regionAnalysisChromosomes.computeRegionFeature(REGION_FEATURE_INTENSITY,originalVoxelMatrix)[numChromosome];
//    float chromosomesIntegratedDensity = regionAnalysisChromosomes.computeRegionFeature(REGION_FEATURE_INTEGRATED_DENSITY,originalVoxelMatrix)[numChromosome];
//    float relativeChromosomesIntensity = ( regionAnalysisChromosomes.computeRegionFeature(REGION_FEATURE_INTENSITY,originalVoxelMatrix)[numChromosome] ) / regionAnalysisChromosomes.computeRegionFeature(REGION_FEATURE_INTENSITY,originalVoxelMatrix).sum();

//    chromosomesDataset.setValue ( "equivalentRadius_vm", numChromosome+totalNumChromosomes, eqRadius );
//    chromosomesDataset.setValue ( "equivalentRadius_tm", numChromosome+totalNumChromosomes, eqRadius_tm );
//    chromosomesDataset.setValue ( "chromosomeVolume_vm", numChromosome+totalNumChromosomes, chromosomeVolume );
//    chromosomesDataset.setValue ( "chromosomeVolume_tm", numChromosome+totalNumChromosomes, chromosomeVolume_tm );
//    chromosomesDataset.setValue ( "distanceToTheBorder", numChromosome+totalNumChromosomes, distanceToBorder );
//    chromosomesDataset.setValue ( "chromosomeRelativeVolume", numChromosome+totalNumChromosomes, chromosomeRelativeVolume );
//    chromosomesDataset.setValue ( "flatness", numChromosome+totalNumChromosomes, flatness );
//    chromosomesDataset.setValue ( "elongation", numChromosome+totalNumChromosomes, elongation );
//    chromosomesDataset.setValue ( "sphericity_vm", numChromosome+totalNumChromosomes, sphericity_vm );
//    chromosomesDataset.setValue ( "sphericity_tm", numChromosome+totalNumChromosomes, sphericity_tm );
//    chromosomesDataset.setValue ( "surfaceArea_vm", numChromosome+totalNumChromosomes, surfaceArea_vm );
//    chromosomesDataset.setValue ( "surfaceArea_tm", numChromosome+totalNumChromosomes, surfaceArea_tm );
//    chromosomesDataset.setValue ( "chromosomesIntensity", numChromosome+totalNumChromosomes, chromosomesIntensity );
//    chromosomesDataset.setValue ( "chromosomesIntegratedDensity", numChromosome+totalNumChromosomes, chromosomesIntegratedDensity );
//    //this last relativeCCsIntensity is the intensity of each cc divided by the total chromosome's intensity
//    chromosomesDataset.setValue ( "relativeChromosomesIntensity", numChromosome+totalNumChromosomes, relativeChromosomesIntensity );

//    individualChromosomesDataset.setValue ( "id", numChromosome, filename );
//    individualChromosomesDataset.setValue ( "idChromosome", numChromosome, numChromosome+1 );
//    individualChromosomesDataset.setValue ( "centroidCoordX", numChromosome, centroid[X] );
//    individualChromosomesDataset.setValue ( "centroidCoordY", numChromosome, centroid[Y] );
//    individualChromosomesDataset.setValue ( "centroidCoordZ", numChromosome, centroid[Z] );
//    individualChromosomesDataset.setValue ( "equivalentRadius_vm", numChromosome, eqRadius );
//    individualChromosomesDataset.setValue ( "equivalentRadius_tm", numChromosome, eqRadius_tm );
//    individualChromosomesDataset.setValue ( "distanceToTheBorder", numChromosome, distanceToBorder );
//    individualChromosomesDataset.setValue ( "chromosomeVolume_vm", numChromosome, chromosomeVolume );
//    individualChromosomesDataset.setValue ( "chromosomeVolume_tm", numChromosome, chromosomeVolume_tm );
//    individualChromosomesDataset.setValue ( "chromosomeRelativeVolume", numChromosome, chromosomeRelativeVolume );
//    individualChromosomesDataset.setValue ( "flatness", numChromosome, flatness );
//    individualChromosomesDataset.setValue ( "elongation", numChromosome, elongation );
//    individualChromosomesDataset.setValue ( "sphericity_vm", numChromosome, sphericity_vm );
//    individualChromosomesDataset.setValue ( "sphericity_tm", numChromosome, sphericity_tm );
//    individualChromosomesDataset.setValue ( "surfaceArea_vm", numChromosome, surfaceArea_vm );
//    individualChromosomesDataset.setValue ( "surfaceArea_tm", numChromosome, surfaceArea_tm );
//    individualChromosomesDataset.setValue ( "chromosomeIntensity", numChromosome, chromosomesIntensity );
//    individualChromosomesDataset.setValue ( "chromosomeIntegratedDensity", numChromosome, chromosomesIntegratedDensity );
//    individualChromosomesDataset.setValue ( "relativeChromosomesIntensity", numChromosome, relativeChromosomesIntensity );

//  }

//  totalNumChromosomes += regionAnalysisChromosomes.numRegions();

  //return chromosomesDataset;
}
