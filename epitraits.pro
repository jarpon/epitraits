TEMPLATE = app
DESTDIR = bin/

CONFIG = debug_and_release release
QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS += -fopenmp
QT =

FREEDLIBS = $(HOME)/src/cpp/freedlibs/

INCLUDEPATH = \
  $$FREEDLIBS/libexception/src \
  $$FREEDLIBS/libaline/src \
  $$FREEDLIBS/libgeneral/src \
  $$FREEDLIBS/libshape/src \
  $$FREEDLIBS/libimage/src \
  $$FREEDLIBS/libfreed/src

LIBS = \
  -L$$FREEDLIBS/libexception/src \
  -L$$FREEDLIBS/libaline/src  \
  -L$$FREEDLIBS/libgeneral/src \
  -L$$FREEDLIBS/libshape/src  \
  -L$$FREEDLIBS/libimage/src  \
  -L$$FREEDLIBS/libfreed/src \
#  -lfreed -limage -lshape -lgeneral -laline -lexception -ltiff
  -limage -lshape -lgeneral -laline -lexception -ltiff

SOURCES += \
    src/main.cpp \
    src/applylabelling.cpp \
    src/spatialmodel.cpp \
    src/spatialdescriptorborder.cpp \
    src/looptest.cpp \
    src/testrealdata.cpp \
    src/analysingchromocenters.cpp \
    src/analysingnuclei.cpp \
    src/fixingdata.cpp \
    src/isolatingnuclei.cpp \
    src/findgenes.cpp \
    src/findchromocenters.cpp \
    src/findnuclei.cpp \
    src/findchromosomes.cpp \
    src/analysingchromosomes.cpp \
    src/findnucleialternative.cpp \
    src/findchromocentersmanually.cpp \
    src/findchromocenters16bits.cpp \
    src/findmorenuclei.cpp \
    src/findnucleoli.cpp \
    src/findnucleicascade.cpp \
    src/analysingnucleoli.cpp \
    src/spatialanalysisnucleoli.cpp \
    src/spatialdescriptorcentroid.cpp \
    src/getobjectsinterdistances.cpp \
    src/spatialmodelboundaryinteraction.cpp \
    src/thresholding.cpp \
    src/otsuthresholding.cpp \
    src/voxelmatrix.cpp \
    src/testsStatiscalTest.cpp \
    src/spatialdescriptorborder2D.cpp \
    src/spatialmodeleval.cpp \
    src/unifyingobjects.cpp \
    src/spatialdescriptorfunctionhh.cpp \
    src/spatialdescriptorfunctiongg.cpp \
    src/findnuclei_old.cpp \
    src/testrealdatamorecompartments.cpp \
    src/spatialdescriptorfunctionm.cpp \
    src/regionanalysis2d.cpp \
    src/regionanalysis3d.cpp \
    src/spatialmodelevaluator2.cpp \
    src/cdftools2.cpp \
    src/comparedistributions.cpp \
    src/trimeshspatialmodeldifferentcompartments.cpp \
    src/spatialmodelcompleterandomness3ddifferentcompartments.cpp \
    src/spatialmodelhardcoredistance3ddifferentcompartments.cpp \
    src/spatialmodelhardcoreborderdistance3ddifferentcompartments.cpp \
    src/spatialmodelborderdistance3ddifferentcompartments.cpp \
    src/spatialmodel2.cpp \
    src/getnucleusprojection.cpp \
    src/spatialdescriptorfunctionz.cpp \
    src/generatepatterns.cpp \
    src/testrealdatawithexistingpatterns.cpp \
    src/spatialdescriptorfunctionlrd.cpp \
    src/spatialdescriptorfunctionalrd.cpp \
    src/spatialdescriptorfunctionnn.cpp \
    src/spatialmodelmaximalrepulsion3d2.cpp \
    src/spatialdescriptorfunctionsrd.cpp \
    src/spatialdescriptorfunctionasrd.cpp \
    src/generatepatternslessobjects.cpp

HEADERS += \
    src/spatialdescriptorborder.h \
    src/spatialdescriptorcentroid.h \
    src/spatialmodelboundaryinteraction.h \
    src/thresholding.h \
    src/regionanalysis2.h \
    src/otsuthresholding.h \
    src/voxelmatrix.h \
    src/spatialdescriptorborder2D.h \
    src/spatialdescriptorfunctionhh.h \
    src/spatialdescriptorfunctiongg.h \
    src/spatialdescriptorfunctionm.h \
    src/regionanalysis2d.h \
    src/regionanalysis3d.h \
    src/spatialmodelevaluator2.h \
    src/cdftools2.h \
    src/trimeshspatialmodeldifferentcompartments.h \
    src/spatialmodelcompleterandomness3ddifferentcompartments.h \
    src/spatialmodelhardcoreborderdistance3ddifferentcompartments.h \
    src/spatialmodelborderdistance3ddifferentcompartments.h \
    src/spatialmodelhardcoredistance3ddifferentcompartments.h \
    src/spatialmodel2.h \
    src/spatialdescriptorfunctionz.h \
    src/spatialdescriptorfunctionlrd.h \
    src/spatialdescriptorfunctionalrd.h \
    src/spatialdescriptorfunctionnn.h \
    src/spatialmodelmaximalrepulsion3d2.h \
    src/spatialdescriptorfunctionsrd.h \
    src/spatialdescriptorfunctionasrd.h

