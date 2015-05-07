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
    src/trimeshspatialmodel.cpp \
    src/spatialmodel.cpp \
    src/spatialmodelevaluator.cpp \
    src/spatialdescriptorborder.cpp \
    src/looptest.cpp \
    src/testrealdata.cpp \
    src/analysingchromocenters.cpp \
    src/analysingnuclei.cpp \
    src/fixingdata.cpp \
    src/isolatingnuclei.cpp \
    src/maximalrepulsion.cpp \
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
    src/maxrepulsionwithdistances.cpp \
    src/findnucleicascade.cpp \
    src/analysingnucleoli.cpp \
    src/spatialanalysisnucleoli.cpp \
    src/spatialdescriptorcentroid.cpp \
    src/getobjectsinterdistances.cpp

HEADERS += \
    src/trimeshspatialmodel.h \
    src/spatialdescriptorborder.h \
    src/maximalrepulsion.h \
    src/maxrepulsionwithdistances.h \
    src/spatialdescriptorcentroid.h

