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
  -lfreed -limage -lshape -lgeneral -laline -lexception -ltiff

SOURCES += \
    src/findnucleus.cpp \
    src/main.cpp \
    src/applylabelling.cpp \
    src/findchromocenters.cpp \
    src/nucleusanalaysis.cpp \
    src/chromocentersanalysis.cpp \
    src/trimeshspatialmodel.cpp \
    src/spatialmodel.cpp \
    src/spatialmodelevaluator.cpp \
    src/spatialdescriptorborder.cpp \
    src/looptest.cpp \
    src/testrealdata.cpp

HEADERS += \
    src/trimeshspatialmodel.h \
    src/spatialdescriptorborder.h

