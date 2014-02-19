#############################################################################
# Makefile for building: bin/epitraits
# Generated by qmake (3.0) (Qt 5.0.2) on: Wed Feb 19 15:43:55 2014
# Project:  epitraits.pro
# Template: app
# Command: /home/javier/Qt5.0.2/5.0.2/gcc_64/bin/qmake -spec linux-g++-64 CONFIG+=debug CONFIG+=declarative_debug CONFIG+=qml_debug -o Makefile epitraits.pro
#############################################################################

MAKEFILE      = Makefile

first: release
install: release-install
uninstall: release-uninstall
QMAKE         = /home/javier/Qt5.0.2/5.0.2/gcc_64/bin/qmake
DEL_FILE      = rm -f
CHK_DIR_EXISTS= test -d
MKDIR         = mkdir -p
COPY          = cp -f
COPY_FILE     = cp -f
COPY_DIR      = cp -f -R
INSTALL_FILE  = install -m 644 -p
INSTALL_PROGRAM = install -m 755 -p
INSTALL_DIR   = $(COPY_DIR)
DEL_FILE      = rm -f
SYMLINK       = ln -f -s
DEL_DIR       = rmdir
MOVE          = mv -f
SUBTARGETS    =  \
		release \
		debug


release: FORCE
	$(MAKE) -f $(MAKEFILE).Release
release-make_first: FORCE
	$(MAKE) -f $(MAKEFILE).Release 
release-all: FORCE
	$(MAKE) -f $(MAKEFILE).Release all
release-clean: FORCE
	$(MAKE) -f $(MAKEFILE).Release clean
release-distclean: FORCE
	$(MAKE) -f $(MAKEFILE).Release distclean
release-install: FORCE
	$(MAKE) -f $(MAKEFILE).Release install
release-uninstall: FORCE
	$(MAKE) -f $(MAKEFILE).Release uninstall
debug: FORCE
	$(MAKE) -f $(MAKEFILE).Debug
debug-make_first: FORCE
	$(MAKE) -f $(MAKEFILE).Debug 
debug-all: FORCE
	$(MAKE) -f $(MAKEFILE).Debug all
debug-clean: FORCE
	$(MAKE) -f $(MAKEFILE).Debug clean
debug-distclean: FORCE
	$(MAKE) -f $(MAKEFILE).Debug distclean
debug-install: FORCE
	$(MAKE) -f $(MAKEFILE).Debug install
debug-uninstall: FORCE
	$(MAKE) -f $(MAKEFILE).Debug uninstall

Makefile: epitraits.pro ../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/linux-g++-64/qmake.conf ../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/features/spec_pre.prf \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/common/shell-unix.conf \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/common/unix.conf \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/common/linux.conf \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/common/gcc-base.conf \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/common/gcc-base-unix.conf \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/common/g++-base.conf \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/common/g++-unix.conf \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/qconfig.pri \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_bootstrap.pri \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_clucene.pri \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_concurrent.pri \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_core.pri \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_dbus.pri \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_declarative.pri \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_designer.pri \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_designercomponents.pri \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_gui.pri \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_help.pri \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_multimedia.pri \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_multimediawidgets.pri \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_network.pri \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_opengl.pri \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_platformsupport.pri \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_printsupport.pri \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_qml.pri \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_qmldevtools.pri \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_qmltest.pri \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_qtmultimediaquicktools.pri \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_quick.pri \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_quickparticles.pri \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_script.pri \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_scripttools.pri \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_sql.pri \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_svg.pri \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_testlib.pri \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_uitools.pri \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_v8.pri \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_webkit.pri \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_webkitwidgets.pri \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_widgets.pri \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_xml.pri \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_xmlpatterns.pri \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/features/qt_functions.prf \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/features/qt_config.prf \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/linux-g++-64/qmake.conf \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/features/spec_post.prf \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/features/exclusive_builds.prf \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/features/default_pre.prf \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/features/unix/default_pre.prf \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/features/resolve_config.prf \
		../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/features/default_post.prf \
		epitraits.pro
	$(QMAKE) -spec linux-g++-64 CONFIG+=debug CONFIG+=declarative_debug CONFIG+=qml_debug -o Makefile epitraits.pro
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/features/spec_pre.prf:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/common/shell-unix.conf:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/common/unix.conf:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/common/linux.conf:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/common/gcc-base.conf:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/common/gcc-base-unix.conf:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/common/g++-base.conf:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/common/g++-unix.conf:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/qconfig.pri:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_bootstrap.pri:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_clucene.pri:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_concurrent.pri:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_core.pri:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_dbus.pri:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_declarative.pri:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_designer.pri:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_designercomponents.pri:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_gui.pri:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_help.pri:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_multimedia.pri:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_multimediawidgets.pri:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_network.pri:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_opengl.pri:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_platformsupport.pri:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_printsupport.pri:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_qml.pri:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_qmldevtools.pri:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_qmltest.pri:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_qtmultimediaquicktools.pri:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_quick.pri:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_quickparticles.pri:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_script.pri:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_scripttools.pri:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_sql.pri:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_svg.pri:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_testlib.pri:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_uitools.pri:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_v8.pri:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_webkit.pri:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_webkitwidgets.pri:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_widgets.pri:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_xml.pri:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/modules/qt_lib_xmlpatterns.pri:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/features/qt_functions.prf:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/features/qt_config.prf:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/linux-g++-64/qmake.conf:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/features/spec_post.prf:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/features/exclusive_builds.prf:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/features/default_pre.prf:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/features/unix/default_pre.prf:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/features/resolve_config.prf:
../../../Qt5.0.2/5.0.2/gcc_64/mkspecs/features/default_post.prf:
epitraits.pro:
qmake: FORCE
	@$(QMAKE) -spec linux-g++-64 CONFIG+=debug CONFIG+=declarative_debug CONFIG+=qml_debug -o Makefile epitraits.pro

qmake_all: FORCE

make_first: release-make_first debug-make_first FORCE
all: release-all debug-all FORCE
clean: release-clean debug-clean FORCE
distclean: release-distclean debug-distclean FORCE
	-$(DEL_FILE) Makefile
FORCE:

$(MAKEFILE).Release: Makefile
$(MAKEFILE).Debug: Makefile
