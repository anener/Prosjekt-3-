TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        gauss_laguerre.cpp \
        gauss_legendre.cpp \
        lib.cpp \
        main.cpp \
        montecarlo.cpp

HEADERS += \
    catch.hpp \
    gauss_laguerre.h \
    gauss_legendre.h \
    lib.h \
    montecarlo.h

QMAKE_CXXFLAGS += -fopenmp
LIBS += -fopenmp


