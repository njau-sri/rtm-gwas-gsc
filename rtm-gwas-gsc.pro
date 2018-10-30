TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -llapack -lblas -lgfortran -lquadmath
QMAKE_CXXFLAGS += -std=c++11 -fopenmp
QMAKE_LFLAGS += -static -fopenmp

SOURCES += \
        main.cpp \
    cmdline.cpp \
    lapack.cpp \
    util.cpp \
    vcf.cpp \
    rtm_gwas_gsc.cpp \
    pheno.cpp

HEADERS += \
    cmdline.h \
    lapack.h \
    split.h \
    util.h \
    vcf.h \
    pheno.h
