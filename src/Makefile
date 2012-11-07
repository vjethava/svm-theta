# This Makefile is used under Linux

SVMDIR=../lib/libsvm-3.1
SVMMATDIR=$(SVMDIR)/matlab

# for Mac
# MATLABDIR ?= /opt/local/matlab

CXX ?= g++
# CC = g++-4.2
# CXX = g++-4.2
CFLAGS = -Wall -Wno-deprecated -Wconversion -O3 -fPIC -I${MATLAB_HOME}/extern/include -I$(SVMDIR) -I$(SVMMATDIR) 
MEX = ${MATLAB_HOME}/bin/mex
MEX_OPTION = CC\#$(CXX) CXX\#$(CXX) LD\#$(CXX) CFLAGS\#"$(CFLAGS)" CXXFLAGS\#"$(CFLAGS)"
# comment the following line if you use MATLAB on 32-bit computer
MEX_OPTION += -largeArrayDims
MEX_EXT = $(shell ${MATLAB_HOME}/bin/mexext)

OCTAVEDIR ?= /usr/include/octave
OCTAVE_MEX = env CC=$(CXX) mkoctfile
OCTAVE_MEX_OPTION = --mex
OCTAVE_MEX_EXT = mex
OCTAVE_CFLAGS = -Wall -O3 -fPIC -I$(OCTAVEDIR) -I..

all:	modelread kernelwrite

modelread: svmReadModel.c 
	$(MEX) $(MEX_OPTION)  svmReadModel.c $(SVMDIR)/svm.o $(SVMMATDIR)/svm_model_matlab.o

kernelwrite: svmWriteKernel.c
	$(MEX) $(MEX_OPTION)  svmWriteKernel.c
clean:
	rm -f *~ *.o *.mex* *.obj
