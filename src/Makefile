# Makefile for the driver tests
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
# 
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 
.PHONY: all

FFLAGS=-g -O3
SOURCES=$(shell ls *.f90)
MODULES=distance.f90 matrixinverse.f90 gaussian.f90 hbmixture.f90 xyz.f90  gaussmix.f90
#OBJECTS=$(SOURCES:.f90=.o)
OBJECTSA=distance.o matrixinverse.o gaussian.o hbmixture.o xyz.o minitest.o
OBJECTSB=distance.o matrixinverse.o gaussian.o gaussmix.o gmm.o
all: modules minitest gmm

modules: $(MODULES)
	gfortran $(FFLAGS) -c $(MODULES)

minitest: $(OBJECTSA)
	gfortran $(FFLAGS) -o minitest $(OBJECTSA)

gmm: $(OBJECTSB)
	gfortran $(FFLAGS) -o gmm $(OBJECTSB)

%.o: %.f90
	gfortran $(FFLAGS) -c $<

clean:
	rm *.o *.mod
