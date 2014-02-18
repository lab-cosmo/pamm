# Makefile for the HB Gaussian Mixture library and tests
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
MODULES=hbmixture.f90 xyz.f90  gaussmix.f90
OBJECTS=$(MODULES:.f90=.o)

all: minitest gmm

minitest: $(OBJECTS) minitest.o
	gfortran $(FFLAGS) -o minitest $(OBJECTS) minitest.o

gmm: $(OBJECTS) gmm.o
	gfortran $(FFLAGS) -o gmm $(OBJECTS) gmm.o

%.o: %.f90
	gfortran $(FFLAGS) -c $<

clean:
	rm *.o *.mod
