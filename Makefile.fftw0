## /***************************************************************************
##   **************************************************************************
##   
##   SOFT: SO(3) Fourier Transforms
##   Version 2.0
## 
##   Copyright (c) 2003, 2004, 2007 Peter Kostelec, Dan Rockmore
##   
##   This file is part of SOFT.
## 
##   SOFT is free software; you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation; either version 3 of the License, or
##   (at your option) any later version.
## 
##   SOFT is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
## 
##   You should have received a copy of the GNU General Public License
##   along with this program.  If not, see <http://www.gnu.org/licenses/>.
##   
##   See the accompanying LICENSE file for details.
##   
##   ************************************************************************
##   ************************************************************************/


## define WALLCLOCK if want walltime and not cputime
## CFLAGS = -Wall -g -DWALLCLOCK
## CFLAGS = -Wall -g
## CFLAGS = -O3 -pg
## CFLAGS = -O3
## CFLAGS = -Wall -g -I./include
CFLAGS = -Wall -O3 -I./include

LIBNAME = libsoft0.a
LDFLAGS = -L. -lsoft0 -lm

COMSRC = $(wildcard common/*.c)
COMOBJ = $(patsubst %.c,%.o,$(COMSRC))

LIBSRC = $(COMSRC) \
	$(wildcard lib0/*.c)
LIBOBJ = $(patsubst %.c,%.o,$(LIBSRC))

EXELIST = test_wigSpec test_genWig test_Wigner_Analysis \
	test_Wigner_Synthesis test_Wigner_Naive test_Wigner_angle
EXELOC = $(patsubst %,bin/%,$(EXELIST))

EXE0LIST = test_soft test_soft_for test_soft_inv \
	test_soft_sym test_soft_sym_for test_soft_sym_inv \
	test_s2_rotate test_s2_rotate_wrap test_s2_rotate_mem \
	test_soft_sym_correlate test_soft_sym_correlate2 \
	test_soft_sym_correlate2_wrap
EXE0LOC = $(patsubst %,bin/%,$(EXE0LIST))

###################################################################
##################################################################
######
######              things that can be made
######
##################################################################
##################################################################

# things that can be made

depend:
	makedepend -I. $(SRC)

clean:
	rm -f {common,lib0,examples,examples0}/*.o

vclean:
	rm -f libsoft0.a $(EXELOC) $(EXE0LOC) {common,lib0,examples,examples0}/*.o

libsoft0.a: $(LIBOBJ)
	ar vru $(LIBNAME) $(LIBOBJ)
	ranlib $(LIBNAME)

lib: libsoft0.a

#
# making the example Wigner routines
#

examples: $(EXELIST) 


#
# making the example SO(3) FT routines
#

examples0: $(EXE0LIST)

tests: examples examples0

all: lib examples examples0

########################################################
#
# making the executables
#
########################################################

#
# discrete Wigner transform routines
#

test_wigSpec: examples/test_wigSpec.o libsoft0.a
	$(CC) $(CFLAGS) examples/test_wigSpec.o \
	$(LDFLAGS) -o bin/test_wigSpec

test_genWig: examples/test_genWig.o libsoft0.a
	$(CC) examples/test_genWig.o \
	$(LDFLAGS) -o bin/test_genWig

test_Wigner_angle: examples/test_Wigner_angle.o libsoft0.a
	$(CC) $(CFLAGS) examples/test_Wigner_angle.o \
	$(LDFLAGS) -o bin/test_Wigner_angle

test_Wigner_Analysis: examples/test_Wigner_Analysis.o libsoft0.a
	$(CC) $(CFLAGS) examples/test_Wigner_Analysis.o \
	$(LDFLAGS) -o bin/test_Wigner_Analysis

test_Wigner_Synthesis: examples/test_Wigner_Synthesis.o libsoft0.a
	$(CC) $(CFLAGS) examples/test_Wigner_Synthesis.o \
	$(LDFLAGS) -o bin/test_Wigner_Synthesis

test_Wigner_Naive: examples/test_Wigner_Naive.o libsoft0.a
	$(CC) $(CFLAGS) examples/test_Wigner_Naive.o \
	$(LDFLAGS) -o bin/test_Wigner_Naive

#
# SO(3) routines that do not depend on FFTW
#

test_soft: examples0/test_soft.o libsoft0.a
	$(CC) examples0/test_soft.o \
	$(LDFLAGS) -o bin/test_soft

test_soft_sym: examples0/test_soft_sym.o libsoft0.a
	$(CC) $(CFLAGS) examples0/test_soft_sym.o \
	$(LDFLAGS) -o bin/test_soft_sym

test_soft_for: examples0/test_soft_for.o libsoft0.a
	$(CC) $(CFLAGS) examples0/test_soft_for.o \
	$(LDFLAGS) -o bin/test_soft_for

test_soft_sym_for: examples0/test_soft_sym_for.o libsoft0.a
	$(CC) $(CFLAGS) examples0/test_soft_sym_for.o \
	$(LDFLAGS) -o bin/test_soft_sym_for

test_soft_inv: examples0/test_soft_inv.o libsoft0.a
	$(CC) $(CFLAGS) examples0/test_soft_inv.o \
	$(LDFLAGS) -o bin/test_soft_inv

test_soft_sym_inv: examples0/test_soft_sym_inv.o libsoft0.a
	$(CC) $(CFLAGS) examples0/test_soft_sym_inv.o \
	$(LDFLAGS) -o bin/test_soft_sym_inv

test_soft_sym_correlate: examples0/test_soft_sym_correlate.o libsoft0.a
	$(CC) $(CFLAGS)  examples0/test_soft_sym_correlate.o \
	$(LDFLAGS) -o bin/test_soft_sym_correlate

test_soft_sym_correlate2: examples0/test_soft_sym_correlate2.o libsoft0.a
	$(CC) $(CFLAGS) examples0/test_soft_sym_correlate2.o \
	$(LDFLAGS) -o bin/test_soft_sym_correlate2

test_soft_sym_correlate2_wrap: examples0/test_soft_sym_correlate2_wrap.o libsoft0.a
	$(CC) $(CFLAGS) examples0/test_soft_sym_correlate2_wrap.o \
	$(LDFLAGS) -o bin/test_soft_sym_correlate2_wrap

test_s2_rotate: examples0/test_s2_rotate.o libsoft0.a
	$(CC) $(CFLAGS) examples0/test_s2_rotate.o \
	$(LDFLAGS) -o bin/test_s2_rotate

test_s2_rotate_wrap: examples0/test_s2_rotate_wrap.o libsoft0.a
	$(CC) $(CFLAGS) examples0/test_s2_rotate_wrap.o \
	$(LDFLAGS) -o bin/test_s2_rotate_wrap

test_s2_rotate_mem: examples0/test_s2_rotate_mem.o libsoft0.a
	$(CC) $(CFLAGS) examples0/test_s2_rotate_mem.o \
	$(LDFLAGS) -o bin/test_s2_rotate_mem




# and now for LOTS OF dependencies ...

# DO NOT DELETE THIS LINE -- make depend depends on it.

