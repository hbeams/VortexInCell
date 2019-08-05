#our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.
DIM=2
VPATH = ../timer ../RectMDArray ../Writers
system := $(shell uname)
FFTWDIR = /usr/local
FFTLIBLOC = ../lib
CXX=clang++
#CXX=g++
FFTINCFLAGS = -I$(FFTWDIR)/include
#CFLAGS = -g -Wall -DDIM=$(DIM)
CFLAGS = -O3 -DDIM=$(DIM)
CFLAGS += -std=c++11 -I. -I../Writers -I../timer -I../RectMDArray $(FFTINCFLAGS)
ifeq ($(CXX) , clang++)
CFLAGS += -stdlib=libc++
endif
LIBS =  -L$(FFTWDIR)/lib -lfftw3 
#CXX = /Users/colella/Desktop/gcc-4.7-bin/usr/local/bin/g++ 

LIBTOOL = ar cru 
ifeq ($(system),Darwin)
  LIBTOOL = libtool -static -s -o
endif

osuffix:=$(DIM)D.o
dsuffix:=$(DIM)D.d
%.o: %.cpp GNUmakefile
	$(CXX) -c $(CFLAGS) $< -o $@
	$(CXX) -MM $(CFLAGS) $< > $*.d

%$(osuffix): %.cpp GNUmakefile
	$(CXX) -c $(CFLAGS) $< -o $@
	$(CXX) -MM $(CFLAGS) $< > $*$(dsuffix)

1DFFTOBJS = FFT1DBRI.o FFTCTBRI.o FFTW1D.o FFT1DRecursive.o 
#add these files when they are implemented.

ifneq ($(DIM),1)
MDFFTOBJS = FFTMD$(osuffix)
libfft$(DIM)D.a: GNUmakefile $(MDFFTOBJS)
	$(LIBTOOL) libfft$(DIM)D.a $(MDFFTOBJS) 
	mkdir -p ../lib;mv libfft$(DIM)D.a ../lib
endif

libfft1D.a: GNUmakefile $(1DFFTOBJS)
	$(LIBTOOL) libfft1D.a $(1DFFTOBJS)
	mkdir -p ../lib;mv libfft1D.a ../lib

clean:
	rm *.o *.d *.a
-include $(MDFFTOBJS:.o=.d)
#mkdir -p $(FFTLIBLOC);mv libfft1D.a ../lib
#mkdir -p $(FFTLIBLOC);mv libfft$(DIM)D.a ../lib