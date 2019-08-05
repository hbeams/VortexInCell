#our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.
ANIMATION = 1
DIM=2

HOME=./../../../../..
RMDA=$(HOME)/src/wrappedCode/RectMDArray
FFT=$(HOME)/src/wrappedCode/fftTools
CONV=$(HOME)/src/wrappedCode/wrappedHockney
BOOST_INC=/usr/include
BOOST_NUMPY_INC=/usr/local/Cellar/boost/1.70.0/include
BOOST_LIB=/usr/local/lib
PYTHON_VERSION=3.7
PYTHON_INCLUDE=/anaconda3/include/python$(PYTHON_VERSION)m
FFTWDIR =/usr/local
VPATH= . $(RMDA) $(FFT) $(CONV) 

CXX=clang++
CFLAGS += -std=c++14 -I$(PYTHON_INCLUDE) -I$(BOOST_NUMPY_INC) -I$(BOOST_INC) -I$(RMDA) -I$(FFT) -I$(CONV) -I$(FFTWDIR)/include
CFLAGS += -D DIM=$(DIM) -D ANIMATION=$(ANIMATION)
CFLAGS += -fPIC 
LIBS:=  -L$(FFTWDIR)/lib -lfftw3 

SRCFILES:= $(notdir $(wildcard $(RMDA)/*.cpp $(FFT)/*.cpp $(CONV)/*.cpp ))

TARGET=Hockney

OBJS:=$(patsubst %.cpp,%.o,$(SRCFILES))

%.o:%.cpp ; $(CXX) -c $(CFLAGS) $< -o $@

objects: GNUmakefile $(OBJS) 

sharedObj:GNUmakefile ; $(CXX) -v -shared -Wl,-export_dynamic *.o $(LIBS) -L$(BOOST_LIB) -lboost_python37 -lboost_numpy37 -L/anaconda3/lib/python$(PYTHON_VERSION)/config-3.7m-darwin -lpython$(PYTHON_VERSION)m -o $(TARGET).so

clean:
	rm -r *.exe $(odir) $(ddir) *.o *.d *.so
