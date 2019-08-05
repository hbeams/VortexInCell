# location of the Python header files
PYTHON_VERSION = 3.7
PYTHON_INCLUDE = /anaconda3/include/python$(PYTHON_VERSION)m
#PYTHON_INCLUDE = /usr/include/python$(PYTHON_VERSION)
# location of the Boost Python include files and library
BOOST_INC = /usr/include
BOOST_LIB = /usr/local/lib
# compile mesh classes
TARGET = hello_ext
$(TARGET).so: $(TARGET).o 
	clang++ -shared -Wl,-export_dynamic $(TARGET).o -L$(BOOST_LIB) -lboost_python37 -L/anaconda3/lib/python$(PYTHON_VERSION)/config-3.7m-darwin -lpython$(PYTHON_VERSION)m -o $(TARGET).so
$(TARGET).o: $(TARGET).C
	clang++ -I$(PYTHON_INCLUDE) -I$(BOOST_INC) -fPIC -c $(TARGET).C

clean:
	rm -r *.o 
