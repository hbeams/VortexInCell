# location of the Python header files
PYTHON_VERSION = 2.7
PYTHON_INCLUDE = /usr/include/python$(PYTHON_VERSION)
# location of the Boost Python include files and library
BOOST_INC = /usr/include
BOOST_LIB = /usr/local/lib
# compile mesh classes
TARGET = hello
$(TARGET).so: $(TARGET).o 
	clang++ -shared -Wl,-export_dynamic $(TARGET).o -L$(BOOST_LIB) -lboost_python27 -L/usr/lib/python$(PYTHON_VERSION)/config -lpython$(PYTHON_VERSION) -o $(TARGET).so
$(TARGET).o: $(TARGET).cpp
	clang++ -I$(PYTHON_INCLUDE) -I$(BOOST_INC) -fPIC -c $(TARGET).cpp

clean:
	rm -r *.o 
