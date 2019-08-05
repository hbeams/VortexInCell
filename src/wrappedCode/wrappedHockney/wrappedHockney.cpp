#include "Hockney.H"
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace p = boost::python;
namespace n = p::numpy;

void convolve( double a_h, int a_M, n::ndarray& input_array)
{
    Hockney hock(a_h, a_M);
    
    DBox bx = hock.getDBox(); 
    RectMDArray<double> writingToo(bx);   
    writingToo.setVal(0.0);
    int stride = Power(2, a_M);
    double* ptr = reinterpret_cast<double*>(input_array.get_data());
    int transposed_tuple[DIM];
    //could use strides and divide by length of double instead
    //see https://github.com/ndarray/Boost.NumPy/blob/master/libs/numpy/example/wrap.cpp
    for (Point p=bx.getLowCorner(); bx.notDone(p); bx.increment(p))
    {
        transposed_tuple[0] = p[1];
        transposed_tuple[1] = p[0];
        Point transposedPoint(transposed_tuple);
        writingToo[p] = ptr[bx.getIndex(transposedPoint)];
    }
    hock.convolveRMDA(writingToo);
    for (Point p=bx.getLowCorner(); bx.notDone(p); bx.increment(p))
    {
        transposed_tuple[0] = p[1];
        transposed_tuple[1] = p[0];
        Point transposedPoint(transposed_tuple);
        ptr[bx.getIndex(transposedPoint)] = writingToo[p];
    }
}
BOOST_PYTHON_MODULE(Hockney)
{
    cout<<"Importing Hockney cpp"<<endl;
//    Py_Initialize();
//    cout<<"Py_Initialize Passed"<<endl;
    n::initialize();
    p::def("convolve", convolve);
}
