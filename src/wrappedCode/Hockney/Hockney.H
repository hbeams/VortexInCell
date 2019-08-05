#ifndef _HOCKNEY_
#define _HOCKNEY_
#include "RectMDArray.H"
#include "FFTMD.H"
#include <array>
#include <memory>
#include <boost/python.hpp>

using namespace std;
class Hockney
{
public:
  Hockney();
  Hockney(double a_h,int a_M);
  void define(double a_h,int a_M);
  void convolve(double* a_rhs_ndArray);
  void convertToRMDA( double* a_rhs_ndArray);
  ~Hockney(){};
protected:
  void convolveRMDA(RectMDArray<double>& a_rhs);
  double m_h;
  int m_M,m_N;
  FFTMD m_fftmd;
  bool m_isDefined;
  DBox m_domainbox;
  void getKernel( RectMDArray<complex<double> >& a_kerArray, double& a_h);
  //void convertToRMDA( RectMDArray<double>& a_rhs, double* a_rhs_ndArray);
};


BOOST_PYTHON_MODULE(Hockney)
{
    using namespace boost::python;
    class_<Hockney>("Hockney", init<double, int>())
    .def("define", &Hockney::define); 
   // .def("convertToRMDA", &Hockney::convertToRMDA);
}
#endif
  
