#ifndef _HOCKNEY_
#define _HOCKNEY_
#include "RectMDArray.H"
#include "FFTMD.H"
#include <array>
#include <memory>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include "ndarray.h"

using namespace std;
class Hockney
{
public:
  Hockney();
  Hockney(double a_h,int a_M);
  void define(double a_h,int a_M);
  void convolve(RectMDArray<double>& a_rhs);
  void convolveRMDA(RectMDArray<double>& a_rhs);
  void convertToRMDA(numpyArray<double> array);
  void convertToRMDA();
  DBox getDBox();
  ~Hockney(){};
protected:
  double m_h;
  int m_M,m_N;
  FFTMD m_fftmd;
  bool m_isDefined;
  DBox m_domainbox;
  void getKernel( RectMDArray<complex<double> >& a_kerArray, double& a_h);
  //void convertToRMDA( RectMDArray<double>& a_rhs, double* a_rhs_ndArray);
};
#endif
  
