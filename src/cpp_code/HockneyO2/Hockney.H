#ifndef _HOCKNEY_
#define _HOCKNEY_
#include "RectMDArray.H"
#include "FFTMD.H"
#include "ConvKernel.H"
#include <array>
#include <memory>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
using namespace std;
class Hockney
{
public:
  Hockney();
  Hockney(shared_ptr<ConvKernel>& a_kerPtr,const double& a_h,int a_M);
  void define(shared_ptr<ConvKernel>& a_kerPtr,const double& a_h,int a_M);
  void convolve(RectMDArray<double>& a_rhs);
  ~Hockney(){};
protected:
  double m_h;
  int m_M,m_N;
  FFTMD m_fftmd;
  shared_ptr<ConvKernel> m_kerPtr;
  bool m_isDefined;
};
#endif
  
