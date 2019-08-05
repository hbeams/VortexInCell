#ifndef _CONVKERNEL_H_
#define  _CONVKERNEL_H_
#include "RectMDArray.H"
#include <assert.h>
#include <cmath>
#include <complex>
#include <vector>
#include <cstdio>
#include <iostream>
using namespace std;
#include "PowerItoI.H"
class ConvKernel
{
public:
  /// Default constructor.
  ConvKernel(){};
  /// Destructor.
  virtual ~ConvKernel(){};
  /// Defines a default-constructed RectMDArray.  
  virtual void getKernel(RectMDArray<complex<double> >& a_srcArray,double& a_h) = 0; 
private:
};  
#endif
