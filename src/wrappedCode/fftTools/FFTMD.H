#ifndef _FFTMD_H_
#define _FFTMD_H_
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <memory>
#include "DBox.H"
#include "RectMDArray.H"
#include "FFT1D.H"
using namespace std;
class FFTMD
{
public:
  FFTMD(){};
  FFTMD(shared_ptr<FFT1D > a_fft1dPtr);
  void define(shared_ptr<FFT1D > a_fft1dPtr);
  void forwardCC(RectMDArray<complex<double> > & a_f) const;
  void inverseCC(RectMDArray<complex<double> > & a_fHat) const;
  void forwardCCcen(RectMDArray<complex<double> > & a_f) const;
  void inverseCCcen(RectMDArray<complex<double> > & a_fHat) const;
  const int& getN() const;
  const int& getM() const;
private:
  int m_N;
  int m_M;
  shared_ptr<FFT1D> m_fft1dPtr;
};
#endif
