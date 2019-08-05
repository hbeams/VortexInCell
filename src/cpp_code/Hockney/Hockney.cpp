#include <cmath>
#include <complex>
#include <vector>
#include <cstdio>
#include <iostream>
#include <stdlib.h>
#include <memory>
#include "PowerItoI.H"
#include "RectMDArray.H"
#include "WriteRectMDArray.H"
#include "DBox.H"
#include "FFT1DW.H"
#include "FFT1D.H"
#include "FFTMD.H"
#include "Hockney.H"
#include "ConvKernel.H"

using namespace std;
void copyReal(RectMDArray<complex<double> >& a_cxarray,RectMDArray<double >& a_real)
{
  DBox d = a_real.getDBox()&a_cxarray.getDBox();
  for (Point pt=d.getLowCorner();d.notDone(pt);d.increment(pt))
    {
      a_real[pt] = real(a_cxarray[pt]);
    }
};
void copyImag(RectMDArray<complex<double> >& a_cxarray,RectMDArray<double >& a_imag)
{
  DBox d = a_imag.getDBox();
  for (Point pt=d.getLowCorner();d.notDone(pt);d.increment(pt))
    {
      a_imag[pt] = imag(a_cxarray[pt]);
    }
};
Hockney::Hockney()
{
  m_isDefined = false;
};
Hockney::Hockney(shared_ptr<ConvKernel>& a_kerPtr,const double& a_h,int a_M)
{ 
  shared_ptr<FFT1DW> p_fftw1d = shared_ptr<FFT1DW>(new FFT1DW(a_M+1));
  shared_ptr<FFT1D> p_fft = dynamic_pointer_cast<FFT1D>(p_fftw1d);
  m_fftmd.define(p_fft);
  m_kerPtr = a_kerPtr;
  m_h = a_h;
  m_M = a_M;
  m_N = Power(2,a_M);
};
void Hockney::define(shared_ptr<ConvKernel>& a_kerPtr,const double& a_h,int a_M)
{ 
  shared_ptr<FFT1DW> p_fftw1d = shared_ptr<FFT1DW>(new FFT1DW(a_M+1));
  shared_ptr<FFT1D> p_fft = dynamic_pointer_cast<FFT1D>(p_fftw1d);
  m_fftmd.define(p_fft);
  m_kerPtr = a_kerPtr;
  m_h = a_h;
  m_M = a_M;
  m_N = Power(2,a_M);
};
void Hockney::convolve(RectMDArray<double>& a_rhs)
{
    streambuf *coutbuf = cout.rdbuf();
    ofstream outputStream;
    outputStream.open("rhs.txt", ofstream::app);
    streambuf *outputStreamBuf = outputStream.rdbuf();
    cout.rdbuf(outputStreamBuf);
    DBox bxOut = a_rhs.getDBox();
    for (Point ptX=bxOut.getLowCorner(); bxOut.notDone(ptX); bxOut.increment(ptX))
    {
        cout<< setprecision(16) <<m_h*ptX[0]<<" "<<m_h*ptX[1]<<" "<<a_rhs[ptX]<<endl;
    }
    cout.rdbuf(coutbuf);
    
    DBox rhsDomain = a_rhs.getDBox();
    Point low = rhsDomain.getLowCorner();
    Point high = rhsDomain.getHighCorner();

    assert(low == getZeros());
    assert(high == getOnes()*m_N);

    low = high*(-1);
    DBox ddomain(low,high);
    RectMDArray<complex<double> > rhsDouble(ddomain);
    complex<double> zero(0.,0.);
    rhsDouble.setVal(zero);
    double scale = 1./pow(m_N*1.,DIM*2)/4;
    for (Point pt = rhsDomain.getLowCorner();rhsDomain.notDone(pt);rhsDomain.increment(pt))
    {
        rhsDouble[pt].real(a_rhs[pt]);
    }
    RectMDArray<complex<double> > kernel(ddomain);
    RectMDArray<double > realOut(ddomain);
    m_kerPtr->getKernel(kernel,m_h);
    m_fftmd.forwardCCcen(rhsDouble);
    m_fftmd.forwardCCcen(kernel);
    for (Point pt = ddomain.getLowCorner();ddomain.notDone(pt);ddomain.increment(pt))
    {
        rhsDouble[pt] *= kernel[pt];
    }
    m_fftmd.inverseCCcen(rhsDouble);
    a_rhs.setVal(0.);
    DBox bx(rhsDomain.getLowCorner(),rhsDomain.getHighCorner() - getOnes());
    for (Point pt = bx.getLowCorner();bx.notDone(pt);bx.increment(pt))
    {
        a_rhs[pt] = real(rhsDouble[pt])*scale;
    }
    
    outputStream.open("soln.txt", ofstream::app);
    outputStreamBuf = outputStream.rdbuf();
    cout.rdbuf(outputStreamBuf);
    bxOut = a_rhs.getDBox();
    for (Point ptX=bxOut.getLowCorner(); bxOut.notDone(ptX); bxOut.increment(ptX))
    {
        cout<< setprecision(16) <<m_h*ptX[0]<<" "<<m_h*ptX[1]<<" "<<a_rhs[ptX]<<endl;
    }
    cout.rdbuf(coutbuf);
}


