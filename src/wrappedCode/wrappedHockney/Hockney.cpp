#include <cmath>
#include <complex>
#include <vector>
#include <cstdio>
#include <iostream>
#include <stdlib.h>
#include <memory>
#include "PowerItoI.H"
#include "RectMDArray.H"
#include "DBox.H"
#include "FFT1DW.H"
#include "FFT1D.H"
#include "FFTMD.H"
#include "Hockney.H"
   
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
Hockney::Hockney(double a_h,int a_M)
{ 
    shared_ptr<FFT1DW> p_fftw1d = shared_ptr<FFT1DW>(new FFT1DW(a_M+1));
    shared_ptr<FFT1D> p_fft = dynamic_pointer_cast<FFT1D>(p_fftw1d);
    m_fftmd.define(p_fft);
    m_h = a_h;
    m_M = a_M;
    m_N = Power(2,a_M);
    int pointInitTupleHigh[2], pointInitTupleLow[2];
    pointInitTupleHigh[0]=m_N; 
    pointInitTupleHigh[1]=m_N; 
    Point high(pointInitTupleHigh);
    pointInitTupleLow[0]=0; 
    pointInitTupleLow[1]=0; 
    Point low(pointInitTupleLow);
    m_domainbox = DBox( low, high);
};
void Hockney::define(double a_h,int a_M)
{ 
    shared_ptr<FFT1DW> p_fftw1d = shared_ptr<FFT1DW>(new FFT1DW(a_M+1));
    shared_ptr<FFT1D> p_fft = dynamic_pointer_cast<FFT1D>(p_fftw1d);
    m_fftmd.define(p_fft);
    m_h = a_h;
    m_M = a_M;
    m_N = Power(2,a_M);
    int pointInitTupleHigh[2], pointInitTupleLow[2];
    pointInitTupleHigh[0]=m_N; 
    pointInitTupleHigh[1]=m_N; 
    Point high(pointInitTupleHigh);
    pointInitTupleLow[0]=0; 
    pointInitTupleLow[1]=0; 
    Point low(pointInitTupleLow);
    m_domainbox = DBox( low, high);
};
void Hockney::getKernel(RectMDArray<complex<double> >& a_kerArray,double &a_h)
{
    DBox domain=m_domainbox;
    Point low = domain.getLowCorner();
    Point high = domain.getHighCorner();
    RectMDArray<double> kerArray(domain);
    complex<double> one(0.,0.);
    a_kerArray.setVal(one);

    for (Point pt = low;domain.notDone(pt);domain.increment(pt))
    {
        double rsq = (pow(pt[0]*m_h,2)
                +  pow(pt[1]*m_h,2));
        if (rsq == 0) 
        {
            a_kerArray[pt].real(0.0);
        }
        else
        {
            a_kerArray[pt].real(-log(sqrt(rsq))/(2*M_PI));
        }
    }

};
void Hockney::convolve(RectMDArray<double>& a_rhs)
{
    DBox rhsDomain = m_domainbox;
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
    getKernel(kernel,m_h);
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
}
void Hockney::convolveRMDA(RectMDArray<double>& a_rhs)
{
    streambuf *coutbuf = cout.rdbuf();
    ofstream outputStream;
    outputStream.open("rhsWrapped.txt", ofstream::app);
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
    ofstream outputStreamRHSDouble;
    outputStreamRHSDouble.open("RHSDoubleWrapped.txt", ofstream::app);
    streambuf *outputStreamBufRHSDouble = outputStreamRHSDouble.rdbuf();
    cout.rdbuf(outputStreamBufRHSDouble);
    for (Point pt = ddomain.getLowCorner();ddomain.notDone(pt);ddomain.increment(pt))
    {
        cout<< setprecision(16) <<pt[0]<<" "<<pt[1]<<" "<<real(rhsDouble[pt])<<" "<<imag(rhsDouble[pt])<<endl;
    }
    cout.rdbuf(coutbuf);
    RectMDArray<complex<double> > kernel(ddomain);
    RectMDArray<double > realOut(ddomain);
    getKernel(kernel,m_h);
    m_fftmd.forwardCCcen(rhsDouble);
    m_fftmd.forwardCCcen(kernel);
    ofstream outputStreamRHSDoubleFFT;
    outputStreamRHSDoubleFFT.open("RHSDoubleFFTWrapped.txt", ofstream::app);
    streambuf *outputStreamBufRHSDoubleFFT = outputStreamRHSDoubleFFT.rdbuf();
    cout.rdbuf(outputStreamBufRHSDoubleFFT);
    for (Point pt = ddomain.getLowCorner();ddomain.notDone(pt);ddomain.increment(pt))
    {
        cout<< setprecision(16) <<pt[0]<<" "<<pt[1]<<" "<<real(rhsDouble[pt])<<" "<<imag(rhsDouble[pt])<<endl;
    }
    cout.rdbuf(coutbuf);
    for (Point pt = ddomain.getLowCorner();ddomain.notDone(pt);ddomain.increment(pt))
    {
        rhsDouble[pt] *= kernel[pt];
    }
    ofstream outputStreamKernel;
    outputStreamKernel.open("KernelWrapped.txt", ofstream::app);
    streambuf *outputStreamBufKernel = outputStreamKernel.rdbuf();
    cout.rdbuf(outputStreamBufKernel);
    for (Point pt = ddomain.getLowCorner();ddomain.notDone(pt);ddomain.increment(pt))
    {
        cout<< setprecision(16) <<pt[0]<<" "<<pt[1]<<" "<<real(kernel[pt])<<" "<<imag(kernel[pt])<<endl;
    }
    cout.rdbuf(coutbuf);
    m_fftmd.inverseCCcen(rhsDouble);
    a_rhs.setVal(0.);
    DBox bx(rhsDomain.getLowCorner(),rhsDomain.getHighCorner() - getOnes());
    for (Point pt = bx.getLowCorner();bx.notDone(pt);bx.increment(pt))
    {
        a_rhs[pt] = real(rhsDouble[pt])*scale;
    }
    ofstream outputStreamSoln;
    outputStreamSoln.open("solnWrapped.txt", ofstream::app);
    streambuf *outputStreamBufSoln = outputStreamSoln.rdbuf();
    cout.rdbuf(outputStreamBufSoln);
//    cout<<"Soln M = "<< m_M<<endl;
    bxOut = a_rhs.getDBox();
    for (Point ptX=bxOut.getLowCorner(); bxOut.notDone(ptX); bxOut.increment(ptX))
    {
        cout<< setprecision(16) <<m_h*ptX[0]<<" "<<m_h*ptX[1]<<" "<<a_rhs[ptX]<<endl;
    }
    cout.rdbuf(coutbuf);
}
//This definitely needs a lot of testing
void Hockney::convertToRMDA(numpyArray<double> array)
{
    Ndarray<double, DIM> a(array);
    RectMDArray<double> writingToo(m_domainbox);   
    writingToo.setVal(0.0);
    for (Point p=m_domainbox.getLowCorner(); m_domainbox.notDone(p); m_domainbox.increment(p))
    {
        writingToo[p] = a[p[0]][p[1]];
        p.print2();
        cout<<writingToo[p]<<endl;
    }
}
void Hockney::convertToRMDA()
{
    cout<<"hello, world. This is Hockney!"<<endl;
}
DBox Hockney::getDBox()
{
    return m_domainbox;
}


