#ifndef _CUTOFFKERNEL_H_
#define  _CUTOFFKERNEL_H_
#include <cmath>
#include <complex>
#include "ConvKernel.H"
#include "WriteRectMDArray.H"
class CutoffKernel:public ConvKernel
{
public:
  //Default constructor.
  CutoffKernel(){};
  virtual ~CutoffKernel(){};
  // Constructor.
  CutoffKernel(double& a_h,double a_delta)
  {
    m_h = a_h;
    m_delta = a_delta;
  }
  virtual void getKernel(RectMDArray<complex<double> >& a_kerArray,double &a_h)
  {
    DBox domain=a_kerArray.getDBox();
    Point low = domain.getLowCorner();
    Point high = domain.getHighCorner();
    RectMDArray<double> kerArray(domain);
    double deltasq = pow(m_delta,2);
    complex<double> one(0.,0.);
    a_kerArray.setVal(one);

    for (Point pt = low;domain.notDone(pt);domain.increment(pt))
        {
          double rsq = (pow(pt[0]*m_h,2)
                          +  pow(pt[1]*m_h,2))/deltasq;
          if (rsq < 1) 

            {
              double r = sqrt(rsq);
              a_kerArray[pt].real((((((2160*r - 9800)*r + 16464)*r 
                                - 11025)*r*r + 2940)*r*r - 739.)/420/(2*M_PI));
            }
          else
            {
              a_kerArray[pt].real(log(rsq)/(4*M_PI));
            }
        }
             
  };
private:
  double m_h, m_delta;
};
#endif
