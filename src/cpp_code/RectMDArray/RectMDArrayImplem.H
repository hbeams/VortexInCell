//#include "CH_Timer.H"
#include <cstdlib>
#include <functional>
#include <algorithm>
#include <iostream>

//Default Constructor
template <class T, unsigned int C> 
RectMDArray<T,C>::RectMDArray():m_rawPtr(NULL){}

//RectMDArray::define()
template <class T, unsigned int C> 
void RectMDArray<T,C>::define(const DBox& a_box)
{
  if (m_rawPtr)
    {
      delete [] m_rawPtr;
    }
  m_rawPtr = new T[a_box.sizeOf()*C];
  m_box = a_box;
}

//Construct a RectMDArray from a DBox. Calls RectMDArray::define()
template <class T, unsigned int C> 
RectMDArray<T,C>::RectMDArray(const DBox& a_box):m_rawPtr(NULL) 
{
  define(a_box);
}

//Destructor
template <class T, unsigned int C> 
RectMDArray<T,C>::~RectMDArray()
{
  delete [] m_rawPtr;
}
template <class T, unsigned int C> bool RectMDArray<T,C>::defined() const
{
  return bool(m_rawPtr);
}

template <class T, unsigned int C> void RectMDArray<T,C>::setVal(const T& a_val)
{
  if (m_rawPtr)
  {
    for (int k = 0; k < dataSize();k++)
    {
      m_rawPtr[k] = a_val;
    }
  }
}
template <class T, unsigned int C> 
void RectMDArray<T,C>::copyTo(RectMDArray<T,C>& a_dest) const
{
  
  DBox bxInt = m_box & a_dest.m_box; // bxInt is the intersection of m_box with a_dest.m_box
  Point hi = bxInt.getHighCorner();
  Point lo = bxInt.getLowCorner();
  int pencilLength = hi[0] - lo[0];
  hi[0] = lo[0]; //flatten the box in the pencil direction
  DBox pencilBox(lo,hi);
  const T* srcVal;
  T* destVal;
  for (Point pt = lo; pencilBox.notDone(pt); pencilBox.increment(pt))
  {
    for (unsigned int comp = 0; comp < C; comp++)
    {
      srcVal = (*this)(pt,comp);
      destVal = a_dest(pt,comp);
      for (int ii = 0; ii <= pencilLength; ii++)
      {
        destVal[ii] = srcVal[ii];
      }
    }
  }
}    
template <class T, unsigned int C> T& RectMDArray<T,C>::operator[](const Point& a_iv)
{
  static_assert(C==1,"operator[] only defined for scalar RectMDArray");
  int k = m_box.getIndex(a_iv);
  assert(k < m_box.sizeOf());
  assert(k >= 0);
  return m_rawPtr[k];
}
template <class T, unsigned int C> const T& RectMDArray<T,C>::operator[](const Point& a_iv) const
{
  static_assert(C==1,"operator[] only defined for scalar RectMDArray");    
  int k = m_box.getIndex(a_iv);
  assert(k < m_box.sizeOf());
  assert(k >= 0);
  return m_rawPtr[k];
}
template <class T, unsigned int C> T& RectMDArray<T,C>::operator()(const Point& a_iv, unsigned int a_comp)
{
  int m=m_box.sizeOf();
  int k = m_box.getIndex(a_iv);
  assert(k < m);
  assert(k >= 0);
  return m_rawPtr[k+m*a_comp];
}
template <class T, unsigned int C> const T& RectMDArray<T,C>::operator()(const Point& a_iv, unsigned int a_comp) const
{
  int m=m_box.sizeOf();
  int k = m_box.getIndex(a_iv);
  assert(k < m);
  assert(k >= 0);
  return m_rawPtr[k+m*a_comp];
}



