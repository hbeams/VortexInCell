#include "DBox.H"
#include <cassert>
///

bool reportMemory = false;
unsigned long long int memory=0;

DBox::DBox()
{
  m_lowCorner = getZeros();
  m_highCorner= getOnes();
  m_highCorner *= -1;
  m_size = 0;
}
///
DBox::DBox(const DBox& a_box)
{
  m_lowCorner = a_box.m_lowCorner;
  m_highCorner = a_box.m_highCorner;
  m_size = a_box.m_size;
}
#ifdef USE_CHOMBO
DBox::DBox(const Box& a_box)
{
  for (int dir = 0; dir < DIM; dir++)
    {
      m_lowCorner[dir] = a_box.smallEnd(dir);
      m_highCorner[dir] = a_box.bigEnd(dir);
      m_size = a_box.numPts();
    }
}
#endif
///
void
DBox::recomputeSize()
{
  m_size = 1;
  for(int idir = 0; idir < DIM; idir++)
    {
      m_size *= size(idir);
    }
}
DBox::DBox(const Point& a_lowCorner,const Point& a_highCorner)
{
  m_lowCorner = a_lowCorner;
  m_highCorner = a_highCorner;
  recomputeSize();
}
//DBox returns a new box, should have put const to make clear.
DBox DBox::operator&(const DBox& a_rightDBox) const
{
  int newLo[DIM];
  int newHi[DIM];
  for (int i = 0 ; i < DIM; i++)
    {
      newLo[i] = m_lowCorner[i];
      if (m_lowCorner[i] < a_rightDBox.m_lowCorner[i])
        {
          newLo[i] = a_rightDBox.m_lowCorner[i];
        }
      newHi[i] = m_highCorner[i];
      if (m_highCorner[i] > a_rightDBox.m_highCorner[i])
        {
          newHi[i] = a_rightDBox.m_highCorner[i];
        }
    }
  Point newLowCorner(newLo);
  Point newHighCorner(newHi); 
  for (int dir = 0; dir < DIM; dir++)
    {
      if (newHi[dir] < newLo[dir])
        {
          DBox ret0;
          return ret0;
        }
    }
  DBox ret(newLowCorner,newHighCorner);
  return ret;
}

void DBox::operator&=(const DBox& a_rightDBox) 
{
  DBox retval = *this & a_rightDBox;
  *this = retval;
}
DBox DBox::shift(int a_direction, int a_offset) const
{
  DBox returnDBox = DBox(*this);
  returnDBox.m_lowCorner += getUnitv(a_direction)*a_offset;
  returnDBox.m_highCorner += getUnitv(a_direction)*a_offset;
  return returnDBox;
}
DBox DBox::shift(const Point& a_pt) const
{
  DBox returnDBox = DBox(*this);
  returnDBox.m_lowCorner += a_pt;
  returnDBox.m_highCorner += a_pt;
  return returnDBox;
}
DBox DBox::grow(int a_offset) const
{
  Point lo = m_lowCorner;
  Point hi = m_highCorner;
  lo -= getOnes()*a_offset;
  hi += getOnes()*a_offset;
  DBox returnDBox(lo, hi);

  return returnDBox;
}
DBox DBox::grow(const Point& a_offset) const
{
  Point lo = m_lowCorner;
  Point hi = m_highCorner;
  lo -= a_offset;
  hi += a_offset;
  DBox returnDBox(lo, hi);

  return returnDBox;
}
DBox DBox::coarsen(int a_nref) const
{
  Point lo = m_lowCorner;
  Point hi = m_highCorner;
  lo /= a_nref;
  hi /= a_nref;
  DBox returnDBox(lo, hi);

  return returnDBox;
}
DBox DBox::coarsen(const Point& a_pt) const
{
  Point lo = m_lowCorner;
  Point hi = m_highCorner;
  lo /= a_pt;
  hi /= a_pt;
  DBox returnDBox(lo, hi);

  return returnDBox;
}
DBox DBox::refine(int a_nref) const
{

  Point lo = m_lowCorner;
  Point hi = m_highCorner;
  lo *= a_nref;
  hi += getOnes();
  hi *= a_nref;
  hi -= getOnes();
  DBox returnDBox(lo, hi);

  return returnDBox;
}
DBox DBox::refineCC(const Point& a_pt) const
{
  return refine(a_pt);
}
DBox DBox::refineCC(int a_nref) const
{
  return refine(a_nref);
}
DBox DBox::refine(const Point& a_pt) const
{
  Point lo = m_lowCorner;
  Point hi = m_highCorner;
  lo *= a_pt;
  hi += getOnes();
  hi *= a_pt;
  hi -= getOnes();
  DBox returnDBox(lo, hi);

  return returnDBox;
}

void DBox::increment(Point& a_pt) const
{
  Point current = a_pt;
  assert(DIM <= 4);

  current[0]++;
#if DIM >= 2
  if (current[0] > m_highCorner[0])
  {
    current[0] = m_lowCorner[0];
    current[1]++;
#if DIM >= 3
    if (current[1] > m_highCorner[1])
    {
      current[1] = m_lowCorner[1];
      current[2]++;
    }
#endif
#if DIM >= 4
    if (current[2] > m_highCorner[2])
    {
      current[2] = m_lowCorner[2];
      current[3]++;
    }
#endif
  }
#endif
  a_pt = current;

}


Point DBox::getPoint(unsigned int k) const
{
  int tuple[DIM];
  for (unsigned char i = 0;i < DIM; i++)
    {      
      int factor = (m_highCorner[i] - m_lowCorner[i] + 1);
      int kred = k%factor;
      tuple[i] = kred + m_lowCorner[i];
      k = (k - kred)/factor;
    }
  Point pt(tuple);
  return pt;
}

bool DBox::contains(const Point& a_pt) const
{
  bool retval = true;
  for(int idir = 0; idir < DIM; idir++)
    {
      if(a_pt[idir] < m_lowCorner[idir])
        {
          retval = false;
        }
      if(a_pt[idir] > m_highCorner[idir])
        {
          retval = false;
        }
    }
  return retval;
}

//Redundant
//***************************************************************************************************************
// bool DBox::hasPoint(const Point& a_point) const
// {
//   //check that point is in between lowCorner and highCorner
//   for (int i = 0; i < DIM; i++)
//   {
//     if (a_point[i] > m_highCorner[i])
//     {
//      return false;
//    }
//    else if (a_point[i] < m_lowCorner[i])
//    {
//      return false;
//    }
//  }
//  return true;
// }
//***************************************************************************************************************

void DBox::print() const 
{
  std::cout<<*this<<std::endl;
}

ostream& operator<<(ostream& os, const DBox& a_box)
{
  os << "[low Corner = " ;
  for (int k = 0;k < DIM; k++)
    {
      os << a_box.getLowCorner()[k]<< " ";
    }
  os << " high Corner = " ;
  for (int k = 0;k < DIM; k++)
    {
      os << a_box.getHighCorner()[k] << " ";
    }
  os<<" size="<<a_box.sizeOf();
  os<<"]";
  return os;
}
void printPoint(const Point& a_pt)
  {
    cout << "(";
    for (int dir = 0; dir < DIM ; dir++)
      {
        cout << a_pt[dir];
        if(dir < DIM-1) cout << ",";
      }
    cout <<   ")"  << endl;
  }
Point DBox::mod(const Point& a_pt) const
{
  int tuple[DIM];
  for (int i = 0;i< DIM; i++)
    {
      int dl = m_highCorner[i] - m_lowCorner[i] + 1;
      tuple[i] = mymod(a_pt[i],dl) + m_lowCorner[i];
    }
  return Point(tuple);
}
  
#ifdef USE_CHOMBO
//BOX CONVERSION METHODS
Box makeBox(const DBox& a_box)
{
    IntVect lo = makeIntVect(a_box.getLowCorner());
    IntVect hi = makeIntVect(a_box.getHighCorner());
    return Box(lo,hi);
}

DBox makeDBox(const Box& a_box)
{
    Point lo = makePoint(a_box.smallEnd());
    Point hi = makePoint(a_box.bigEnd());
    return DBox(lo,hi);
}
#endif
