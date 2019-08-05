#ifndef _DBOX_H_
#define _DBOX_H_
#include "Point.H"
#ifdef USE_CHOMBO
#include "Box.H"
#endif
#include <iostream>
extern bool reportMemory;
extern unsigned long long int memory;

/// An interval in DIM dimensional space.
/**
  A DBox is a region specified by two corners, highCorner and lowCorner. In 2 dimensions, a DBox is a rectangular region. In 3 dimensions, a DBox is a rectangular prism, etc. 
*/
class DBox
{
public:

  /// Default constructor.
  /**
    The DBox's upper and lower corner default to (-1,-1,...,-1) and (0,0,...,0) respectively. The default size is -1 (empty).
  */
  DBox();

  /// Constructor for nontrivial DBox using two Points. 
  DBox(const Point& a_lowCorner,
      const Point& a_highCorner);

  /// Copy constructor.
  DBox(const DBox& a_DBox); 

#ifdef USE_CHOMBO
  /// Construction of a DBox from a Chombo Box.
  DBox(const Box& a_box);

#endif

  ///  Computes the intersection of the DBox with a_rightDBox.  DBox AND operation. 
  /**
    Returns a new DBox which is the intersection of *this and *a_rightDBox. *this is unchanged
  */
  DBox operator&(const DBox& a_rightDBox) const;

  /// Computes intersection of DBox and a_rightDBox in place.
  /**
    Modifies *this to be the intersection of *this and a_rightDBox.
  */
  void operator&=(const DBox& a_rightDBox);

  /// Computes shifted DBox by a_offset in direction a_direction.
  /**
    Does not change *this. Returns the shifted box.
  */
  DBox shift(int a_direction,
            int a_offset) const;

  /// Computes a DBox shifted according to the vector a_pt
  /**
    Does not change *this. Returns the shifted box. Each point is shifted by adding a_pt. 
   */
  DBox shift(const Point& a_pt) const;

  /// Grow in all of the coordinate directions by a_numpoints.
  /**
    Does not change *this. Returns the grown box. Returned DBox satisfies: m_upperCorner_new = m_upperCorner_old + getOnes()*a_numpoints, m_lowerCorner_new = m_lowerCorner_old - getOnes()*a_numpoints.
   */
  DBox grow(int a_numpoints) const;

  /// Grow in each coordinate direction by an amount given by the component of a_pt.
  /**
    Does not change *this. Returns the grown box. Returned DBox satisfies: m_upperCorner_new = m_upperCorner_old + a_pt, m_lowerCorner_new = m_lowerCorner_old - a_pt.
   */
  DBox grow(const Point& a_pt) const;

  /// coarsen in all of the coordinate directions by a_numpoints.
  /**
     does not change *this. returns the coarsened box
   */
  DBox coarsen(int a_numpoints) const;

  /// coarsen in each coordinate direction by an amount given by the component of a_pt.
  /**
     does not change *this. returns the coarsened box
   */
  DBox coarsen(const Point& a_pt) const;

  /// refine in all of the coordinate directions by a_numpoints.
  /**
     does not change *this. returns the refined box
   */
  DBox refine(int a_numpoints) const;

  /// refine in each coordinate direction by an amount given by the component of a_pt.
  /**
     does not change *this. returns the refined box
   */
  DBox refine(const Point& a_pt) const;

  /// refine in each coordinate direction by an amount given by the component of a_pt.
  /**
     does not change *this. returns the refined box
   */
  DBox refineCC(const Point& a_pt) const;

  // refine in all of the coordinate directions by a_numpoints.
  /**
     does not change *this. returns the refined box
   */
  DBox refineCC(int a_nref) const;

  /// Access a DBox's lowCorner.
  /**
    Called by a constant DBox. Returned object cannot be changed. 
  */
  inline const Point& getLowCorner() const {return m_lowCorner;};

  /// Access a DBox's highCorner.
  /**
    Called by a constant DBox. Returned object cannot be changed. 
  */
  inline const Point& getHighCorner() const {return m_highCorner;}; 

  /// Access a DBox's lowCorner.
  /**
    Called by a non-constant DBox. Returned object can be changed.
  */
  inline const Point& getLowCorner()  {return m_lowCorner;};

  /// Access a DBox's highCorner.
  /**
    Called by a non-constant DBox. Returned object can be changed.
  */
  inline const Point& getHighCorner() {return m_highCorner;}; 

  //REDUNDANT (and not called anywhere in the code)
  // /// Computes a linear index corresponding to the input array. \todo Needs implementation.
  // int linearIndex(const Point& a_pt) const;

  /// Computes the size, or "volume", of the DBox.
  const unsigned int& sizeOf() const{return m_size;};

  /// Defines equality between DBoxes.
  /**
    Two DBoxes are considered equal if they have identical (==) lowCorner and highCorner
  */
  bool operator==(const DBox& a_rhsDBox) const
  {return ((m_lowCorner == a_rhsDBox.m_lowCorner) 
           && (m_highCorner == a_rhsDBox.m_highCorner));}

  /// Checks if the DBox contains the Point a_pt.
  bool contains(const Point& a_pt) const;

  /// Checks if the DBox contains another DBox a_rhs.
  bool contains(const DBox& a_rhs) const{
    return this->contains(a_rhs.m_lowCorner) && this->contains(a_rhs.m_highCorner);
  }

  // REDUNDANT
  //***************************************************************************************************************
  // /// Checks if the DBox contains another DBox which is augmented by a_extent.
  // /**
  //   This function will test if the DBox with lowCorner = a_rhs.m_lowerCorner + a_extent.m_lowerCorner and highCorner = // a_rhs.m_highCorner + a_extent.m_highCorner resides within *this. \todo Determine the purpose of this function
  // */
  // bool contains(const DBox& a_rhs, const DBox& a_extent)
  // {
  //   return this->contains(a_rhs.m_lowCorner +a_extent.m_lowCorner) && 
  //     this->contains(a_rhs.m_highCorner+a_extent.m_highCorner);
  // }
  //***************************************************************************************************************

  // REDUNDANT
  //***************************************************************************************************************
  // /// Checks if the DBox contains the Point a_pt.
  // bool hasPoint(const Point& a_point) const;
  //***************************************************************************************************************

  /// Finds periodic image of input that is contained in the DBox.
  Point mod(const Point& a_pt) const;

  /// Get linear index of a Point in a DBox.
  /**
    Indices have values from 0 to m_size-1. 
  */
  inline unsigned int  getIndex(const Point& a_pt) const;

  /// Returns false if a_pt(k) > m_highCorner(k) for k in [0,Dim)
  inline bool notDone(const Point& a_pt) const;

  /// Returns the number of points that lie along a DBox's edge in direction a_dim.
  /**
    This is a measure of the "edge length" of the DBox. Note that the value returned is the number of points along the DBox's edge, NOT m_highCorner[a_dim] - m_lowCorner[a_dim]. 
  */
  unsigned int size(unsigned char a_dim) const
  {
    return m_highCorner[a_dim]-m_lowCorner[a_dim] + 1;
  }

  /// iteration through the points in a DBox. a_pt is incremented to the next point in the DBox.
  /**
    Currently, increment only works when DIM <= 3. When calling increment on m_highCorner, the resulting Point will no longer be inside *this (contains() will return FALSE). 
  */
  void increment(Point& a_pt) const;

  /// Get Point corresponding to a linear index in [0, ... sizeOf()-1] inside the DBox
  Point getPoint(unsigned int k) const;

  /// Prints the DBox to the command line. 
  void print() const; 

  /// Returns True if the DBox has non-positive "volume"
  bool isEmpty() const 
  {bool ret = (m_size < 1); return ret;}

private:  
  void recomputeSize(); ///< Used to reevaluate the size of the box when it is changed.
  Point m_lowCorner; ///< Point object containing the lower bounds of the DBox.
  Point m_highCorner;  ///< Point object containing the upper bounds of the DBox.
  unsigned int m_size; ///< "Volume" of the box.
};

/// Defines "<<" operator for DBoxes. Used for printing. 
ostream& operator<<(ostream& os, const DBox& a_box);

/////======================IMPLEMENTATION=============================/////

inline  bool DBox::notDone(const Point& a_pt) const
{
  bool retval = true;
  for(int idir = 0; idir < DIM; idir++)
    {
      retval = retval && (a_pt[idir] <= m_highCorner[idir]);
    }
  return retval;
}
inline unsigned int DBox::getIndex(const Point& a_pt) const
{
  unsigned int factor = 1;
  unsigned int linIndex = a_pt[0] - m_lowCorner[0];
  for (unsigned char i = 1;i < DIM;i++)
    { 
      factor = factor*(m_highCorner[i-1] - m_lowCorner[i-1]+1);
      linIndex = linIndex + (a_pt[i] - m_lowCorner[i])*factor;
    }
  return linIndex;
}
inline int mymod(int i,int j) 
{
  if (i < 0)
    {
      int k = (-1-i)%j;
      return j-1-k;
    }
  else
    {
      return i%j;
    }
}
void printPoint(const Point& a_pt);

#ifdef USE_CHOMBO
//BOX CONVERSION METHODS
Box makeBox(const DBox& a_box);

DBox makeDBox(const Box& a_box);

#endif

#endif
