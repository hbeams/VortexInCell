#ifndef _RECTMDARRAY_H_
#define _RECTMDARRAY_H_
#include "DBox.H"
#ifdef USE_CHOMBO
#include "BaseFab.H"
#endif
#include <memory>
#include <cassert>
#include <functional>
//#include "PowerItoI.H"


template <class T=double, unsigned int C=1>
  class RectMDArray {
  public:
  /// Default constructor.
    RectMDArray();

  /// Constructs a RectMDArray over the Dbox a_box. Data is USUALLY initialized as zero.
    /**
      When using this constructor, it is recommended that the user initialize the data manually:
      @code
      Dbox B = Dbox(getZeros(),getOnes()*7);  //Domain
      RectMDArray<double> R = RectMDArray<double>(B); //Construct array
      R.setVal(0);  //Initialize array values as 0. 
      @endcode
    */
    RectMDArray(const DBox& a_box);

  /// Defines a default-constructed RectMDArray.
    /**
      Called by Constructor RectMDArray(const Dbox& a_box). Note that this can be called by an object that is already been defined. If that is the case, you need to delete the memory from the previous definition and allocate a new chunk of memory. This also means that the 
    */
    void define(const DBox& a_box);

  

  /// Destructor.
    ~RectMDArray(); 

  /// Sets all values in a RectMDArray to a constant value.
  /**
    \param a_val A constant value. 
  */
    void setVal(const T& a_val);

  /// Gets box over which array is defined. 
  
    DBox getDBox() const {return m_box;};

  /// Copy on Intersection.
  /**
    Copy the part of *this's data which intersects with the domain of a_dest's domain. 
  */
    void copyTo(RectMDArray<T,C>& a_dest) const;
  

  /// Indexing operator. only works for scalar RectMDArray objects
  /**
    Returns the data stored in *this corresponding to Point a_iv. Assertion error if a_iv is not in the domain of *this. 
    \param a_iv A Point in the domain of *this
  */
    inline T& operator[](const Point& a_iv);

  /// Indexing operator for const RectMDArrays. only works for scalar RectMDArray objects.
  /**
    Returns the data stored in *this corresponding to Point a_iv. Assertion error if a_iv is not in the domain of *this. 
    \param a_iv A Point in the domain of *this
  */
    inline const T& operator[](const Point& a_iv) const;

  /// Indexing operator for vector-valued RectMDArray objects with a single index. Assertion failure if returned type is not scalar.
  /**
    \param a_iv A Point in the domain of *this
    \param a_comp Integer corresponding to the desired component of the data at a_iv.
  */
    inline T& operator()(const Point& a_iv, unsigned int a_comp);

  /// Indexing operator for constant vector-valued RectMDArray objects with a single index. Assertion failure if returned type is not scalar.
  /**
    \param a_iv A Point in the domain of *this
    \param a_comp Integer corresponding to the desired component of the data at a_iv.
  */
  inline const T& operator()(const Point& a_iv, unsigned int a_comp) const;
  inline unsigned int dataSize()
  {return m_box.sizeOf()*C;};
  bool defined() const;
private:
  DBox m_box;
  T* m_rawPtr;
  // making copy contructor and assignment operator private.
  // We'll change this later.
  RectMDArray(const RectMDArray<T,C>& a_srcArray);
  RectMDArray<T,C> operator=(const RectMDArray<T,C>& a_srcArray);
};
#include "RectMDArrayImplem.H"  
#endif
