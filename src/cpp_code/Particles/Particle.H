#ifndef _PARTICLE_H_
#define _PARTICLE_H_
#include <cstdio>
#include <cmath>
#include <cassert>
#include <array>
using namespace std;
/// Class representing the increment of a single point.
/**
ParticleShift is a vector of these, one for each particle.
*/
class DX
{
public:
  /// Default contructor sets values to zero.
  DX()
  {
    for(unsigned int i=0; i<DIM; i++) 
      {
	m_x[i]=0;
      }
  }
  /// Single data member is a DIM-tuple. 
  array<double, DIM> m_x;
  /// Increment is used to implement ParticleShift::increment. 
  inline void increment(double a_scale, const DX& a_rhs)
  {
    for(unsigned int i=0; i<DIM; i++) 
      {
	m_x[i]+=a_scale*a_rhs.m_x[i];
      }
  }
  /// Used to implement ParticleShift::scale.
  /**
     Does not have a data member corresponding to strength since the vorticity doesn't change.
   */ 
  inline void operator*=(double a_scale)
  {
    for(unsigned int i=0; i<DIM; i++) 
      {
	m_x[i]*=a_scale;
      }
  }
};
/// Represents a single particle in a ParticleSet.
class Particle
{
public:
  /// m_x stores the spatial location of the particle, and strength the vorticity.
  array<double, DIM> m_x;
  double strength;
  /// used to implement ParticleSet::increment.
  inline void increment(const DX& a_shift)
  {
    for(unsigned int i=0; i<DIM; i++) 
      {
	m_x[i]+=a_shift.m_x[i];
      }
  }
};

#endif
