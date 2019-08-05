#ifndef _PARTICLESET_H_
#define _PARTICLESET_H_

#include <vector>
#include "Particle.H"
#include "RectMDArray.H"
#include "Hockney.H"

using namespace std;

class ParticleSet;
/// Represents an increment in a ParticleSet conforming to the dX class in RK4. 
class ParticleShift
{
public:
  /// vector of increments of particles.
  vector<DX> m_particles;
  ParticleShift(){};
  ~ParticleShift(){};
  /// Sets length of m_particles, initializes to zero via default constructor of DX.
  void init(const ParticleSet& a_particles);
  /// m_particles[k] += a_rhs.m_particles[k]*a_scale.
  void increment(
                 double a_scale, 
                 const ParticleShift& a_rhs);
  /// m_particles[k] *= a_scale
  void operator*=(double a_scale);
  /// reinitializes the values m_particles[k] to zero. Not used in RK4.
  void setToZero();
};
/// Represents a collection of particles, conforming to the X class in RK4. Also contains the persistent data required to compute the right-hand side in ParticlesVelocities.
class ParticleSet
{
public:
  /// constructor. Called by main.
  /**
     a_kerptr: required to construct Hockney.
     a_box: DBox over which convolution is performed.
     a_dx: mesh spacing of PIC grid.
     a_lowCorner: DIM-tuple giving the location in physical space of the lowCorner of the PIC grid.
     a_M: N = 2^a_M is the number of grid points in the PIC grid.
   */
  ParticleSet(
              shared_ptr<ConvKernel>& a_kerptr,
              DBox& a_box,
              double& a_dx, 
              array<double, DIM>& a_lowCorner,
              int a_M);

  ParticleSet(){};
  ~ParticleSet(){};
  /// Main container for particle data.
  vector<Particle> m_particles;
  /// Mesh spacing.
  double m_dx;
  /// PIC grid.
  DBox m_box;
  /// Location in physical space of m_box.m_lowCorner;
  array<double,DIM> m_lowCorner;
  /// Hockney class that performs discrete convolution.
  Hockney m_hockney;
  /// inplements increment function required for RK$.
  void increment(const ParticleShift& a_shift);

};
#endif
