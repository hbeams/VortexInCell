#include "ParticleSet.H"

ParticleSet::ParticleSet(
              shared_ptr<ConvKernel>& a_kerptr,
              DBox& a_box,
              double& a_dx, 
              array<double, DIM>& a_lowCorner,
              int a_M):
  m_particles{},
  m_dx{a_dx},
  m_box{a_box},
  m_lowCorner{a_lowCorner}
{
  m_hockney.define(a_kerptr, a_dx, a_M);
}
void ParticleSet::increment(const ParticleShift& a_shift)
{
  for (int j = 0; j<a_shift.m_particles.size(); j++)
    {
      m_particles[j].increment(a_shift.m_particles[j]);
    }
}
