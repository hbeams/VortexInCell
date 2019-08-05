#include "ParticleVelocities.H"

ParticleVelocities::ParticleVelocities()
{}
void ParticleVelocities::operator()(ParticleShift& a_k, 
                     const double& a_time, const double& dt, 
                     ParticleSet& a_state)
{
  //First need to compute a_state + a_k
  vector<Particle> t_particles = a_state.m_particles;
  for (int j = 0; j<a_k.m_particles.size(); j++)
    {
      t_particles[j].increment(a_k.m_particles[j]);
    }
  //Then use interpolation made up of the interpolating function given
  RectMDArray<double> omegaG(a_state.m_box);
  omegaG.setVal(0.0);
  DBox Ubox = a_state.m_box.grow(-1);
  RectMDArray<double, DIM> UG(Ubox);
  UG.setVal(0.0);
  double h = a_state.m_dx;
  Point e0 = getUnitv(0);
  Point e1 = getUnitv(1);
  for (int k =0; k<t_particles.size(); k++)
    {
      array<int,DIM> ipos;
      array<double, DIM> sk;
      for (int l = 0; l < DIM; l++)
        {
	  double pos = t_particles[k].m_x[l];
	  ipos[l] = floor(pos/h);
	  sk[l] = (pos - ipos[l]*h)/h;
	}
      Point pt(ipos);
      omegaG[pt] += t_particles[k].strength*(1 - sk[0])*(1 - sk[1]);  
      omegaG[pt + e0] += t_particles[k].strength*sk[0]*(1 - sk[1]);
      omegaG[pt + e1] += t_particles[k].strength*sk[1]*(1 - sk[0]);
      omegaG[pt + e0 + e1] += t_particles[k].strength*sk[0]*sk[1];
      }
  //Call Hockney
  a_state.m_hockney.convolve(omegaG);
  //Finite Difference
  for (Point p=Ubox.getLowCorner(); Ubox.notDone(p); Ubox.increment(p))
        {
	  UG( p, 0)=(omegaG[p+e1] - omegaG[p-e1])/(2*h);
	  UG( p, 1)= -1.0*(omegaG[p+e0] - omegaG[p-e0])/(2*h);
	}
  //Interpolate back and return particle fields in a_k
  a_k.setToZero();
  for (int k =0; k<t_particles.size(); k++)
    {
      array<int,DIM> ipos;
      array<double,DIM> sk;
      for (int l = 0; l < DIM; l++)
        {
	  double pos = t_particles[k].m_x[l];
	  ipos[l] = floor(pos/h);
	  sk[l] = (pos - ipos[l]*h)/h;
	}
      Point pt(ipos);
      for (int l = 0; l < DIM; l++)
	{
	  a_k.m_particles[k].m_x[l] = UG(pt, l)*(1-sk[0])*(1 - sk[1]) + UG(pt + e0, l)*sk[0]*(1 - sk[1]) + UG(pt + e1, l)*sk[1]*(1 - sk[0]) + UG(pt + e0 + e1, l)*sk[0]*sk[1];
	}
    }
  a_k *= dt;
}
