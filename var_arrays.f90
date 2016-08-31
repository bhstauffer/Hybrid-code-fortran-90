module Var_Arrays
      use dimensions
      implicit none
      real::      b0(nx,ny,nz,3), &     !ambient mag field
                  b1(nx,ny,nz,3), &     !1st order mag field   
                  b12(nx,ny,nz,3), &    !b1 at previous time step
                  b1p2(nx,ny,nz,3), &   !temp b1 at time level m+1
                  bt(nx,ny,nz,3), &     !total mag field, mc covarient
                  btmf(nx,ny,nz,3), &   !main cell contravarient bt field
                  btc(nx,ny,nz,3), &    !btmf at cell center for particle move
                  np(nx,ny,nz), &       !particle ion density at level n, n+1/2
                  np3(nx,ny,nz,3), &    
                  vp(Ni_max,3), &       !particle velocity at t level n+1/2
                  vp1(Ni_max,3), &      !particle velocity at t level n
                  vplus(Ni_max,3), &    !v+ used in velocity update
                  vminus(Ni_max,3), &   !v- used in velocity update
                  up(nx,ny,nz,3), &     !particle flow at n, n+1/2
                  xp(Ni_max,3), &       !coordinates of ion particles
                  aj(nx,ny,nz,3), &     !curlB/(alpha*n)
                  nu(nx,ny,nz), &       !collision frequency
                  Ep(Ni_max,3), &       !Ion particle electric field
                  E(nx,ny,nz,3), &        !E field from electron mom. eq.
                  temp_p(nx,ny,nz), &   !temperature
                  mnp(nx,ny,nz), &      !mass density
                  beta, beta_p(Ni_max), &       !variable for particle scaling
                  mrat(Ni_max), &       !mass array for mulit-ion species
                  m_arr(Ni_max), &
                  np_t(nx,ny,nz), &
                  np_b(nx,ny,nz), &
                  up_t(nx,ny,nz,3), &
                  up_b(nx,nz,nz,3), &
                  input_p(3), &
                  input_E, input_Eb, bndry_Eflux, prev_Etot, &
                  grav(nx,ny,nz), &            !gravity term
                  gradP(nx,ny,nz,3)            !electron pressure gradient
      
      integer(4):: Ni_tot
      
      !Location (indices) of particles in the grid
      
      integer:: ijkp(Ni_max,3)
      logical:: in_bounds(Ni_max)
      real:: mix_ind(Ni_max)
      
      !Weight variables for trilinear interpolation
      
      real:: wght(Ni_max,8)
                  
                  
      integer:: np_t_flg(Ni_max), np_b_flg(Ni_max)
      
end module Var_Arrays
                  