module inputs
      use dimensions
      use mpi
!      use var_arrays, only: Ni_tot_0
      implicit none
      save
      
      real:: b0_init, nf_init,dt_frac, vsw, vth, Ni_tot_frac, dx_frac, &
            nu_init_frac,lambda_i,m_pu, mO, ppc, nu_init, ion_amu, load_rate, amp, &
            height_stretch, zsf, etemp0, mion, va
      real, parameter:: amu=1.6605e-27!, mion = 3.841e-26
      integer:: mp, nt, nout, loc, grad, nrgrd, boundx
      integer(4):: Ni_tot_0

      real, parameter:: q=1.6e-19         !electron charge

!       Grid paramters
      
      real:: dx,dy,delz,dt,dtsub_init

!       time stepping paramters

      integer, parameter:: ntsub = 10    !number of subcycle timesteps
      
!       Output directory

      character(50):: out_dir
      
!     logical variable for restart

      integer, parameter:: mrestart = 6000       !use -1 for no save
      
!       Neutral cloud expansion characteristics
      real:: vtop, vbottom
      
!       Max number of ions to be produced
!      integer(4):: Ni_tot_0
      
!       Misc. constants

      real, parameter:: pi = 3.14159
      real, parameter:: mu0 = pi*4.0e-7, &
                        epsilon = 8.85e-12, &
                        rtod=180.0/pi, &
                        km_to_m = 1.0e3, &
                        kboltz = 1.38e-29, &
                        tempf0 = 50*11600.0
                        
      real:: np_top, np_bottom, b0_top, b0_bottom, Lo, vth_top, vth_bottom, vth_max, &
             m_top, m_bottom, m_heavy, np_bottom_proton, f_proton_top, mBa
             
             
      real, parameter:: beta_particle = 1.0     !beta value of intial particles       
      real, parameter:: beta_pu = 10.0            !beta value of pickup ions
      real:: omega_p                            !ion gyrofrequency
      
!       Electron ion collision frequency
      real, parameter:: lww2 = 1.00             !artificial diffusion for the magnetic field update
      real, parameter:: lww1 = (1-lww2)/6.0     !divide by six for nearest neighbor
      
!       Density scaling paramter, alpha, and ion particle array dims

      real:: alpha
      
      
      contains
      
            subroutine readInputs()
                 implicit none
                  
                 open(unit=100, file= 'inputs.dat', status='old')
                  
                 read(100,*) b0_init
                 write(*,*) 'b0_init...........',b0_init
                 read(100,*) ion_amu
                 write(*,*) 'amu...............',ion_amu
                 read(100,*) m_pu
                 write(*,*) 'm_pu..............',m_pu
                 read(100,*) nf_init
                 write(*,*) 'nf_init...........',nf_init
                 read(100,*) amp
                 write(*,*) 'amplitude.........', amp
                 read(100,*) dt_frac
                 write(*,*) 'dt_frac...........',dt_frac
                 read(100,*) nt
                 write(*,*) 'nt................',nt
                 read(100,*) nout
                 write(*,*) 'nout..............',nout
                 read(100,*) vsw
                 write(*,*) 'vsw...............',vsw
                 read(100,*) vth
                 write(*,*) 'vth...............',vth
                 read(100,*) Ni_tot_frac
                 write(*,*) 'Ni_tot_frac.......',Ni_tot_frac
                 read(100,*) dx_frac
                 write(*,*) 'dx_frac...........',dx_frac
                 read(100,*) grad
                 write(*,*) 'scale height......', grad
                 read(100,*) height_stretch
                 write(*,*) 'start stretching...', height_stretch
                 read(100,*) zsf
                 write(*,*) 'z stretch factor...', zsf
                 read(100,*) nu_init_frac
                 write(*,*) 'nu_init_frac......',nu_init_frac
                 read(100,*) ppc
                 write(*,*) 'part per cell.....',ppc
                 read(100,*) loc
                 write(*,*) 'location of transform...',loc
                 read(100,*) load_rate
                 write(*,*) 'mass loading rate.....', load_rate
                 read(100,*) etemp0
                 write(*,*) 'electon temperature (eV)...', etemp0
                 read(100,*) boundx
                 write(*,*) 'boundary condition......', boundx
                 read(100,*) out_dir
                 write(*,*) 'output dir........',out_dir
                 
                 close(100)
                 
                
            end subroutine readInputs
            
            
            
            subroutine initparameters()
                  implicit none
                  
                  
                  
                  mion = amu*ion_amu!3.841e-26
                  write(*,*) 'mion...',mion
                  
                  omega_p = q*b0_init/mion
                  
                  lambda_i = (3e8/sqrt((nf_init*amp/1e9)*q*q/(8.85e-12*mion)))/1e3
                                    
                  dx= lambda_i*dx_frac
                  dy=lambda_i*dx_frac           !units in km
                  delz = lambda_i*dx_frac       !dz at release coordinates
                  
                  dt= dt_frac*mion/(q*b0_init)  !main time step
                  dtsub_init = dt/ntsub         !subcycle time step
                  vtop = vsw
                  vbottom = -vsw
                  
                  Ni_tot_0 = int(Ni_max*Ni_tot_frac)
                  write(*,*) 'Ni_tot_0...',Ni_tot_0, Ni_max, Ni_tot_frac
                  
                  mO = mion     !mass of H (kg)
                  
                  mBa = m_pu*mO !mass of Ba (kg)
                  
                  m_heavy = 1.0
                  np_top = nf_init
                  np_bottom = nf_init/m_heavy
                  f_proton_top = 0.5
                  b0_top= 1.0*b0_init
                  b0_bottom = b0_init
                  vth_top = vth
                  vth_bottom = vth
                  vth_max = 3*vth
                  m_top = mion
                  m_bottom = mion
                  Lo = 4.0*dx           !gradient scale length of boundary
                  
                  nu_init = nu_init_frac*omega_p
                  
                  alpha = (mu0/1.0e3)*q*(q/mion)  !mH...determines particle scaling
                  
                  nrgrd = nint(grad*height_stretch)
                  
            end subroutine initparameters
            
            
            subroutine check_inputs()
                  use mult_proc, only: my_rank
                  use grid, only: dz_grid
                  use var_arrays, only: np
                  implicit none
                  integer:: i,j,k
                  
                  real*8:: ak, btot, a1, a2, womega, phi, deltat, cwpi
                  
      ! Check input paramters
      
                  if (my_rank .eq. 0) then
                        write(*,*) mu0,q,ion_amu,mion
                        write(*,*) 'alpha...', alpha
                        write(*,*) 'c/wpi...', lambda_i,dx,dy,delz
                        write(*,*) 'dt......', dt
                        write(*,*) 'dt_sub...', dtsub_init
                        
                        
                        do i=1,nx
                        do j=1,ny
                        do k=1,nz
                        
                        ak = 2.0/dz_grid(k)
                        btot = b0_init*q/mion
                        a1 = ak**2*btot/(alpha*np(i,j,k))
                        a2 = (ak*btot)**2/(alpha*np(i,j,k))
                        womega = 0.5*(a1+ sqrt(a1**2 + 4.0*a2))
                        phi = womega/ak
                        deltat = dz_grid(k)/phi
                        if (deltat/dtsub_init .le. 2.0) then
                              write(*,*) 'deltat/dtsub....', deltat/dtsub_init
                              write(*,*) 'Field time stp too long...', i,j,k
!                              stop
                        endif
                        enddo
                        enddo
                        enddo
                        
                        write(*,*) 'Courant check (>1)...', deltat/dtsub_init
                       ! write(*,*) 'deltat, dtsub_init...', deltat, dtsub_init
                      
                        
                        write(*,*) '  '
                        write(*,*) 'Bottom paramters...'
                        write(*,*) '  '
                        va = b0_init/sqrt(mu0*m_bottom*np_bottom/1e9)/1e3
                        
                        write(*,*) 'Alfven veloctiy ......', va
                        write(*,*) 'Thermal velocity......', vth_top
                        write(*,*) 'Mach number...........', vbottom/(va+vth_bottom)
                        
                        write(*,*) 'Thermal gyroradius....', m_bottom*vth_bottom/(q*b0_init),m_bottom*vth_bottom/(q*b0_init)/dx
                        cwpi = 3.0e8/sqrt((np_bottom/1.0e9)*q*q/(epsilon*m_bottom))
                        write(*,*) 'Ion inertial length...', cwpi/1e3,cwpi/1e3/dx
                        
!                        write(*,*) 'Particles per cell...', Ni_tot_sys/(nx*nz)

                        write(*,*) '  '
                        write(*,*) 'Top parameters...'
                        write(*,*) '  '
                        
                        va = b0_init/sqrt(mu0*m_top*np_top/1e9)/1e3
                        
                        write(*,*) 'Alfven velocity......', va
                        write(*,*) 'Thermal velocity.....', vth_top
                        write(*,*) 'Mach number..........', vtop/(va+ vth_top)
                        
                        write(*,*) 'Thermal gyroradius...', m_top*vth_top/(q*b0_init), m_top*vth_top/(q*b0_init)/dx
                        cwpi = 3.0e8/sqrt((np_top/1e9)*q*q/(epsilon*m_top))
                        write(*,*) 'Ion inertial length..', cwpi/1e3,cwpi/1e3/dx
                        
                        
                  endif
                  
                  
            end subroutine check_inputs
            
end module inputs
                        
