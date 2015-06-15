module chem_rates
      implicit none
      contains

      subroutine Neut_Center(cx,cy,cz)
! Locates the cartesian coords of the center of the neutral cloud at a given time t.
            use dimensions
            use grid, only: ri,rj,rk,dz_grid,qx,qy,qz
            use inputs, only: dx,dy
            implicit none
            real, intent(out):: cx,cy,cz
            real:: x0,y0,z0
            
            x0 = dx/2.0
            y0 = dy/2.0
            z0 = dz_grid(nz/2)/2.0

            cx = qx(ri) + x0
            cy = qy(rj) + x0
            cz = qz(rj) + z0

      end subroutine Neut_Center
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Ionize_sw_mp(m_tstep)
! Ionizeds the neutral cloud with a 28s time constant and fill particle arrays, np, vp, up.
            use dimensions
            use MPI
            use misc
            use gutsp
            use grid, only: qx,qy,qz
            use inputs, only: PI,nf_init,omega_p,beta_pu,dx,dy,delz,dt,vsw,mion,km_to_m
            use var_arrays, only: input_p,up,Ni_tot,input_E,ijkp,m_arr,mrat,beta,beta_p,wght,np,vp,vp1,xp
            implicit none
            integer, intent(in):: m_tstep
            real:: ndot,dNr,theta
            integer:: i,j,k,l,m,dNi,l1,ierr
            
            ndot = 1e-6*nf_init*omega_p         !normalized ionization rate
            dNr = ndot*dt*beta*beta_pu*dx*delz*nz
            
            if (dNr .lt. 1.0) then
                  if (dNr .gt. pad_ranf()) then
                        dNr = 1.0
                        write(*,*) 'new ions...', dNr
                  endif
            endif
            
            dNi = nint(dNr)
            
            l1 = Ni_tot +1
            
            do l= l1, l1+dNi-1
                  theta = pad_ranf()*2*PI
                  
                  vp(l,1) = vsw+vsw*cos(theta) !+dvx
                  vp(l,2) = vsw*sin(theta) !+dvz 
                  vp(l,3) = 0.0

                  xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
                  xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
                  xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz-1)-qz(1))
                  
                  i=1
                  do while (xp(l,1) .gt. qx(i))
                        i=1+1
                  enddo
                  i=i-1
                  ijkp(l,1) = i
                  
                  ijkp(l,2) = floor(xp(l,2)/dy)
                  
                  k=1
                  do while (xp(l,3) .gt. qz(k))
                        k=k+1
                  enddo
                  k=k-1
                  ijkp(l,3) = k
                  
                  mrat(l) = 1.0
                  m_arr(l) = mion
                  beta_p(l) = beta_pu
                  
                  do m=1,3
                        vp1(l,m) = vp(l,m)
                        input_E = input_E + 0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 / &
                              beta*beta_p(l)
                        input_p(m) = input_p(m) + m_arr(l)*vp(l,m) / beta*beta_p(l)
                  enddo
                  
            enddo
                  
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
                  
            Ni_tot = Ni_tot + dNi
                  
            call get_interp_weights()
            call update_np()
            call update_up(vp)
            
      end subroutine Ionize_sw_mp
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Mass_load_Io(m_tstep)
! Mass load Io with ions witha a specific weight interp
            use dimensions
            use MPI
            use misc
            use gutsp
            use mult_proc, only: procnum, my_rank
            use grid, only: qx,qy,qz,dz_grid
            use inputs, only: PI,dt,mion,nf_init,mu0,b0_init,dx,dy,beta_pu,m_pu,ion_amu,km_to_m, load_rate
            use var_arrays, only: input_p,up,Ni_tot,input_E,ijkp,m_arr,mrat,beta,beta_p,wght,np,vp,vp1,xp
            implicit none
            integer, intent(in):: m_tstep
            real:: ddni, theta2, rand1, mdot, sca1 !scaling parameter
            integer:: i,j,k,l,m,l1,dNi,flg,ierr
            
     !       if ((m_tstep .gt. 20) .and. (m_tstep .lt. 600)) then
                  ! default load_rate = 0.1
                  sca1 = exp(-(dt*m_tstep-2)**2/load_rate**2)
                  
                  mdot = sca1*sqrt(mion*nf_init/mu0/1e9)*b0_init*dx*dy*1e6
                  if (my_rank == 0) then
                  write(*,*) 'mdot.........', mdot
                  endif
                  
                 
                  
                  ddni = dt*mdot*beta*beta_pu/(procnum*1.67e-27*m_pu)
                  if (my_rank == 0) then
                  write(*,*) 'ddni......', ddni
                  endif
                  !stop
                  
                  if (ddni .lt. 1.0) then
                        if (ddni .gt. pad_ranf()) then
                              ddni = 1.0
                        else
                              ddni = 0.0
                        endif
                  endif
                  
                  dNi = nint(ddni)
                  !write(*,*) 'dNi per processor......', dNi
                  !stop
                  l1 = Ni_tot + 1
                  
                  do l= l1, l1+ dNi-1
                        beta_p(l) = beta_pu
                        m_arr(l) = m_pu*1.67e-27
                        mrat(l) = ion_amu/m_pu
!                        vp(l,1) = 0
!                        vp(l,2) = 0
!                        vp(l,3) = 0
                        theta2 = pad_ranf()*2*PI
                        vp(l,1) = 0.0!57.0*cos(theta2)
                        vp(l,2) = 0.0!57.0*sin(theta2)
                        vp(l,3) = 0.0
                        
                        xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
                        xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
                        
                        flg=0
                        do while (flg .eq. 0)
                              xp(l,3) = qz(nz/2-20) + (1.0-pad_ranf())*(qz(nz/2+20)-qz(nz/2-20))
                              rand1=pad_ranf()
                              if (exp(-(xp(l,3)-qz(nz/2))**2/(5*dz_grid(nz/2)**2)) .gt. rand1) then
                                    flg = 1
                              endif
                        enddo
                        
                        i=1
                        do while (xp(l,1) .gt. qx(i))
                              i=1+1
                        enddo
                        i=i-1
                        ijkp(l,1) = i
                  
                        ijkp(l,2) = floor(xp(l,2)/dy)
                  
                        k=1
                        do while (xp(l,3) .gt. qz(k))
                              k=k+1
                        enddo
                        k=k-1
                        ijkp(l,3) = k
                  
                        do m=1,3
                              vp1(l,m) = vp(l,m)
                              input_E = input_E + 0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 / &
                                    beta*beta_p(l)
                              input_p(m) = input_p(m) + m_arr(l)*vp(l,m)/ &
                                    beta*beta_p(l)
                        enddo
                        
                  enddo      
                 
                  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
                  Ni_tot = Ni_tot + dNi
                 
                  call get_interp_weights()
                  call update_np()
                  call update_up(vp)
                 
      !      endif
            
      end subroutine Mass_load_Io
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine res_chex()
            use dimensions
            use misc
            use inputs, only: dt,m_pu
            use var_arrays, only: Ni_tot,m_arr,mrat,vp
            implicit none
            real, parameter:: sigma_chex = 8e-26
            real:: vrel,nn,chex_tau,chex_prob
            integer:: l
            
!            call Neut_Center(cx,cy,cz)
            
            do l=1, Ni_tot
                  vrel = sqrt(vp(l,1)**2+vp(l,2)**2+vp(l,3)**2)
                  nn = 10000e15
                  chex_tau = 1.0/(nn*sigma_chex*vrel)
                  chex_prob = dt/chex_tau
                  
                  if (pad_ranf() .lt. chex_prob) then
                        write(*,*) 'chex...',l,chex_tau,chex_prob
                        vp(l,:) = 0.0
                        mrat(l) = 16.0/m_pu
                        m_arr(l) = 1.67e-27*m_pu
                        
                  endif
            enddo
            
      end subroutine res_chex
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module chem_rates