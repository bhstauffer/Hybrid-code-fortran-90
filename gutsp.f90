module gutsp
      implicit none
      contains

      subroutine remove_ion(ion_l)
! Removes particles from the simulation that have gone out of bounds
            use dimensions
            use inputs, only: km_to_m
            use var_arrays, only: xp,vp,vp1,Ni_tot,input_E,ijkp,beta,beta_p,m_arr,mrat,wght
            implicit none
            integer, intent(in):: ion_l
            integer:: l,m
            
            do m=1,3    !remove ion energy from total input energy
                  input_E = input_E - 0.5*m_arr(l)*(vp(ion_l,m)*km_to_m)**2 &
                        / (beta * beta_p(l))
            enddo

            do l = ion_l, Ni_tot-1
                  do m = 1,3
                        xp(l,m) = xp(l+1,m)
                        vp(l,m) = vp(l+1,m)
                        vp(1,m) = vp(l+1,m)
                        ijkp(l,m) = ijkp(l+1,m)
                  enddo
                  beta_p(l) = beta_p(l+1)
                  m_arr(l) = m_arr(l+1)
                  mrat(l) = mrat(l+1)
            enddo

            do m=1,8
                  do l= ion_l,Ni_tot-1
                        wght(l,m) = wght(l+1,m)
                  enddo
            enddo

            Ni_tot = Ni_tot - 1

      end subroutine remove_ion

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine check_min_den()
            use dimensions
            use boundary
            use grid, only: dz_grid,qx,qy,qz
            use inputs, only: dx,dy,mion,beta_particle
            use var_arrays, only: np,xp,vp,up,Ni_tot,ijkp,beta,beta_p,m_arr,mrat
            implicit none
            real:: den_part, minden
            integer:: i,j,k,l,m,kk,ipart,npart,ii,jj
            
            den_part = 1.0/(beta*dx**3)
            minden=2.0*den_part
            
            do i=2,nx-1
                  do j=2,ny-1
                        do k=2,nz-1
!                              ak = PI/dx
!                              btot = sqrt(bt(i,j,k,1)**2 + bt(i,j,k,2)**2 + bt(i,j,k,3)**2)
!                              a1 = ak**2*btot/(alpha*np(i,j,k))
!                              a2 = (ak*btot)**2/(alpha*np(i,j,k))
!                              womega = 0.5*(a1 + sqrt(a1**2 + 4*a2)
!                              phi = womega/ak
!                              deltat = 0.1*dx/phi
                              
                              if (np(i,j,k) .le. minden) then
                                    npart = nint(minden/np(i,j,k))
                                    do ipart = 1, npart
                                          l = Ni_tot + 1
                                          do m=1,3
                                                vp(l,m) = up(i,j,k,m)
                                          enddo
                                          xp(l,1) = qx(i) + (0.5-pad_ranf())*dx
                                          xp(l,2) = qy(j) + (0.5-pad_ranf())*dy
                                          xp(l,3) = qz(k) + (0.5-pad_ranf())*dz_grid(k)
                                          
!                                          ijkp(l,1) = nint(xp(l,1)/dx)
!                                          ijkp(l,2) = nint(xp(l,2)/dy)
                                          call get_pindex(ii,jj,kk,l)
                                          kk = 1
                                          do while (xp(l,3) .gt. qz(kk))
                                                ijkp(l,3) = kk
                                                kk=kk+1
                                          enddo
                                          kk= ijkp(l,3)
                                          
                                          if (xp(l,3) .gt. (qz(kk) + (dz_grid(kk)/2))) then
                                                ijkp(l,3) = kk+1
                                          endif
                                          mrat(l) = 1.0
                                          m_arr(l) = mion
                                          beta_p(l) = beta_particle
                                          Ni_tot = Ni_tot +1
                                    enddo
                              endif
                        enddo
                  enddo
            enddo
            
      end subroutine check_min_den
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine extrapol_up()
! This subroutine does the provisional extrapolation of the particle
! bulk flow velocity to time level n, and replaces up_n-3/2 
! with up_n-1/2
            use dimensions
            use var_arrays, only: up,vp,vp1,np,Ni_tot,beta,beta_p,wght
            implicit none
            real:: v_at_n(Ni_max,3)
            integer:: l,m
            
            do l=1,Ni_tot
                  do m=1,3
                        v_at_n(l,m) = 1.5*vp(l,m) -0.5*vp1(l,m)
                  enddo
            enddo
            
            call update_up(v_at_n)
            
      end subroutine extrapol_up
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_Ep()
            use dimensions
            use grid_interp
            use var_arrays, only: Ep,aj,up,btc,Ni_tot,ijkp,mrat,wght,grav, gradP
            use inputs, only: mion
            implicit none
            real:: ajc(nx,ny,nz,3), &     !aj at cell center
!                   upc(nx,ny,nz,3), &   !up at cell center
                   gravc(nx,ny,nz), & !gravity at cell center
                   aa(3),bb(3),cc(3),aj3(3),up3(3),btc3(3), grav3, gradP3(3)    !dummy variables
            integer:: l,i,j,k,m,ip,jp,kp
            
            
            call face_to_center(aj,ajc)
!            call face_to_center(up,upc)


            ! grav term
            call grav_to_center(grav,gravc)
            
            do l=1, Ni_tot
                  i=ijkp(l,1)
                  j=ijkp(l,2)
                  k=ijkp(l,3)
                  
                  ip=i+1
                  jp=j+1
                  kp=k+1            
                  
                  do m=1,3
                        aj3(m) = ajc(i,j,k,m)*wght(l,1) + ajc(ip,j,k,m)*wght(l,2) &
                              + ajc(i,j,kp,m)*wght(l,3) + ajc(ip,j,kp,m)*wght(l,4) &
                              + ajc(i,jp,k,m)*wght(l,5) + ajc(ip,jp,k,m)*wght(l,6) &
                              + ajc(i,jp,kp,m)*wght(l,7) + ajc(ip,jp,kp,m)*wght(l,8)

                        up3(m) = up(i,j,k,m)*wght(l,1) + up(ip,j,k,m)*wght(l,2) &
                              + up(i,j,kp,m)*wght(l,3) + up(ip,j,kp,m)*wght(l,4) &
                              + up(i,jp,k,m)*wght(l,5) + up(ip,jp,k,m)*wght(l,6) &
                              + up(i,jp,kp,m)*wght(l,7) + up(ip,jp,kp,m)*wght(l,8)

!                        uf3(m) = ufc(i,j,k,m)*wght(l,1) + ufc(ip,j,k,m)*wght(l,2) 
!                              + ufc(i,j,kp,m)*wght(l,3) + ufc(ip,j,kp,m)*wght(l,4) &
!                              + ufc(i,jp,k,m)*wght(l,5) + ufc(ip,jp,k,m)*wght(l,6) &
!                              + ufc(i,jp,kp,m)*wght(l,7) + ufc(ip,jp,kp,m)*wght(l,8)

                        btc3(m) = btc(i,j,k,m)*wght(l,1) & 
                              + btc(ip,j,k,m)*wght(l,2) &
                              + btc(i,j,kp,m)*wght(l,3) &
                              + btc(ip,j,kp,m)*wght(l,4) &
                              + btc(i,jp,k,m)*wght(l,5) &
                              + btc(ip,jp,k,m)*wght(l,6) &
                              + btc(i,jp,kp,m)*wght(l,7) &
                              + btc(ip,jp,kp,m)*wght(l,8)
                               
                               
                       !electron pressure term
                       gradP3(m) = gradP(i,j,k,m)*wght(l,1) &
                               + gradP(ip,j,k,m)*wght(l,2) &
                               + gradP(i,j,kp,m)*wght(l,3) &
                               + gradP(ip,j,kp,m)*wght(l,4) &
                               + gradP(i,jp,k,m)*wght(l,5) &
                               + gradP(ip,jp,k,m)*wght(l,6) &
                               + gradP(i,jp,kp,m)*wght(l,7) &
                               + gradP(ip,jp,kp,m)*wght(l,8) 
                   
                  enddo 
                  do m=1,3 
                        aa(m) = aj3(m) - up3(m)
                        bb(m) = btc3(m)
                  enddo
                  ! Add in gravity term
                  grav3 = gravc(i,j,k)*wght(l,1) & 
                              + gravc(ip,j,k)*wght(l,2) &
                              + gravc(i,j,kp)*wght(l,3) &
                              + gravc(ip,j,kp)*wght(l,4) &
                              + gravc(i,jp,k)*wght(l,5) &
                              + gravc(ip,jp,k)*wght(l,6) &
                              + gravc(i,jp,kp)*wght(l,7) &
                              + gravc(ip,jp,kp)*wght(l,8)
                  
                  !Cross product
                  cc(1) = aa(2)*bb(3) - aa(3)*bb(2)
                  cc(2) = aa(3)*bb(1) - aa(1)*bb(3)
                  cc(3) = aa(1)*bb(2) - aa(2)*bb(1)
                  
                  
                  do m=1,2
                        Ep(l,m) = cc(m) - gradP3(m) !add in electron pressure term
                        Ep(l,m) = Ep(l,m) * mrat(l)
                  enddo
                  Ep(l,3) = cc(3) - gradP3(3) !add in electron pressure term
                  Ep(l,3) = Ep(l,3) * mrat(l) + grav3*mrat(l)  ! Second term is for gravity
!                  write(*,*) 'Electric field..............', Ep(l,m)*mrat(l)
!                  write(*,*) 'Gravity field...............', grav3*mrat(l), gravc(2,2,2), sum(wght(l,:))
!                  stop
                 

                  
            enddo
            !write(*,*) 'electric field, gravity....', maxval(Ep(:,:)), maxval(gravc(:,:,:))
            
      end subroutine get_Ep
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_vplus_vminus()
            use dimensions
            use inputs, only: dt
            use var_arrays, only: Ep,btc,vp,vplus,vminus,Ni_tot,ijkp,mrat,wght
            implicit none
            real:: a1,a2,a3,a_d,B2,dt2,Bx,By,Bz,vminus_x_B(3),vminus_dot_B,btc3(3)
            integer:: l,i,j,k,ip,jp,kp,m
            
            do l=1, Ni_tot
                  do m=1,3
                        vminus(l,m) = vp(l,m) + 0.5*dt*Ep(l,m)
                  enddo
            enddo
            
            do l = 1, Ni_tot
                  i=ijkp(l,1)
                  j=ijkp(l,2)
                  k=ijkp(l,3)
                  
                  ip=i+1
                  jp=j+1
                  kp=k+1
                  
                  do m=1,3
                        btc3(m) = btc(i,j,k,m)*wght(l,1) &
                              + btc(ip,j,k,m)*wght(l,2) & 
                              + btc(i,j,kp,m)*wght(l,3) &
                              + btc(ip,j,kp,m)*wght(l,4) &
                              + btc(i,jp,k,m)*wght(l,5) &
                              + btc(ip,jp,k,m)*wght(l,6) &
                              + btc(i,jp,kp,m)*wght(l,7) &
                              + btc(ip,jp,kp,m)*wght(l,8) 
                  enddo
                  
                  vminus_x_B(1) = vminus(l,2)*btc3(3)*mrat(l) - &
                        vminus(l,3)*btc3(2)*mrat(l)   
                  vminus_x_B(2) = vminus(l,3)*btc3(1)*mrat(l) - &
                        vminus(l,1)*btc3(3)*mrat(l)   
                  vminus_x_B(3) = vminus(l,1)*btc3(2)*mrat(l) - &
                        vminus(l,2)*btc3(1)*mrat(l)   

                  vminus_dot_B = vminus(l,1)*btc3(1)*mrat(l) + &
                        vminus(l,2)*btc3(2)*mrat(l) + &
                        vminus(l,3)*btc3(3)*mrat(l)   

                  Bx = btc3(1)*mrat(l) 
                  By = btc3(2)*mrat(l) 
                  Bz = btc3(3)*mrat(l) 
      
                  B2 = Bx*Bx + By*By + Bz*Bz
                  dt2 = dt*dt

                  a_d = 1.0/(1.0 + (0.25*B2*dt2))
                  a1 = (1 - (0.25*B2*dt2))*a_d
                  a2 = dt*a_d
                  a3 = 0.5*dt2*a_d
                  
                  do m=1,3
                        vplus(l,m) = a1*vminus(l,m) + a2*vminus_x_B(m) + &
                              a3*vminus_dot_B*btc3(m)*mrat(l)
                  enddo
            enddo
            
      end subroutine get_vplus_vminus
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine improve_up()
! The routine calculates v at time level n, and the associated bulk
! flow velocity up using the v+, v- technique.  The new up at
! time level n replaces the provisional extrapolation for up.
            use dimensions
            use var_arrays, only:vp1,vplus,vminus,up,np,Ni_tot,beta,beta_p,wght
            implicit none
            integer:: l,m
            
            do l=1, Ni_tot
                  do m=1,3
                        vp1(l,m) = 0.5*(vplus(l,m) + vminus(l,m))
                  enddo
            enddo
            
            call update_up(vp1)
            
      end subroutine improve_up
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_vp_final()
            use dimensions
            use inputs, only: dt
            use var_arrays, only: Ep,vp,vp1,vplus,Ni_tot
            implicit none

            integer:: l,m,ierr
            
            do l=1, Ni_tot
                  do m=1,3
                        vp1(l,m) = vp(l,m)      !to be used in extrapol_up for n-3/2
                        vp(l,m) = vplus(l,m) + 0.5*dt*Ep(l,m)
                  enddo
            enddo
            
      end subroutine get_vp_final
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine move_ion_half()
            use dimensions
            use boundary
            use inputs, only: dt,boundx
            use grid, only: qx,qy,qz
            use var_arrays, only: xp,vp,Ni_tot
            implicit none
            real:: dth
            integer:: l
            
            dth = dt/2.0
            
            if (boundx .eq. 1) then
            do l=1, Ni_tot
                  xp(l,1) = xp(l,1) + dth*vp(l,1)
                  xp(l,2) = xp(l,2) + dth*vp(l,2)
                  xp(l,3) = xp(l,3) + dth*vp(l,3)
                  
                  
                !  Periodic boundary conditions
                  
                        if (xp(l,1) .gt. qx(nx-1)) then
                              xp(l,1) = qx(1) + (xp(l,1) - qx(nx-1))
                        else if (xp(l,1) .le. qx(1)) then
                              xp(l,1) = qx(nx-1) -(qx(1)-xp(l,1))
                        endif
                        if (xp(l,2) .gt. qy(ny-1)) then
                              xp(l,2) = qy(1) + (xp(l,2) - qy(ny-1))
                        else if (xp(l,2) .le. qy(1)) then
                              xp(l,2) = qy(ny-1) -(qy(1)-xp(l,2))
                        endif
                        if (xp(l,3) .gt. qz(nz-1)) then
                              xp(l,3) = qz(1) + (xp(l,3) - qz(nz-1))
                        else if (xp(l,3) .le. qz(1)) then
                              xp(l,3) = qz(nz-1) -(qz(1)-xp(l,3))
                        endif
                  
            enddo
            else
                  do l=1,Ni_tot
                        xp(l,1) = xp(l,1) + dth*vp(l,1)
                        xp(l,2) = xp(l,2) + dth*vp(l,2)
                        xp(l,3) = xp(l,3) + dth*vp(l,3)
                  enddo      
                  
                  call particle_boundary()
            endif      
            
      end subroutine move_ion_half
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine check_min_den_boundary()
            use dimensions
            use grid, only: qx,qy,qz,dz_grid
            use mult_proc, only: my_rank, procnum
            use inputs, only: dx,np_top,vth_max,vsw,dx,dy,mion,m_top,km_to_m,Lo,vth_top,np_bottom,beta_particle
            use var_arrays, only: np,xp,vp,Ni_tot,input_E,ijkp,beta,beta_p,m_arr,mrat
            use boundary
            implicit none
            real:: den_part, minden,v,f,rnd,vx,vy,vz
            integer:: i,j,k,l,m,kk,flg,npart,ipart,ii,jj
            
            den_part = 1/(beta*dx**3)
            k = nx-1                    !top boundary
            minden = np_top-den_part
            
            do i=2,nx-1
                  do j=2,ny-1
                        if (np(i,j,k) .le. minden) then
                              npart = nint((np_top - np(i,j,k))/den_part)
                              if (my_rank .eq. nint(pad_ranf()*procnum)) then
                                    do ipart = 1,npart
                                          l = Ni_tot + 1        !beginning array element  for new borns
                                          
                                          flg = 0
                                          do while (flg .eq. 0)
                                                v = (2*vth_max*pad_ranf()) - vth_max
                                                f = exp(-(v**2)/vth_top**2)
                                                rnd = pad_ranf()
                                                if (f .ge. rnd) then
                                                      flg = 1
                                                      vx = v
                                                end if
                                          enddo
                                          flg = 0
                                          do while (flg .eq. 0)
                                                v = (2*vth_max*pad_ranf()) - vth_max
                                                f = exp(-(v**2)/vth_top**2)
                                                rnd = pad_ranf()
                                                if (f .ge. rnd) then
                                                      flg = 1
                                                      vy = v
                                                end if
                                          enddo
                                          flg = 0
                                          do while (flg .eq. 0)
                                                v = (2*vth_max*pad_ranf()) - vth_max
                                                f = exp(-(v**2)/vth_top**2)
                                                rnd = pad_ranf()
                                                if (f .ge. rnd) then
                                                      flg = 1
                                                      vz = v
                                                end if
                                          enddo
                                          
                                          vp(l,1) = vsw*(tanh((qz(nz)-qz(nz/2))/(Lo)))+vx 
                                          vp(l,2) = vy !*sin(phi)*sin(theta)
                                          vp(l,3) = vz !*cos(theta)
                     
                                          xp(l,1) = qx(i) + (0.5-pad_ranf())*dx
                                          xp(l,2) = qy(j) + (0.5-pad_ranf())*dy
                                          xp(l,3) = qz(k) + 2.0*(0.5-pad_ranf())*dz_grid(k)

                     
!                                          ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
!                                          ijkp(l,2) = nint(xp(l,2)/dy)
                                          call get_pindex(ii,jj,kk,l)
                                          kk = 1
                                          do while (xp(l,3) .gt. qz(kk))
                                                ijkp(l,3) = kk
                                                kk = kk+1
                                          enddo
                                          kk=ijkp(l,3)
                                          if (xp(l,3) .gt. (qz(kk) + (dz_grid(kk)/2))) then
                                                ijkp(l,3) = kk+1
                                          endif 
                                          mrat(l) = mion/m_top
                                          m_arr(l) = m_top
                                          beta_p(l) = beta_particle
                                          
                                          !Add energy
                                          do m=1,3
                                                input_E = input_E + &
                                                      0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 / (beta*beta_p(l))
                                          enddo
                                          
                                          Ni_tot = Ni_tot +1
                                    enddo
                              endif
                        endif
                  enddo
            enddo
            
            k=1 !bottom boundary
            minden = np_bottom-den_part
            
            do i=2,nx-1
                  do j=2,ny-1
                        if (np(i,j,k) .le. minden) then
                              npart = nint((np_bottom - np(i,j,k))/den_part)
                              if (my_rank .eq. nint(pad_ranf()*procnum)) then
                                    do ipart = 1,npart
                                          l = Ni_tot + 1        !beginning array element  for new borns
                                          
                                          flg = 0
                                          do while (flg .eq. 0)
                                                v = (2*vth_max*pad_ranf()) - vth_max
                                                f = exp(-(v**2)/vth_top**2)
                                                rnd = pad_ranf()
                                                if (f .ge. rnd) then
                                                      flg = 1
                                                      vx = v
                                                end if
                                          enddo
                                          flg = 0
                                          do while (flg .eq. 0)
                                                v = (2*vth_max*pad_ranf()) - vth_max
                                                f = exp(-(v**2)/vth_top**2)
                                                rnd = pad_ranf()
                                                if (f .ge. rnd) then
                                                      flg = 1
                                                      vy = v
                                                end if
                                          enddo
                                          flg = 0
                                          do while (flg .eq. 0)
                                                v = (2*vth_max*pad_ranf()) - vth_max
                                                f = exp(-(v**2)/vth_top**2)
                                                rnd = pad_ranf()
                                                if (f .ge. rnd) then
                                                      flg = 1
                                                      vz = v
                                                end if
                                          enddo
                                          
                                          vp(l,1) = vsw*(tanh((qz(nz)-qz(nz/2))/(Lo)))+vx 
                                          vp(l,2) = vy !*sin(phi)*sin(theta)
                                          vp(l,3) = vz !*cos(theta)
                     
                                          xp(l,1) = qx(i) + (0.5-pad_ranf())*dx
                                          xp(l,2) = qy(j) + (0.5-pad_ranf())*dy
                                          xp(l,3) = qz(k) + 2.0*(0.5-pad_ranf())*dz_grid(k)

                     
!                                          ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
!                                          ijkp(l,2) = nint(xp(l,2)/dy)
                                          call get_pindex(ii,jj,kk,l)
                                          kk = 1
                                          do while (xp(l,3) .gt. qz(kk))
                                                ijkp(l,3) = kk
                                                kk = kk+1
                                          enddo
                                          kk=ijkp(l,3)
                                          if (xp(l,3) .gt. (qz(kk) + (dz_grid(kk)/2))) then
                                                ijkp(l,3) = kk+1
                                          endif 
                                          mrat(l) = mion/m_top
                                          m_arr(l) = m_top
                                          beta_p(l) = beta_particle
                                          
                                          !Add energy
                                          do m=1,3
                                                input_E = input_E + &
                                                      0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 / (beta*beta_p(l))
                                          enddo
                                          
                                          Ni_tot = Ni_tot +1
                                    enddo
                              endif
                        endif
                  enddo
            enddo
            
      end subroutine check_min_den_boundary
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine check_min_den_boundary_1()
            use dimensions
            use grid, only: qx,qy,qz,dz_grid
            use mult_proc, only: my_rank, procnum
            use inputs, only: dx,np_top,vth_max,vsw,dx,dy,mion,m_top,km_to_m,Lo,vth_top,np_bottom,vth,beta_particle
            use boundary
            use var_arrays, only: np,xp,vp,Ni_tot,input_E,ijkp,beta,beta_p,m_arr,mrat
            implicit none
            real:: den_part, minden,v,f,rnd,vx,vy,vz
            integer:: i,j,k,l,m,kk,flg,npart,ipart,ii,jj
            
            den_part = 1/(beta*dx**3)
            k = nz-1    !top boundary
            minden = np_top-den_part
            
            do i=2,nx-1
                  do j= 2, ny-1
                        if (np(i,j,k) .le. minden) then
                              npart = nint((np_top - np(i,j,k))/den_part)
                              if (my_rank .eq. nint(pad_ranf()*procnum)) then
                                    do ipart = 1, npart
                                          l = Ni_tot + 1
                                          flg = 0
                                          do while (flg .eq. 0)
                                                vx = (600*pad_ranf()) - 300
                                                vy = (600*pad_ranf()) - 300
                                                vz = (600*pad_ranf()) - 300
                                                v = sqrt(vx**2+vy**2+vz**2)
                                                f = exp(-(v**2)/vth**2)
                                                rnd = pad_ranf()
                                                if (f .ge. rnd) then
                                                      flg = 1
                                                endif
                                          enddo
                                          
                                          vp(l,1) = vsw*(tanh((qz(nz)-qz(nz/2))/(Lo)))+vx 
                                          vp(l,2) = vy !*sin(phi)*sin(theta)
                                          vp(l,3) = vz !*cos(theta)
                     
                                          xp(l,1) = qx(i) + (0.5-pad_ranf())*dx
                                          xp(l,2) = qy(j) + (0.5-pad_ranf())*dy
                                          xp(l,3) = qz(k) + 2.0*(0.5-pad_ranf())*dz_grid(k)

                     
!                                          ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
!                                          ijkp(l,2) = nint(xp(l,2)/dy)
                                          call get_pindex(ii,jj,kk,l)
                                          kk=1
                                          do while (xp(l,3) .gt. qz(kk))
                                                ijkp(l,3) = kk
                                                kk = kk+1
                                          enddo
                                          kk=ijkp(l,3)
                                          if (xp(l,3) .gt. (qz(kk) + (dz_grid(kk)/2))) then
                                                ijkp(l,3) = kk+1
                                          endif 
                                          mrat(l) = mion/m_top
                                          m_arr(l) = m_top
                                          beta_p(l) = beta_particle 
                                          
                                         !Add energy
                                          do m=1,3
                                                input_E = input_E + &
                                                      0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 / (beta*beta_p(l))
                                          enddo
                                          
                                          Ni_tot = Ni_tot +1
                                    enddo
                              endif
                        endif
                  enddo
            enddo
            
            k=2 ! bottom boundary
            minden = np_bottom-den_part
            
            do i=2,nx-1
                  do j= 2, ny-1
                        if (np(i,j,k) .le. minden) then
                              npart = nint((np_bottom - np(i,j,k))/den_part)
                              if (my_rank .eq. nint(pad_ranf()*procnum)) then
                                    do ipart = 1, npart
                                          l = Ni_tot + 1
                                          flg = 0
                                          do while (flg .eq. 0)
                                                vx = (600*pad_ranf()) - 300
                                                vy = (600*pad_ranf()) - 300
                                                vz = (600*pad_ranf()) - 300
                                                v = sqrt(vx**2+vy**2+vz**2)
                                                f = exp(-(v**2)/vth**2)
                                                rnd = pad_ranf()
                                                if (f .ge. rnd) then
                                                      flg = 1
                                                endif
                                          enddo
                                          
                                          vp(l,1) = vsw*(tanh((qz(nz)-qz(nz/2))/(Lo)))+vx 
                                          vp(l,2) = vy !*sin(phi)*sin(theta)
                                          vp(l,3) = vz !*cos(theta)
                     
                                          xp(l,1) = qx(i) + (0.5-pad_ranf())*dx
                                          xp(l,2) = qy(j) + (0.5-pad_ranf())*dy
                                          xp(l,3) = qz(k) + 2.0*(0.5-pad_ranf())*dz_grid(k)

                     
!                                          ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
!                                          ijkp(l,2) = nint(xp(l,2)/dy)
                                          call get_pindex(ii,jj,kk,l)
                                          kk=1
                                          do while (xp(l,3) .gt. qz(kk))
                                                ijkp(l,3) = kk
                                                kk = kk+1
                                          enddo
                                          kk=ijkp(l,3)
                                          if (xp(l,3) .gt. (qz(kk) + (dz_grid(kk)/2))) then
                                                ijkp(l,3) = kk+1
                                          endif 
                                          mrat(l) = mion/m_top
                                          m_arr(l) = m_top
                                          beta_p(l) = beta_particle 
                                          
                                         !Add energy
                                          do m=1,3
                                                input_E = input_E + &
                                                      0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 / (beta*beta_p(l))
                                          enddo
                                          
                                          Ni_tot = Ni_tot +1
                                    enddo
                              endif
                        endif
                  enddo
            enddo         
            
      end subroutine check_min_den_boundary_1
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_interp_weights_2()
! Weights are used for trilinear interpolation to/from main cell
! centers to particle positions.  For each particle there are 8
! grid points associated with the interpolation.  These 8 points
! are determined by the location of the particle within the main
! cell.  There are 8 sets of 8 grid points for each cell.
            use dimensions
            use inputs, only: dx,dy
            use grid, only: qx,qy,qz
            use var_arrays, only: xp,Ni_tot,ijkp,wght
            implicit none
            real:: vol,x1,x2,y1,y2,z1,z2
            integer:: l,i,j,k
            
            do l=1, Ni_tot
                  i=ijkp(l,1)
                  j=ijkp(l,2)
                  k=ijkp(l,3)
                  
                  if ((xp(l,1) .le. qx(i)) .and. (xp(l,2) .le. qy(j)) .and. (xp(l,3) .le. qz(k))) then
                        vol = dx*dy*(qz(k)-qz(k-1))
                        x1=abs(xp(l,1)-qx(i))
                        x2=abs(xp(l,1)-qx(i+1))
                        y1=abs(xp(l,2)-qy(j-1))
                        y2=abs(xp(l,2)-qy(j))
                        z1=abs(xp(l,3)-qz(k-1))
                        z2=abs(xp(l,3)-qz(k))
                        wght(l,1) = x2*y2*z2/vol
                        wght(l,2) = x1*y2*z2/vol
                        wght(l,3) = x2*y2*z1/vol
                        wght(l,4) = x1*y2*z1/vol
                        wght(l,5) = x2*y1*z2/vol
                        wght(l,6) = x1*y1*z2/vol
                        wght(l,7) = x2*y1*z1/vol
                        wght(l,8) = x1*y1*z1/vol
                  endif
                  if ((xp(l,1) .gt. qx(i)) .and. (xp(l,2) .le. qy(j)) .and. (xp(l,3) .le. qz(k))) then
                        vol = dx*dy*(qz(k)-qz(k-1))
                        x1=abs(xp(l,1)-qx(i))
                        x2=abs(xp(l,1)-qx(i+1))
                        y1=abs(xp(l,2)-qy(j-1))
                        y2=abs(xp(l,2)-qy(j))
                        z1=abs(xp(l,3)-qz(k-1))
                        z2=abs(xp(l,3)-qz(k))
                        wght(l,1) = x2*y2*z2/vol
                        wght(l,2) = x1*y2*z2/vol
                        wght(l,3) = x2*y2*z1/vol
                        wght(l,4) = x1*y2*z1/vol
                        wght(l,5) = x2*y1*z2/vol
                        wght(l,6) = x1*y1*z2/vol
                        wght(l,7) = x2*y1*z1/vol
                        wght(l,8) = x1*y1*z1/vol
                  endif
                  if ((xp(l,1) .le. qx(i)) .and. (xp(l,2) .le. qy(j)) .and. (xp(l,3) .gt. qz(k))) then
                        vol = dx*dy*(qz(k+1)-qz(k))
                        x1=abs(xp(l,1)-qx(i-1))
                        x2=abs(xp(l,1)-qx(i))
                        y1=abs(xp(l,2)-qy(j-1))
                        y2=abs(xp(l,2)-qy(j))
                        z1=abs(xp(l,3)-qz(k))
                        z2=abs(xp(l,3)-qz(k+1))
                        wght(l,1) = x2*y2*z2/vol
                        wght(l,2) = x1*y2*z2/vol
                        wght(l,3) = x2*y2*z1/vol
                        wght(l,4) = x1*y2*z1/vol
                        wght(l,5) = x2*y1*z2/vol
                        wght(l,6) = x1*y1*z2/vol
                        wght(l,7) = x2*y1*z1/vol
                        wght(l,8) = x1*y1*z1/vol
                  endif
                  if ((xp(l,1) .gt. qx(i)) .and. (xp(l,2) .le. qy(j)) .and. (xp(l,3) .gt. qz(k))) then
                        vol = dx*dy*(qz(k+1)-qz(k))
                        x1=abs(xp(l,1)-qx(i))
                        x2=abs(xp(l,1)-qx(i+1))
                        y1=abs(xp(l,2)-qy(j-1))
                        y2=abs(xp(l,2)-qy(j))
                        z1=abs(xp(l,3)-qz(k))
                        z2=abs(xp(l,3)-qz(k+1))
                        wght(l,1) = x2*y2*z2/vol
                        wght(l,2) = x1*y2*z2/vol
                        wght(l,3) = x2*y2*z1/vol
                        wght(l,4) = x1*y2*z1/vol
                        wght(l,5) = x2*y1*z2/vol
                        wght(l,6) = x1*y1*z2/vol
                        wght(l,7) = x2*y1*z1/vol
                        wght(l,8) = x1*y1*z1/vol
                  endif
                  if ((xp(l,1) .le. qx(i)) .and. (xp(l,2) .gt. qy(j)) .and. (xp(l,3) .le. qz(k))) then
                        vol = dx*dy*(qz(k)-qz(k-1))
                        x1=abs(xp(l,1)-qx(i-1))
                        x2=abs(xp(l,1)-qx(i))
                        y1=abs(xp(l,2)-qy(j))
                        y2=abs(xp(l,2)-qy(j+1))
                        z1=abs(xp(l,3)-qz(k-1))
                        z2=abs(xp(l,3)-qz(k))
                        wght(l,1) = x2*y2*z2/vol
                        wght(l,2) = x1*y2*z2/vol
                        wght(l,3) = x2*y2*z1/vol
                        wght(l,4) = x1*y2*z1/vol
                        wght(l,5) = x2*y1*z2/vol
                        wght(l,6) = x1*y1*z2/vol
                        wght(l,7) = x2*y1*z1/vol
                        wght(l,8) = x1*y1*z1/vol
                  endif
                  if ((xp(l,1) .gt. qx(i)) .and. (xp(l,2) .gt. qy(j)) .and. (xp(l,3) .le. qz(k))) then
                        vol = dx*dy*(qz(k)-qz(k-1))
                        x1=abs(xp(l,1)-qx(i))
                        x2=abs(xp(l,1)-qx(i+1))
                        y1=abs(xp(l,2)-qy(j))
                        y2=abs(xp(l,2)-qy(j+1))
                        z1=abs(xp(l,3)-qz(k-1))
                        z2=abs(xp(l,3)-qz(k))
                        wght(l,1) = x2*y2*z2/vol
                        wght(l,2) = x1*y2*z2/vol
                        wght(l,3) = x2*y2*z1/vol
                        wght(l,4) = x1*y2*z1/vol
                        wght(l,5) = x2*y1*z2/vol
                        wght(l,6) = x1*y1*z2/vol
                        wght(l,7) = x2*y1*z1/vol
                        wght(l,8) = x1*y1*z1/vol
                  endif
                  if ((xp(l,1) .le. qx(i)) .and. (xp(l,2) .gt. qy(j)) .and. (xp(l,3) .gt. qz(k))) then
                        vol = dx*dy*(qz(k+1)-qz(k))
                        x1=abs(xp(l,1)-qx(i-1))
                        x2=abs(xp(l,1)-qx(i))
                        y1=abs(xp(l,2)-qy(j))
                        y2=abs(xp(l,2)-qy(j+1))
                        z1=abs(xp(l,3)-qz(k))
                        z2=abs(xp(l,3)-qz(k+1))
                        wght(l,1) = x2*y2*z2/vol
                        wght(l,2) = x1*y2*z2/vol
                        wght(l,3) = x2*y2*z1/vol
                        wght(l,4) = x1*y2*z1/vol
                        wght(l,5) = x2*y1*z2/vol
                        wght(l,6) = x1*y1*z2/vol
                        wght(l,7) = x2*y1*z1/vol
                        wght(l,8) = x1*y1*z1/vol
                  endif
                  if ((xp(l,1) .gt. qx(i)) .and. (xp(l,2) .gt. qy(j)) .and. (xp(l,3) .gt. qz(k))) then
                        vol = dx*dy*(qz(k+1)-qz(k))
                        x1=abs(xp(l,1)-qx(i))
                        x2=abs(xp(l,1)-qx(i+1))
                        y1=abs(xp(l,2)-qy(j))
                        y2=abs(xp(l,2)-qy(j+1))
                        z1=abs(xp(l,3)-qz(k))
                        z2=abs(xp(l,3)-qz(k+1))
                        wght(l,1) = x2*y2*z2/vol
                        wght(l,2) = x1*y2*z2/vol
                        wght(l,3) = x2*y2*z1/vol
                        wght(l,4) = x1*y2*z1/vol
                        wght(l,5) = x2*y1*z2/vol
                        wght(l,6) = x1*y1*z2/vol
                        wght(l,7) = x2*y1*z1/vol
                        wght(l,8) = x1*y1*z1/vol
                  endif
            enddo
            
      end subroutine get_interp_weights_2
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_interp_weights()
! Weights are used for trilinear interpolation to/from main cell
! centers to particle positions.  For each particle there are 8
! grid points associated with the interpolation.  These 8 points
! are determined by the location of the particle within the main
! cell.  There are 8 sets of 8 grid points for each cell.
            use dimensions
            use grid, only: qx,qy,qz
            use inputs, only: dx,dy
            use var_arrays, only: xp,Ni_tot,ijkp,wght
            implicit none
            real:: vol,x1,x2,y1,y2,z1,z2
            integer:: i,j,k,l
            
            do l=1, Ni_tot
!                  i=1

!                  do while (xp(l,1) .gt. qx(i))
!                        i=i+1
!                  enddo

!                  i=i-1
!                  ijkp(l,1)=i
!                  j = floor(xp(l,2)/dy)
!                  ijkp(l,2) = j
!                  k=1
!                  do while (xp(l,3) .gt. qz(k))
!                        k=k+1
!                  enddo
!                  k=k-1
!                  ijkp(l,3) = k
                  call get_pindex(i,j,k,l)
                 
                  
                  vol = 1.0/((qx(i+1)-qx(i))*(qy(j+1)-qy(j))*(qz(k+1)-qz(k)))
                  x1=abs(xp(l,1)-qx(i))
                  x2=abs(xp(l,1)-qx(i+1))
                  y1=abs(xp(l,2)-qy(j))
                  y2=abs(xp(l,2)-qy(j+1))
                  z1=abs(xp(l,3)-qz(k))
                  z2=abs(xp(l,3)-qz(k+1))
                 
                  wght(l,1) = x2*y2*z2*vol
                  wght(l,2) = x1*y2*z2*vol
                  wght(l,3) = x2*y2*z1*vol
                  wght(l,4) = x1*y2*z1*vol
                  wght(l,5) = x2*y1*z2*vol
                  wght(l,6) = x1*y1*z2*vol
                  wght(l,7) = x2*y1*z1*vol
                  wght(l,8) = x1*y1*z1*vol
                 
                
            enddo
            
            
      end subroutine get_interp_weights
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine update_np()
! Weight density to eight nearest grid points
            use dimensions
            use MPI
            use grid, only: dx_grid,dy_grid,dz_grid
            use boundary
            use var_arrays, only: np,Ni_tot,ijkp,beta,beta_p,wght
            implicit none
            real:: volb, recvbuf(nx*ny*nz)
            integer:: i,j,k,l,ip,jp,kp,ierr,count
            
            count = nx*ny*nz
            
            do i=1,nx
                  do j=1,ny
                        do k=1,nz
                              np(i,j,k) = 0.0
                        enddo
                  enddo
            enddo
            do l=1, Ni_tot
                  i=ijkp(l,1)
                  j=ijkp(l,2)
                  k=ijkp(l,3)
                  
                  ip=i+1
                  jp=j+1
                  kp=k+1
                  
                  volb = 1.0/(dx_grid(i)*dy_grid(j)*dz_grid(k)*beta*beta_p(l))
!                  volb = beta_p(l)/(dx_grid(i)*dy_grid(j)*dz_grid(k)*beta)
                  
                  np(i,j,k) = np(i,j,k) + wght(l,1)*volb
                  np(ip,j,k) = np(ip,j,k) + wght(l,2)*volb
                  np(i,j,kp) = np(i,j,kp) + wght(l,3)*volb
                  np(ip,j,kp) = np(ip,j,kp) + wght(l,4)*volb
                  np(i,jp,k) = np(i,jp,k) + wght(l,5)*volb
                  np(ip,jp,k) = np(ip,jp,k) + wght(l,6)*volb
                  np(i,jp,kp) = np(i,jp,kp) + wght(l,7)*volb
                  np(ip,jp,kp) = np(ip,jp,kp) + wght(l,8)*volb
                  
                  
            enddo
                  
            !Used for periodic boundary conditions
            call add_boundary_scalar(np)
            
!            np(nx-1,:,:) = np(nx-1,:,:) + np(1,:,:)
!            np(:,ny-1,:) = np(:,ny-1,:) + np(:,1,:)
!            np(:,:,nz-1) = np(:,:,nz-1) + np(:,:,1)
            
            
            
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(np(:,:,:),recvbuf,count,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
            
            np(:,:,:) = reshape(recvbuf,(/nx,ny,nz/))
            
            call boundary_scalar(np)
!            call periodic_scalar(np)
            
      end subroutine update_np
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      subroutine update_rho()
! Weight density to eight neares grid points
            use dimensions
            use MPI
            use boundary
            use inputs, only: mion
            use grid, only: qx,qy,qz
            use var_arrays, only: mnp,Ni_tot,ijkp,beta,beta_p,mrat,wght
            implicit none
            real:: volb, recvbuf(nx*ny*nz)
            integer:: i,j,k,l,ip,jp,kp,count,ierr
            
            count = nx*ny*nz
            
            do i=1,nx
                  do j=1,ny
                        do k=1,nz
                              mnp(i,j,k) = 0.0
                        enddo
                  enddo
            enddo
            
            do l = 1, Ni_tot
                  i=ijkp(l,1)
                  j=ijkp(l,2)
                  k=ijkp(l,3)
                  
                  ip=i+1
                  jp=j+1
                  kp=k+1
                  
                  volb = (qx(ip)-qx(i))*(qy(jp)-qy(j))*(qz(kp)-qz(k))*beta*beta_p(l)
                  
                  mnp(i,j,k) = mnp(i,j,k) + (wght(l,1)/mrat(l))/volb
                  mnp(ip,j,k) = mnp(ip,j,k) + (wght(l,2)/mrat(l))/volb
                  mnp(i,j,kp) = mnp(i,j,kp) + (wght(l,3)/mrat(l))/volb
                  mnp(ip,j,kp) = mnp(ip,j,kp) + (wght(l,4)/mrat(l))/volb
                  mnp(i,jp,k) = mnp(i,jp,k) + (wght(l,5)/mrat(l))/volb
                  mnp(ip,jp,k) = mnp(ip,jp,k) + (wght(l,6)/mrat(l))/volb
                  mnp(i,jp,kp) = mnp(i,jp,kp) + (wght(l,7)/mrat(l))/volb
                  mnp(ip,jp,kp) = mnp(ip,jp,kp) + (wght(l,8)/mrat(l))/volb
                  
            enddo
            
            mnp(:,:,:) = mion*mnp(:,:,:) !mass density
            
            !Used for periodic boundary conditions
            call add_boundary_scalar(mnp)
            
!            mnp(nx-1,:,:) = mnp(nx-1,:,:)+mnp(1,:,:)
!            mnp(:,ny-1,:) = mnp(:,ny-1,:)+mnp(:,1,:)
!            mnp(:,:,nz-1) = mnp(:,:,nz-1)+mnp(:,:,1)
            
          
            
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(mnp(:,:,:),recvbuf,count,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
            
            mnp(:,:,:) = reshape(recvbuf,(/nx,ny,nz/))
            call boundary_scalar(mnp)
!            call periodic_scalar(mnp)
            
      end subroutine update_rho
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine update_mixed(mixed,mix_cnt)
! Weight density to eight nearest grid points
            use dimensions
            use MPI
            use var_arrays, only: Ni_tot,ijkp,mix_ind
            implicit none
            real, intent(out):: mixed(nx,ny,nz), mix_cnt(nx,ny,nz)
            real:: recvbuf(nx*ny*nz)
            integer:: i,j,k,l,count,ierr
            
            count = nx*ny*nz
            
            do i=1,nx
                  do j=1,ny
                        do k=1,nz
                              mixed(i,j,k) = 0.0
                              mix_cnt(i,j,k) = 0.0
                        enddo
                  enddo
            enddo
            
            do l = 1, Ni_tot
                  i=ijkp(l,1)
                  j=ijkp(l,2)
                  k=ijkp(l,3)
                  
                  mixed(i,j,k) = mixed(i,j,k) + mix_ind(l)
                  mix_cnt(i,j,k) = mix_cnt(i,j,k) + 1
                  
            enddo
            
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(mixed(:,:,:),recvbuf,count,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
            
            mixed(:,:,:) = reshape(recvbuf,(/nx,ny,nz/))
            
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(mix_cnt(:,:,:),recvbuf,count,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
            
            mix_cnt(:,:,:) = reshape(recvbuf,(/nx,ny,nz/))
            
            mixed(:,:,:) = mixed(:,:,:)/mix_cnt(:,:,:)
            
      end subroutine update_mixed
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine update_up(vp)
            use dimensions
            use MPI
            use boundary
            use grid, only: qx,qy,qz
            use var_arrays, only: np,up,Ni_tot,ijkp,beta,beta_p,wght
            implicit none
            real, intent(in):: vp(Ni_max,3)
            real:: ct(nx,ny,nz,3), recvbuf(nx*ny*nz*3),volb,nvolb
            integer:: i,j,k,m,l,ip,jp,kp,count,ierr
            
            count=nx*ny*nz*3
            
            do i=1,nx
                  do j=1,ny
                        do k=1,nz
                              do m=1,3
                                    up(i,j,k,m) = 0.0
                                    ct(i,j,k,m) = 0.0
                              enddo
                        enddo
                  enddo
            enddo
            
            do l=1,Ni_tot
                  i=ijkp(l,1)
                  j=ijkp(l,2)
                  k=ijkp(l,3)
                  
                  ip=i+1
                  jp=j+1
                  kp=k+1
                  
                  volb = (qx(ip)-qx(i))*(qy(jp)-qy(j))*(qz(kp)-qz(k))*beta*beta_p(l)
                  
!                  if(np(i,j,k) .gt. 0.0) then                  
                  nvolb = 1.0/(np(i,j,k)*volb)
                  ct(i,j,k,1) = ct(i,j,k,1) + vp(l,1)*wght(l,1)*nvolb
                  ct(i,j,k,2) = ct(i,j,k,2) + vp(l,2)*wght(l,1)*nvolb
                  ct(i,j,k,3) = ct(i,j,k,3) + vp(l,3)*wght(l,1)*nvolb
!                  endif

!                  if (np(ip,j,k) .gt. 0.0) then
                  nvolb = 1.0/(np(ip,j,k)*volb)
                  ct(ip,j,k,1) = ct(ip,j,k,1) + vp(l,1)*wght(l,2)*nvolb
                  ct(ip,j,k,2) = ct(ip,j,k,2) + vp(l,2)*wght(l,2)*nvolb
                  ct(ip,j,k,3) = ct(ip,j,k,3) + vp(l,3)*wght(l,2)*nvolb
!                  endif

!                  if (np(i,j,kp) .gt. 0.0) then
                  nvolb = 1.0/(np(i,j,kp)*volb)
                  ct(i,j,kp,1) = ct(i,j,kp,1) + vp(l,1)*wght(l,3)*nvolb
                  ct(i,j,kp,2) = ct(i,j,kp,2) + vp(l,2)*wght(l,3)*nvolb
                  ct(i,j,kp,3) = ct(i,j,kp,3) + vp(l,3)*wght(l,3)*nvolb
!                  endif

!                  if (np(ip,j,kp) .gt. 0.0) then
                  nvolb = 1.0/(np(ip,j,kp)*volb)
                  ct(ip,j,kp,1) = ct(ip,j,kp,1) + vp(l,1)*wght(l,4)*nvolb
                  ct(ip,j,kp,2) = ct(ip,j,kp,2) + vp(l,2)*wght(l,4)*nvolb
                  ct(ip,j,kp,3) = ct(ip,j,kp,3) + vp(l,3)*wght(l,4)*nvolb
!                  endif

!                  if (np(i,jp,k) .gt. 0.0) then
                  nvolb = 1.0/(np(i,jp,k)*volb)
                  ct(i,jp,k,1) = ct(i,jp,k,1) + vp(l,1)*wght(l,5)*nvolb
                  ct(i,jp,k,2) = ct(i,jp,k,2) + vp(l,2)*wght(l,5)*nvolb
                  ct(i,jp,k,3) = ct(i,jp,k,3) + vp(l,3)*wght(l,5)*nvolb
!                  endif

!                  if (np(ip,jp,k) .gt. 0.0) then
                  nvolb = 1.0/(np(ip,jp,k)*volb)
                  ct(ip,jp,k,1) = ct(ip,jp,k,1) + vp(l,1)*wght(l,6)*nvolb
                  ct(ip,jp,k,2) = ct(ip,jp,k,2) + vp(l,2)*wght(l,6)*nvolb
                  ct(ip,jp,k,3) = ct(ip,jp,k,3) + vp(l,3)*wght(l,6)*nvolb
!                  endif

!                  if (np(i,jp,kp) .gt. 0.0) then
                  nvolb = 1.0/(np(i,jp,kp)*volb)
                  ct(i,jp,kp,1) = ct(i,jp,kp,1) + vp(l,1)*wght(l,7)*nvolb
                  ct(i,jp,kp,2) = ct(i,jp,kp,2) + vp(l,2)*wght(l,7)*nvolb
                  ct(i,jp,kp,3) = ct(i,jp,kp,3) + vp(l,3)*wght(l,7)*nvolb
!                  endif

!                  if (np(ip,jp,kp) .gt. 0.0) then
                  nvolb = 1.0/(np(ip,jp,kp)*volb)
                  ct(ip,jp,kp,1) = ct(ip,jp,kp,1) + vp(l,1)*wght(l,8)*nvolb
                  ct(ip,jp,kp,2) = ct(ip,jp,kp,2) + vp(l,2)*wght(l,8)*nvolb
                  ct(ip,jp,kp,3) = ct(ip,jp,kp,3) + vp(l,3)*wght(l,8)*nvolb
!                  endif


!         nvolb = np(i,j,k)*volb
!         ct(i,j,k,1) = ct(i,j,k,1) + vp(l,1)*wght(l,1)/nvolb
!         ct(i,j,k,2) = ct(i,j,k,2) + vp(l,2)*wght(l,1)/nvolb
!         ct(i,j,k,3) = ct(i,j,k,3) + vp(l,3)*wght(l,1)/nvolb

!         nvolb = np(ip,j,k)*volb
!         ct(ip,j,k,1) = ct(ip,j,k,1) + vp(l,1)*wght(l,2)/nvolb
!         ct(ip,j,k,2) = ct(ip,j,k,2) + vp(l,2)*wght(l,2)/nvolb
!         ct(ip,j,k,3) = ct(ip,j,k,3) + vp(l,3)*wght(l,2)/nvolb

!         nvolb = np(i,j,kp)*volb
!         ct(i,j,kp,1) = ct(i,j,kp,1) + vp(l,1)*wght(l,3)/nvolb
!         ct(i,j,kp,2) = ct(i,j,kp,2) + vp(l,2)*wght(l,3)/nvolb
!         ct(i,j,kp,3) = ct(i,j,kp,3) + vp(l,3)*wght(l,3)/nvolb

!         nvolb = np(ip,j,kp)*volb
!         ct(ip,j,kp,1) = ct(ip,j,kp,1) + vp(l,1)*wght(l,4)/nvolb
!         ct(ip,j,kp,2) = ct(ip,j,kp,2) + vp(l,2)*wght(l,4)/nvolb
!         ct(ip,j,kp,3) = ct(ip,j,kp,3) + vp(l,3)*wght(l,4)/nvolb

!         nvolb = np(i,jp,k)*volb
!         ct(i,jp,k,1) = ct(i,jp,k,1) + vp(l,1)*wght(l,5)/nvolb
!         ct(i,jp,k,2) = ct(i,jp,k,2) + vp(l,2)*wght(l,5)/nvolb
!         ct(i,jp,k,3) = ct(i,jp,k,3) + vp(l,3)*wght(l,5)/nvolb

!         nvolb = np(ip,jp,k)*volb
!         ct(ip,jp,k,1) = ct(ip,jp,k,1) + vp(l,1)*wght(l,6)/nvolb
!         ct(ip,jp,k,2) = ct(ip,jp,k,2) + vp(l,2)*wght(l,6)/nvolb
!         ct(ip,jp,k,3) = ct(ip,jp,k,3) + vp(l,3)*wght(l,6)/nvolb

!         nvolb = np(i,jp,kp)*volb
!         ct(i,jp,kp,1) = ct(i,jp,kp,1) + vp(l,1)*wght(l,7)/nvolb
!         ct(i,jp,kp,2) = ct(i,jp,kp,2) + vp(l,2)*wght(l,7)/nvolb
!         ct(i,jp,kp,3) = ct(i,jp,kp,3) + vp(l,3)*wght(l,7)/nvolb

!         nvolb = np(ip,jp,kp)*volb
!         ct(ip,jp,kp,1) = ct(ip,jp,kp,1) + vp(l,1)*wght(l,8)/nvolb
!         ct(ip,jp,kp,2) = ct(ip,jp,kp,2) + vp(l,2)*wght(l,8)/nvolb
!         ct(ip,jp,kp,3) = ct(ip,jp,kp,3) + vp(l,3)*wght(l,8)/nvolb

            enddo
            
            !Used for periodic boundary conditions
            call add_boundary_vector(ct)
            
!            ct(nx-1,:,:,:) = ct(nx-1,:,:,:)+ct(1,:,:,:)
!            ct(:,ny-1,:,:) = ct(:,ny-1,:,:)+ct(:,1,:,:)
!            ct(:,:,nz-1,:) = ct(:,:,nz-1,:)+ct(:,:,1,:)
            
           
            call boundary_vector(ct)
            !call periodic(ct)
            
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(ct(:,:,:,:),recvbuf,count,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
            
            ct(:,:,:,:) = reshape(recvbuf,(/nx,ny,nz,3/))
            up=ct
             
!            do i=1,nx-1                 !interpolate back to contravarient positions  (now in separate routine)
!                  do j=1,ny-1
!                        do k=1,nz-1
!                              up(i,j,k,1) = 0.5*(ct(i,j,k,1)+ct(i+1,j,k,1))
!                              up(i,j,k,2) = 0.5*(ct(i,j,k,2)+ct(i,j+1,k,2))
!                              up(i,j,k,3) = 0.5*(ct(i,j,k,3)+ct(i,j,k+1,3))
!                        enddo
!                  enddo
!            enddo
            call boundary_vector(up)      
!            call periodic(up)
            
      end subroutine update_up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine up_fld_mv()
            use dimensions
            use boundary
            use var_arrays, only: up
            real:: up_temp(nx,ny,nz,3)
            integer:: i,j,k
            up_temp=up
            do i=1,nx-1                 !interpolate back to contravarient positions.
                  do j=1,ny-1
                        do k=1,nz-1
                              up(i,j,k,1) = 0.5*(up_temp(i,j,k,1)+up_temp(i+1,j,k,1))
                              up(i,j,k,2) = 0.5*(up_temp(i,j,k,2)+up_temp(i,j+1,k,2))
                              up(i,j,k,3) = 0.5*(up_temp(i,j,k,3)+up_temp(i,j,k+1,3))
                        enddo
                  enddo
            enddo
            
            call boundary_vector(up)
            
      end subroutine up_fld_mv
            
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine update_up_2()
            use dimensions
            use MPI
            use mult_proc, only: procnum
            use boundary
            use var_arrays, only: vp,up,Ni_tot,ijkp,wght
            implicit none
            real:: recvbuf(nx*ny*nz*3),ct(nx,ny,nz,3),cnt(nx,ny,nz)
            integer:: i,j,k,l,m,ip,jp,kp,count,ierr
            
            count = nx*ny*nz*3
            
            do i=1,nx
                  do j=1,ny
                        do k=1,nz
                              up(i,j,k,m) = 0.0
                              ct(i,j,k,m) = 0.0
                        enddo
                  enddo
            enddo
            
            cnt(:,:,:) = 0.0
            
            do l=1,Ni_tot
                  i=ijkp(l,1)
                  j=ijkp(l,2)
                  k=ijkp(l,3)
                  
                  ip=i+1
                  jp=j+1
                  kp=k+2
                  
!         nvolb = 1.0
!         if (np(i,j,k) .gt. 0.0) then
!         nvolb = np(i,j,k)*volb
         cnt(i,j,k) = cnt(i,j,k) + wght(l,1)
         ct(i,j,k,1) = ct(i,j,k,1) + vp(l,1)*wght(l,1) 
         ct(i,j,k,2) = ct(i,j,k,2) + vp(l,2)*wght(l,1) 
         ct(i,j,k,3) = ct(i,j,k,3) + vp(l,3)*wght(l,1) 
         
!         endif

!         if (np(ip,j,k) .gt. 0.0) then
!         nvolb = np(ip,j,k)*volb
         cnt(ip,j,k) = cnt(ip,j,k) + wght(l,2)
         ct(ip,j,k,1) = ct(ip,j,k,1) + vp(l,1)*wght(l,2) 
         ct(ip,j,k,2) = ct(ip,j,k,2) + vp(l,2)*wght(l,2) 
         ct(ip,j,k,3) = ct(ip,j,k,3) + vp(l,3)*wght(l,2) 
!         endif

!         if (np(i,j,kp) .gt. 0.0) then
!         nvolb = np(i,j,kp)*volb
         cnt(i,j,kp) = cnt(i,j,kp) + wght(l,3)
         ct(i,j,kp,1) = ct(i,j,kp,1) + vp(l,1)*wght(l,3) 
         ct(i,j,kp,2) = ct(i,j,kp,2) + vp(l,2)*wght(l,3) 
         ct(i,j,kp,3) = ct(i,j,kp,3) + vp(l,3)*wght(l,3) 
!         endif

!         if (np(ip,j,kp) .gt. 0.0) then
!         nvolb = np(ip,j,kp)*volb
         cnt(ip,j,kp) = cnt(ip,j,kp) + wght(l,4) 
         ct(ip,j,kp,1) = ct(ip,j,kp,1) + vp(l,1)*wght(l,4) 
         ct(ip,j,kp,2) = ct(ip,j,kp,2) + vp(l,2)*wght(l,4) 
         ct(ip,j,kp,3) = ct(ip,j,kp,3) + vp(l,3)*wght(l,4) 
!         endif

!         if (np(i,jp,k) .gt. 0.0) then
!         nvolb = np(i,jp,k)*volb
         cnt(i,jp,k) = cnt(i,jp,k) + wght(l,5)
         ct(i,jp,k,1) = ct(i,jp,k,1) + vp(l,1)*wght(l,5) 
         ct(i,jp,k,2) = ct(i,jp,k,2) + vp(l,2)*wght(l,5) 
         ct(i,jp,k,3) = ct(i,jp,k,3) + vp(l,3)*wght(l,5) 
!         endif

!         if (np(ip,jp,k) .gt. 0.0) then
!         nvolb = np(ip,jp,k)*volb
         cnt(ip,jp,k) = cnt(ip,jp,k) + wght(l,6)
         ct(ip,jp,k,1) = ct(ip,jp,k,1) + vp(l,1)*wght(l,6) 
         ct(ip,jp,k,2) = ct(ip,jp,k,2) + vp(l,2)*wght(l,6) 
         ct(ip,jp,k,3) = ct(ip,jp,k,3) + vp(l,3)*wght(l,6) 
!         endif

!         if (np(i,jp,kp) .gt. 0.0) then
!         nvolb = np(i,jp,kp)*volb
         cnt(i,jp,kp) = cnt(i,jp,kp) + wght(l,7)
         ct(i,jp,kp,1) = ct(i,jp,kp,1) + vp(l,1)*wght(l,7) 
         ct(i,jp,kp,2) = ct(i,jp,kp,2) + vp(l,2)*wght(l,7) 
         ct(i,jp,kp,3) = ct(i,jp,kp,3) + vp(l,3)*wght(l,7) 
!         endif

!         if (np(ip,jp,kp) .gt. 0.0) then
!         nvolb = np(ip,jp,kp)*volb
         cnt(ip,jp,kp) = cnt(ip,jp,kp) + wght(l,8)
         ct(ip,jp,kp,1) = ct(ip,jp,kp,1) + vp(l,1)*wght(l,8) 
         ct(ip,jp,kp,2) = ct(ip,jp,kp,2) + vp(l,2)*wght(l,8) 
         ct(ip,jp,kp,3) = ct(ip,jp,kp,3) + vp(l,3)*wght(l,8) 
!         endif

            enddo
            
            !Used for periodic boundary conditions
            call add_boundary_vector(ct)
            call add_boundary_scalar(cnt)
            
!            ct(nx-1,:,:,:) = ct(nx-1,:,:,:)+ct(1,:,:,:)
!            cnt(nx-1,:,:) = cnt(nx-1,:,:)+cnt(1,:,:)

!            ct(:,ny-1,:,:) = ct(:,ny-1,:,:)+ct(:,1,:,:)
!            cnt(:,ny-1,:) = cnt(:,ny-1,:)+cnt(:,1,:)

!            ct(:,:,nz-1,:) = ct(:,:,nz-1,:)+ct(:,:,1,:)
!            cnt(:,:,nz-1) = cnt(:,:,nz-1)+cnt(:,:,1)
            
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            
            where(cnt(:,:,:) .gt. 0.0)
                  ct(:,:,:,1) = ct(:,:,:,1)/cnt(:,:,:)/procnum
                  ct(:,:,:,2) = ct(:,:,:,2)/cnt(:,:,:)/procnum
                  ct(:,:,:,3) = ct(:,:,:,3)/cnt(:,:,:)/procnum
            endwhere
            
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            
            call MPI_ALLREDUCE(ct(:,:,:,:),recvbuf,count,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
            
            ct(:,:,:,:) = reshape(recvbuf,(/nx,ny,nz,3/))
            
            call boundary_vector(ct)
!            call periodic(ct)
            
            do i=1,nx-1
                  do j=1,ny-1
                        do k=1,nz-1
                              up(i,j,k,1) = 0.5*(ct(i,j,k,1)+ct(i+1,j,k,1))
                              up(i,j,k,2) = 0.5*(ct(i,j,k,2)+ct(i,j+1,k,2))
                              up(i,j,k,3) = 0.5*(ct(i,j,k,3)+ct(i,j,k+1,3))
                        enddo
                  enddo
            enddo
            
            call boundary_vector(up)
!            call periodic(up)
            
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            
      end subroutine update_up_2
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine separate_np(flg)
! Weight density to eight nearest grid points
            use dimensions
            use MPI
            use boundary
            use grid, only: qx,qy,qz
            use var_arrays, only: Ni_tot,ijkp,beta,beta_p,wght,np
            implicit none
            integer, intent(in):: flg(Ni_max)
            real:: recvbuf(nx*ny*nz), volb
            integer:: i,j,k,l,ip,jp,kp,count,ierr
            
            count = nx*ny*nz
            
            do i=1,nx
                  do j=1,ny
                        do k=1,nz
                              np(i,j,k) = 0.0
                        enddo
                  enddo
            enddo
            
            do l=1,Ni_tot
                  i=ijkp(l,1)
                  j=ijkp(l,2)
                  k=ijkp(l,3)
                  
                  ip=i+1
                  jp=j+1
                  kp=k+1
                  
                  volb = (qx(ip)-qx(i))*(qy(jp)-qy(j))*(qz(kp)-qz(k))*beta*beta_p(l)
                  
                  np(i,j,k) = np(i,j,k) + flg(l)*wght(l,1)/volb
                  np(ip,j,k) = np(ip,j,k) + flg(l)*wght(l,2)/volb
                  np(i,j,kp) = np(i,j,kp) + flg(l)*wght(l,3)/volb
                  np(ip,j,kp) = np(ip,j,kp) + flg(l)*wght(l,4)/volb
                  np(i,jp,k) = np(i,jp,k) + flg(l)*wght(l,5)/volb
                  np(ip,jp,k) = np(ip,jp,k) + flg(l)*wght(l,6)/volb
                  np(i,jp,kp) = np(i,jp,kp) + flg(l)*wght(l,7)/volb
                  np(ip,jp,kp) = np(ip,jp,kp) + flg(l)*wght(l,8)/volb
                  
            enddo
            
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(np(:,:,:),recvbuf,count,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
            
            np(:,:,:) = reshape(recvbuf,(/nx,ny,nz/))
            
            !use for periodic boundary conditions
            call add_boundary_scalar(np)
            
!            np(nx-1,:,:) = np(nx-1,:,:)+np(1,:,:)
!            np(:,ny-1,:) = np(:,ny-1,:)+np(:,1,:)
!            np(:,:,nz-1) = np(:,:,nz-1)+np(:,:,1)
            
            call boundary_scalar(np)
!            call periodic_scalar(np)
            
      end subroutine separate_np
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine separate_up(flg)
            use dimensions
            use MPI
            use boundary
            use grid, only: qx,qy,qz
            use var_arrays, only: up,Ni_tot,ijkp,beta,beta_p,wght, vp, np
            implicit none
            integer, intent(in):: flg(Ni_max)
            real:: recvbuf(nx*ny*nz*3),volb,nvolb,ct(nx,ny,nz,3)
            integer:: i,j,k,l,m,ip,jp,kp,count,ierr
            
            count = nx*ny*nz*3
            do i=1,nx
                  do j=1,ny
                        do k=1,nz
                              do m=1,3
                                    up(i,j,k,m) = 0.0
                                    ct(i,j,k,m) = 0.0
                              enddo
                        enddo
                  enddo
            enddo
            
            do l=1,Ni_tot
                  i=ijkp(l,1)
                  j=ijkp(l,2)
                  k=ijkp(l,3)
                  
                  ip=i+1
                  jp=j+1
                  kp=k+1
                  
                  volb = (qx(ip)-qx(i))*(qy(jp)-qy(j))*(qz(kp)-qz(k))*beta*beta_p(l)
                  
         if (np(i,j,k) .gt. 0.0) then
         nvolb = np(i,j,k)*volb
         ct(i,j,k,1) = ct(i,j,k,1) + flg(l)*vp(l,1)*wght(l,1)/nvolb
         ct(i,j,k,2) = ct(i,j,k,2) + flg(l)*vp(l,2)*wght(l,1)/nvolb
         ct(i,j,k,3) = ct(i,j,k,3) + flg(l)*vp(l,3)*wght(l,1)/nvolb
         endif

         if (np(ip,j,k) .gt. 0.0) then
         nvolb = np(ip,j,k)*volb
         ct(ip,j,k,1) = ct(ip,j,k,1) + flg(l)*vp(l,1)*wght(l,2)/nvolb
         ct(ip,j,k,2) = ct(ip,j,k,2) + flg(l)*vp(l,2)*wght(l,2)/nvolb
         ct(ip,j,k,3) = ct(ip,j,k,3) + flg(l)*vp(l,3)*wght(l,2)/nvolb
         endif

         if (np(i,j,kp) .gt. 0.0) then
         nvolb = np(i,j,kp)*volb
         ct(i,j,kp,1) = ct(i,j,kp,1) + flg(l)*vp(l,1)*wght(l,3)/nvolb
         ct(i,j,kp,2) = ct(i,j,kp,2) + flg(l)*vp(l,2)*wght(l,3)/nvolb
         ct(i,j,kp,3) = ct(i,j,kp,3) + flg(l)*vp(l,3)*wght(l,3)/nvolb
         endif

         if (np(ip,j,kp) .gt. 0.0) then
         nvolb = np(ip,j,kp)*volb
         ct(ip,j,kp,1) = ct(ip,j,kp,1) + flg(l)*vp(l,1)*wght(l,4)/nvolb
         ct(ip,j,kp,2) = ct(ip,j,kp,2) + flg(l)*vp(l,2)*wght(l,4)/nvolb
         ct(ip,j,kp,3) = ct(ip,j,kp,3) + flg(l)*vp(l,3)*wght(l,4)/nvolb
         endif

         if (np(i,jp,k) .gt. 0.0) then
         nvolb = np(i,jp,k)*volb
         ct(i,jp,k,1) = ct(i,jp,k,1) + flg(l)*vp(l,1)*wght(l,5)/nvolb
         ct(i,jp,k,2) = ct(i,jp,k,2) + flg(l)*vp(l,2)*wght(l,5)/nvolb
         ct(i,jp,k,3) = ct(i,jp,k,3) + flg(l)*vp(l,3)*wght(l,5)/nvolb
         endif

         if (np(ip,jp,k) .gt. 0.0) then
         nvolb = np(ip,jp,k)*volb
         ct(ip,jp,k,1) = ct(ip,jp,k,1) + flg(l)*vp(l,1)*wght(l,6)/nvolb
         ct(ip,jp,k,2) = ct(ip,jp,k,2) + flg(l)*vp(l,2)*wght(l,6)/nvolb
         ct(ip,jp,k,3) = ct(ip,jp,k,3) + flg(l)*vp(l,3)*wght(l,6)/nvolb
         endif

         if (np(i,jp,kp) .gt. 0.0) then
         nvolb = np(i,jp,kp)*volb
         ct(i,jp,kp,1) = ct(i,jp,kp,1) + flg(l)*vp(l,1)*wght(l,7)/nvolb
         ct(i,jp,kp,2) = ct(i,jp,kp,2) + flg(l)*vp(l,2)*wght(l,7)/nvolb
         ct(i,jp,kp,3) = ct(i,jp,kp,3) + flg(l)*vp(l,3)*wght(l,7)/nvolb
         endif

         if (np(ip,jp,kp) .gt. 0.0) then
         nvolb = np(ip,jp,kp)*volb
         ct(ip,jp,kp,1) = ct(ip,jp,kp,1) + flg(l)*vp(l,1)*wght(l,8)/nvolb
         ct(ip,jp,kp,2) = ct(ip,jp,kp,2) + flg(l)*vp(l,2)*wght(l,8)/nvolb
         ct(ip,jp,kp,3) = ct(ip,jp,kp,3) + flg(l)*vp(l,3)*wght(l,8)/nvolb
         endif
         
            enddo
            
            !Used for periodic boundary conditions
            call add_boundary_vector(ct)
            
!            ct(nx-1,:,:,:) = ct(nx-1,:,:,:)+ct(1,:,:,:)
!            ct(:,ny-1,:,:) = ct(:,ny-1,:,:)+ct(:,1,:,:)
!            ct(:,:,nz-1,:) = ct(:,:,nz-1,:)+ct(:,:,1,:)
            
            call boundary_vector(ct)
!            call periodic(ct)
            
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(ct(:,:,:,:),recvbuf,count,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
            
            ct(:,:,:,:) = reshape(recvbuf,(/nx,ny,nz,3/))
            
            do i=1,nx-1         !interpolate back to contravarient positions
                  do j=1,ny-1
                        do k=1,nz-1
                              up(i,j,k,1) = 0.5*(ct(i,j,k,1)+ct(i+1,j,k,1))
                              up(i,j,k,2) = 0.5*(ct(i,j,k,2)+ct(i,j+1,k,2))
                              up(i,j,k,3) = 0.5*(ct(i,j,k,3)+ct(i,j,k+1,3))
                        enddo
                  enddo
            enddo
            
            call boundary_scalar(up)
!            call periodic(up)
            
      end subroutine separate_up
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_temperature()
            use dimensions
            use MPI
            use boundary
            use inputs, only: mion
            use grid, only: qx,qy,qz
            use var_arrays, only: vp,np,temp_p,Ni_tot,ijkp,beta,beta_p,mrat,wght
            implicit none
            real:: recvbuf(nx*ny*nz*3),up2(nx,ny,nz,3),up_ave(nx,ny,nz,3),ct(nx,ny,nz,3),volb,nvolb,mvp(Ni_max,3)
            integer:: i,j,k,l,m,ip,jp,kp,count,ierr
            
            count = nx*ny*nz*3
            
            up2(:,:,:,:) = 0.0
            up_ave(:,:,:,:) = 0.0
            ct(:,:,:,:) = 0.0
            
            do m=1,3
                  mvp(1:Ni_tot,m) = vp(1:Ni_tot,m)/sqrt(mrat(:))
            enddo
            
            do l=1,Ni_tot
                  i=ijkp(l,1)
                  j=ijkp(l,2)
                  k=ijkp(l,3)
                  
                  ip=i+1
                  jp=j+1
                  kp=k+1
                  
                  volb = (qx(ip)-qx(i))*(qy(jp)-qy(j))*(qz(kp)-qz(k))*beta*beta_p(l)
                  
         if (np(i,j,k) .gt. 0.0) then
         nvolb = np(i,j,k)*volb
         ct(i,j,k,1) = ct(i,j,k,1) + mvp(l,1)**2*wght(l,1)/nvolb  
         ct(i,j,k,2) = ct(i,j,k,2) + mvp(l,2)**2*wght(l,1)/nvolb
         ct(i,j,k,3) = ct(i,j,k,3) + mvp(l,3)**2*wght(l,1)/nvolb
         endif

         if (np(ip,j,k) .gt. 0.0) then
         nvolb = np(ip,j,k)*volb
         ct(ip,j,k,1) = ct(ip,j,k,1) + mvp(l,1)**2*wght(l,2)/nvolb
         ct(ip,j,k,2) = ct(ip,j,k,2) + mvp(l,2)**2*wght(l,2)/nvolb
         ct(ip,j,k,3) = ct(ip,j,k,3) + mvp(l,3)**2*wght(l,2)/nvolb
         endif

         if (np(i,j,kp) .gt. 0.0) then
         nvolb = np(i,j,kp)*volb
         ct(i,j,kp,1) = ct(i,j,kp,1) + mvp(l,1)**2*wght(l,3)/nvolb
         ct(i,j,kp,2) = ct(i,j,kp,2) + mvp(l,2)**2*wght(l,3)/nvolb
         ct(i,j,kp,3) = ct(i,j,kp,3) + mvp(l,3)**2*wght(l,3)/nvolb
         endif

         if (np(ip,j,kp) .gt. 0.0) then
         nvolb = np(ip,j,kp)*volb
         ct(ip,j,kp,1) = ct(ip,j,kp,1) + mvp(l,1)**2*wght(l,4)/nvolb
         ct(ip,j,kp,2) = ct(ip,j,kp,2) + mvp(l,2)**2*wght(l,4)/nvolb
         ct(ip,j,kp,3) = ct(ip,j,kp,3) + mvp(l,3)**2*wght(l,4)/nvolb
         endif

         if (np(i,jp,k) .gt. 0.0) then
         nvolb = np(i,jp,k)*volb
         ct(i,jp,k,1) = ct(i,jp,k,1) + mvp(l,1)**2*wght(l,5)/nvolb
         ct(i,jp,k,2) = ct(i,jp,k,2) + mvp(l,2)**2*wght(l,5)/nvolb
         ct(i,jp,k,3) = ct(i,jp,k,3) + mvp(l,3)**2*wght(l,5)/nvolb
         endif

         if (np(ip,jp,k) .gt. 0.0) then
         nvolb = np(ip,jp,k)*volb
         ct(ip,jp,k,1) = ct(ip,jp,k,1) + mvp(l,1)**2*wght(l,6)/nvolb
         ct(ip,jp,k,2) = ct(ip,jp,k,2) + mvp(l,2)**2*wght(l,6)/nvolb
         ct(ip,jp,k,3) = ct(ip,jp,k,3) + mvp(l,3)**2*wght(l,6)/nvolb
         endif

         if (np(i,jp,kp) .gt. 0.0) then
         nvolb = np(i,jp,kp)*volb
         ct(i,jp,kp,1) = ct(i,jp,kp,1) + mvp(l,1)**2*wght(l,7)/nvolb
         ct(i,jp,kp,2) = ct(i,jp,kp,2) + mvp(l,2)**2*wght(l,7)/nvolb
         ct(i,jp,kp,3) = ct(i,jp,kp,3) + mvp(l,3)**2*wght(l,7)/nvolb
         endif

         if (np(ip,jp,kp) .gt. 0.0) then
         nvolb = np(ip,jp,kp)*volb
         ct(ip,jp,kp,1) = ct(ip,jp,kp,1) + mvp(l,1)**2*wght(l,8)/nvolb
         ct(ip,jp,kp,2) = ct(ip,jp,kp,2) + mvp(l,2)**2*wght(l,8)/nvolb
         ct(ip,jp,kp,3) = ct(ip,jp,kp,3) + mvp(l,3)**2*wght(l,8)/nvolb
         endif
         
            enddo
            
            !Used for periodic boundary conditions
            call add_boundary_vector(ct)
!            ct(nx-1,:,:,:) = ct(nx-1,:,:,:)+ct(1,:,:,:)
!            ct(:,ny-1,:,:) = ct(:,ny-1,:,:)+ct(:,1,:,:)
!            ct(:,:,nz-1,:) = ct(:,:,nz-1,:)+ct(:,:,1,:)
            
            call boundary_vector(ct)
            
!            call periodic(ct)
            
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(ct(:,:,:,:),recvbuf,count,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
            
            ct(:,:,:,:) = reshape(recvbuf,(/nx,ny,nz,3/))
            
            do i=1,nx-1
                  do j=1,ny-1
                        do k=1,nz-1
                              up2(i,j,k,1) = 0.5*(ct(i,j,k,1)+ct(i+1,j,k,1))
                              up2(i,j,k,2) = 0.5*(ct(i,j,k,2)+ct(i,j+1,k,2))
                              up2(i,j,k,3) = 0.5*(ct(i,j,k,3)+ct(i,j,k+1,3))
                        enddo
                  enddo
            enddo
            
            call boundary_vector(up2)
!            call periodic(up2)
            
            ct(:,:,:,:) = 0.0
            
            do l=1,Ni_tot
                  i=ijkp(l,1)
                  j=ijkp(l,2)
                  k=ijkp(l,3)
                  
                  ip=i+1
                  jp=j+1
                  kp=k+1
                  
                  volb = (qx(ip)-qx(i))*(qy(jp)-qy(j))*(qz(kp)-qz(k))*beta*beta_p(l)
                  
         if (np(i,j,k) .gt. 0.0) then
         nvolb = np(i,j,k)*volb
         ct(i,j,k,1) = ct(i,j,k,1) + mvp(l,1)*wght(l,1)/nvolb  
         ct(i,j,k,2) = ct(i,j,k,2) + mvp(l,2)*wght(l,1)/nvolb
         ct(i,j,k,3) = ct(i,j,k,3) + mvp(l,3)*wght(l,1)/nvolb
         endif

         if (np(ip,j,k) .gt. 0.0) then
         nvolb = np(ip,j,k)*volb
         ct(ip,j,k,1) = ct(ip,j,k,1) + mvp(l,1)*wght(l,2)/nvolb
         ct(ip,j,k,2) = ct(ip,j,k,2) + mvp(l,2)*wght(l,2)/nvolb
         ct(ip,j,k,3) = ct(ip,j,k,3) + mvp(l,3)*wght(l,2)/nvolb
         endif

         if (np(i,j,kp) .gt. 0.0) then
         nvolb = np(i,j,kp)*volb
         ct(i,j,kp,1) = ct(i,j,kp,1) + mvp(l,1)*wght(l,3)/nvolb
         ct(i,j,kp,2) = ct(i,j,kp,2) + mvp(l,2)*wght(l,3)/nvolb
         ct(i,j,kp,3) = ct(i,j,kp,3) + mvp(l,3)*wght(l,3)/nvolb
         endif

         if (np(ip,j,kp) .gt. 0.0) then
         nvolb = np(ip,j,kp)*volb
         ct(ip,j,kp,1) = ct(ip,j,kp,1) + mvp(l,1)*wght(l,4)/nvolb
         ct(ip,j,kp,2) = ct(ip,j,kp,2) + mvp(l,2)*wght(l,4)/nvolb
         ct(ip,j,kp,3) = ct(ip,j,kp,3) + mvp(l,3)*wght(l,4)/nvolb
         endif

         if (np(i,jp,k) .gt. 0.0) then
         nvolb = np(i,jp,k)*volb
         ct(i,jp,k,1) = ct(i,jp,k,1) + mvp(l,1)*wght(l,5)/nvolb
         ct(i,jp,k,2) = ct(i,jp,k,2) + mvp(l,2)*wght(l,5)/nvolb
         ct(i,jp,k,3) = ct(i,jp,k,3) + mvp(l,3)*wght(l,5)/nvolb
         endif

         if (np(ip,jp,k) .gt. 0.0) then
         nvolb = np(ip,jp,k)*volb
         ct(ip,jp,k,1) = ct(ip,jp,k,1) + mvp(l,1)*wght(l,6)/nvolb
         ct(ip,jp,k,2) = ct(ip,jp,k,2) + mvp(l,2)*wght(l,6)/nvolb
         ct(ip,jp,k,3) = ct(ip,jp,k,3) + mvp(l,3)*wght(l,6)/nvolb
         endif

         if (np(i,jp,kp) .gt. 0.0) then
         nvolb = np(i,jp,kp)*volb
         ct(i,jp,kp,1) = ct(i,jp,kp,1) + mvp(l,1)*wght(l,7)/nvolb
         ct(i,jp,kp,2) = ct(i,jp,kp,2) + mvp(l,2)*wght(l,7)/nvolb
         ct(i,jp,kp,3) = ct(i,jp,kp,3) + mvp(l,3)*wght(l,7)/nvolb
         endif

         if (np(ip,jp,kp) .gt. 0.0) then
         nvolb = np(ip,jp,kp)*volb
         ct(ip,jp,kp,1) = ct(ip,jp,kp,1) + mvp(l,1)*wght(l,8)/nvolb
         ct(ip,jp,kp,2) = ct(ip,jp,kp,2) + mvp(l,2)*wght(l,8)/nvolb
         ct(ip,jp,kp,3) = ct(ip,jp,kp,3) + mvp(l,3)*wght(l,8)/nvolb
         endif
         
            enddo
            !Used for periodic boundary conditions
            call add_boundary_vector(ct)
!            ct(nx-1,:,:,:) = ct(nx-1,:,:,:)+ct(1,:,:,:)
!            ct(:,ny-1,:,:) = ct(:,ny-1,:,:)+ct(:,1,:,:)
!            ct(:,:,nz-1,:) = ct(:,:,nz-1,:)+ct(:,:,1,:)
            
            call boundary_vector(ct)
!            call periodic(ct)
            
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(ct(:,:,:,:),recvbuf,count,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
            
            ct(:,:,:,:) = reshape(recvbuf,(/nx,ny,nz,3/))
            
            do i=1,nx-1
                  do j=1,ny-1
                        do k=1,nz-1
                              up_ave(i,j,k,1) = 0.5*(ct(i,j,k,1)+ct(i+1,j,k,1))
                              up_ave(i,j,k,2) = 0.5*(ct(i,j,k,2)+ct(i,j+1,k,2))
                              up_ave(i,j,k,3) = 0.5*(ct(i,j,k,3)+ct(i,j,k+1,3))
                        enddo
                  enddo
            enddo
            
            do i=1,nx
                  do j=1,ny
                        do k=1,nz
!                              temp_p(i,j,k) = (1./3.)*1e6*mion*(sqrt((up2(i,j,k,1) &
!                                    - up_ave(i,j,k,1)**2)**2 + &
!                                    (up2(i,j,k,2) - up_ave(i,j,k,2)**2)**2 + & 
!                                    (up2(i,j,k,3) - up_ave(i,j,k,3)**2)**2))  
                              temp_p(i,j,k) = (1./3.)*1e6*mion*( &
                                    up2(i,j,k,1) - up_ave(i,j,k,1)**2 + &
                                    up2(i,j,k,2) - up_ave(i,j,k,2)**2 + & 
                                    up2(i,j,k,3) - up_ave(i,j,k,3)**2)
                        enddo
                  enddo
            enddo

            call boundary_scalar(temp_p)
            
      end subroutine get_temperature
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine check_index()
            use dimensions
            use var_arrays, only: Ni_tot,ijkp,m_arr
            implicit none
            integer:: l
            
            do l=1,Ni_tot
                  if (ijkp(l,1) .gt. nx) then
                        write(*,*) 'i OB...', ijkp(l,:),m_arr(l)
                  endif
                  
                  if (ijkp(l,1) .lt. 1) then
                        write(*,*) 'i OB...', ijkp(l,:),m_arr(l)
                  endif
                  
                  if (ijkp(l,2) .gt. ny) then
                        write(*,*) 'j OB...',ijkp(l,:),m_arr(l)
                  endif

                  if (ijkp(l,2) .lt. 1) then
                        write(*,*) 'j OB...',ijkp(l,:),m_arr(l)
                  endif


                  if (ijkp(l,3) .gt. nz) then
                        write(*,*) 'k OB...',ijkp(l,:),m_arr(l)
                  endif

                  if (ijkp(l,3) .lt. 1) then
                        write(*,*) 'k OB...',ijkp(l,:),m_arr(l)
                  endif
                  
            enddo
            
      end subroutine check_index
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_pindex(i,j,k,l)
            use dimensions
            use inputs, only: dx,dy
            use grid, only: qx,qy,qz
            use var_arrays, only: ijkp,xp
            implicit none
            integer, intent(in):: l
            integer, intent(out):: i,j,k
            integer:: hi,mid
!            i=1               
!            do while (xp(l,1) .gt. qx(i))
!                  i=i+1
!            enddo
!            i=i-1
            i=1
            hi = nx
            do
                  mid = (i+hi)/2
                  if (xp(l,1) .lt. qx(mid)) then
                        hi=mid
                  else
                        i=mid
                  endif
                  if (i+1 .ge. hi) exit
            enddo

            ijkp(l,1)=i
            j = floor(xp(l,2)/dy)
            ijkp(l,2) = j
!            k=1
!            do while (xp(l,3) .gt. qz(k))
!                  k=k+1
!            enddo
!            k=k-1
            k=1
            hi = nz
            do
                  mid = (k+hi)/2
                  if (xp(l,3) .lt. qz(mid)) then
                        hi=mid
                  else
                        k=mid
                  endif
                  if (k+1 .ge. hi) exit
            enddo
            ijkp(l,3) = k


            
      end subroutine get_pindex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
                  
end module gutsp
