module part_init
      implicit none
      save
      contains

      subroutine Energy_diag(Evp,Euf,EB1,EB1x,EB1y,EB1z,EE,EeP)
            use dimensions
            use mpi
            use mult_proc, only: my_rank
            use grid, only: dx_cell,dy_cell,dz_cell
            use inputs, only: mion,q,mu0,mO,km_to_m,epsilon
            use var_arrays, only: vp,b0,b1,E,nu,up,np,Ni_tot,beta,beta_p,input_E,prev_Etot,bndry_Eflux,m_arr
            implicit none
            real, intent(out):: Euf,EB1,EB1x,EB1y,EB1z,EE,EeP,Evp
            real:: denf,m_q,recvbuf,total_E,aveEvp,norm_E,vol
            real:: S_Evp,S_input_E
            integer:: count, ierr
            integer:: i,j,k,m,l
            
            count = 1
            m_q = mion/q
            
            Euf = 0.0
            EB1 = 0.0
            EB1x = 0.0
            EB1y = 0.0
            EB1z = 0.0
            EE = 0.0
            EeP = 0.0
            
            do i=1,nx-1
!                  do j = 1,ny-1
                   j=2
                        do k=1,nz-1
                              vol= dx_cell(i)*dy_cell(j)*dz_cell(k)*km_to_m**3
                              EB1x=EB1x + (vol/(2.0*mu0))*(m_q*b1(i,j,k,1))**2
                              EB1y=EB1y + (vol/(2.0*mu0))*(m_q*b1(i,j,k,2))**2
                              EB1z=EB1z + (vol/(2.0*mu0))*(m_q*b1(i,j,k,3))**2
                                    do m=1,3
                                          denf = np(i,j,k)/(km_to_m**3)
                                          Euf = Euf + 0.5*mO*denf*vol*(up(i,j,k,m)*km_to_m)**2
                                          EB1 = EB1 + (vol/(2.0*mu0))*(m_q*(b1(i,j,k,m)-b0(i,j,k,m)))**2
                                          EE = EE + (epsilon*vol/2.0)*(m_q*E(i,j,k,m)*km_to_m)**2
                                    enddo
                        enddo
!                  enddo
            enddo
            
            Evp = 0.0
            do l=1, Ni_tot
                  do m=1,3
                        Evp = Evp + 0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 / (beta*beta_p(l))
                  enddo
            enddo
            
            
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            
            call MPI_ALLREDUCE(Evp,recvbuf,count,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
            S_Evp = recvbuf
            
            call MPI_ALLREDUCE(input_E,recvbuf,count,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
            S_input_E = recvbuf
            
            total_E = S_Evp + EE + EB1
            aveEvp = S_Evp/S_input_E
            
            
            if (my_rank .eq. 0) then
                  write(*,*) 'Normalized energy.................',total_E/S_input_E
                  write(*,*) 'Normalized energy (bndry).........',total_E/(S_input_E+bndry_Eflux)
!                  write(*,*) '          particles...............',S_Evp/S_input_E       
!                  write(*,*) '          field...................',(EE+EB1)/S_input_E
!                  write(501) m
!                  write(501) S_Evp/S_input_E
!                  write(501) (EE+EB1)/S_input_E
            endif
            
            norm_E = total_E/S_input_E
            prev_Etot = norm_E

      end subroutine Energy_diag
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine load_Maxwellian(pbeta,Ni_tot_1,mass,mratio)
            use dimensions
            use boundary
            use inputs, only: PI, vsw, dx, dy, km_to_m, beta_particle, kboltz, mion, amp, grad, nf_init,b0_init,mu0
            use grid, only: qx,qy,qz,dz_grid
            use gutsp
            use var_arrays, only: np,vp,vp1,xp,input_p,up,Ni_tot,input_E,ijkp,m_arr,mrat,beta,beta_p,wght,grav,temp_p
            implicit none
            integer(4), intent(in):: Ni_tot_1 
            real, intent(in):: mratio, mass, pbeta
                                  
            integer:: disp
            real:: vx, vy, vz, Temp, Tempcalc, vth
            integer:: l,m,i,j,k
            
            disp = 0 !Displacement of gradient
!            amp = 100.0
!            grad = 100.0 ! density gradient (larger = more gradual
            

!            va = b0_init/sqrt(mu0*mion*nf_init/1e9)/1e3
      !      vth = sqrt(pbeta*va**2)
            
            do l = Ni_tot_1,Ni_tot
                  xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
                  xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
                  xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz-1)-qz(1))
                  m_arr(l) = mass
                  mrat(l) = mratio

!                  beta_p(l) = 1.0/(beta_particle+beta_particle*amp*exp(-((xp(l,3)-qz(nz/2-disp))/ &
!                        (grad*dz_grid(nz/2-disp)))**2))
                  beta_p(l) = beta_particle

                  call get_pindex(i,j,k,l)
!                  vth2=sqrt(vth*vth*beta_p(l)) !thermal speed dependent on np to set up pressure balance for density gradient

!                  vth2=va*sqrt(pl_beta(ijkp(l,1),ijkp(l,2),ijkp(l,3)))
                 

                  
                  vx = 40.0*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf()) !remember to add in vsw to get the flow velocity
                  vy = 40.0*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
                  vz = 40.0*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
                  
!                  ii = ijkp(l,1)
!                  kk = ijkp(l,3)
                  
!                  vp(l,1) = -0.0*(exp(-(xp(l,3)-qz(nz/2))**2/(10.*delz)**2)
!               x        *exp(-(xp(l,1)-qx(nx/2))**2/(10.*dx)**2))+vx
                  vp(l,1) = vx+57.0*exp(-(xp(l,3)-qz(nz/2))**2/(5*dz_grid(nz/2))**2) !Gaussian velocity perturbation (20)
                  vp(l,2) = vy 
                  vp(l,3) = vz 
                  
                  do m=1,3
                        vp1(l,m) = vp(l,m)
                        input_E = input_E + 0.5*m_arr(l)*(vp(l,m)*km_to_m)**2/(beta * beta_p(l))
                        input_p(m) = input_p(m) + m_arr(l) * vp(l,m) / (beta * beta_p(l))
                  enddo
                  
            enddo
            
            call get_interp_weights()
            call update_np()
            call update_up(vp)

            
            ! Add a centrifugal gravity term to keep the plasma confined to the torus.  Use T * dn/dz = nmg.  
            ! Depends on the density gradient.  Currently set as a gaussian.
            
!            Temp = vth**2/(3*kboltz)*mion*1.48*10-23!8.61738e-5
!            write(*,*) 'vth.................', vth
!            write(*,*) 'boltzman............', kboltz
!            write(*,*) 'temperature(analytic)..', Temp
!            call get_temperature()
!            Tempcalc = sum(temp_p(2,2,1:(nz-1)))/1e6/(nz-1) !in kg km^2/s^2
!            write(*,*) 'temperature (2,2,100)..', temp_p(2,2,2:10)/1.6e-19
!            stop
            
            do i=1,nx
            do j=1,ny
            do k=1,nz
                  ! Gravity is based on the analytical expression for the density profile (look at beta_p)
                  ! np = const/(beta*beta_p), and grav = const * (dn/dx) / n
                  

                        
                       
                        
!                  grav(i,j,k) = -2.0*Tempcalc/(mion*(grad*dz_grid(nz/2-disp))**2 &
!                        *(1.0+amp*exp(-((qz(k)-qz(nz/2-disp))/(grad*dz_grid(nz/2-disp)))**2))) &
!                        *amp*(qz(k)-qz(nz/2-disp))*exp(-((qz(k)-qz(nz/2-disp))/(grad*dz_grid(nz/2-disp)))**2)   
                 
       
                 grav(i,j,k) = 0.0
                  
                  
            enddo
            enddo
            enddo
           ! write(*,*) 'gravity...', grav(2,2,nz/2+50), grav(2,2,nz/2-50)
           ! stop

           call count_ppc() 
      end subroutine load_Maxwellian
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine load_RT_pad(pbeta,Ni_tot_1,mass,mratio)
            use dimensions
            use mult_proc
            use MPI
            use boundary
            use inputs, only: PI, vsw, dx, dy, km_to_m, beta_particle, kboltz, mion, amp, grad, delz, Lo, nf_init, va
            use grid, only: qx,qy,qz,dz_grid
            use gutsp
            use var_arrays, only: vp,vp1,xp,input_p,Ni_tot,input_E,m_arr,mrat,beta,beta_p,grav,mix_ind
            implicit none
            integer(4), intent(in):: Ni_tot_1 
            real, intent(in):: mratio, mass, pbeta
            real:: npb(nz),npt(nz), np_tot(nz)  !need to establish analytical density profile                     
           

            integer:: disp
            real:: vth, vx, vy, vz, grav0, rho(nz),vbal(nz), np_top,np_bottom
            integer:: l,m,i,j,k,flg
            real:: rnd
            integer:: Ni_tot_p_t,Ni_tot_2

            vth = 50.0!sqrt(pbeta*va**2)
            grav0= .4*vth**2/Lo
            do i=1,nx
               do j=1,ny
                  grav(i,j,:)=-grav0*(tanh((qz(:)-qz(nz/2) - 60*delz)/(10*delz))-tanh((qz(:)-qz(nz/2)+60*delz)/(10*delz)))
               enddo
            enddo
            np_top = nf_init
            np_bottom = nf_init*amp
            npb(:) = np_bottom*((1.0-tanh((qz(:)-qz(nz/2))/(Lo)))/2.0) 
            npt(:) = np_top
            np_tot = npb+npt
            
!            rho = (amp-1)/2*tanh((qz(:)-qz(nz/2))/Lo)+(amp+1)/2

            vbal(1) = vth
            do k = 2,nz 
               vbal(k) = sqrt((vbal(1))**2 - (vth**2/np_top)*(npb(k)-npb(1)) + &
                    2*sum(grav(1,1,1:k)*np_tot(1:k))*delz/np_top)
!               if (my_rank .eq. 0) &
!                    write(*,*) 'vth_bal...',k,vbal(k),grav(1,1,k),grav0
            enddo
!            stop

!            do k=1,nz
!                  vbal(k) = sqrt(3*sum(rho(k:nz)*grav(nx/2,2,k:nz)*dz_grid(k:nz))/rho(k))
!!                  if (my_rank .eq. 0) write(*,*) 'vbal...',k,vbal(k)
!            enddo

            
            do l = 1,Ni_tot
               xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
               xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))

               flg = 0
               do while (flg .eq. 0)
                  xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz-1)-qz(1))
                  rnd = ((1.0-tanh((xp(l,3)-qz(nz/2))/(Lo)))/2.0)  !for bottom
                  if (pad_ranf() .le. rnd) flg = 1
                  m_arr(l) = mion
                  mrat(l) = 1.0
               enddo

               beta_p(l) = 1.0
               call get_pindex(i,j,k,l)
               
               
               if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
               if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
               
               
               vp(l,1) = (vth)*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
               vp(l,2) = (vth)*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
               vp(l,3) = (vth)*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())+&
                    0.0*exp(-(xp(l,3)-qz(nz/2))**2/(5*dz_grid(nz/2))**2)*cos(2*PI*(qx(nx-1)-xp(l,1))/(qx(nx-1)-qx(1)))
               
               do m=1,3
                  vp1(l,m) = vp(l,m)
                  input_E = input_E + 0.5*m_arr(l)*(vp(l,m)*km_to_m)**2/(beta * beta_p(l))
                  input_p(m) = input_p(m) + m_arr(l) * vp(l,m) / (beta * beta_p(l))
               enddo
            enddo


            Ni_tot_2 = Ni_tot + 1
            
            ! Ni_tot = beta / procnum * integral(N_top dv)
            ! Ni_tot2 = beta/ procnum * integral( N_bot dz ) * [x]*[y]
            ! integral of (1-tanh(z-z0/Lo))/2 = (z-L0*ln(cosh(z-z0/Lo)))/2
            ! Ni_tot 2 = Ni_tot * 1/2 * (amp-1) * (z-L0*ln(cosh(z-z0/Lo))) | (qz(1) to qz(nz-1)) / (qz(nz-1)-qz(1))
            Ni_tot_p_t = nint( Ni_tot / 2 * amp * ( (qz(nz-1) - Lo * log(cosh((qz(nz-1)-qz(nz/2))/Lo))) - &
                  (qz(1) - Lo * log(cosh((qz(1) - qz(nz/2))/Lo))) ) / (qz(nz-1)-qz(1)) )
                  !write(*,*) 'one.....', Ni_tot_p_t
            !Ni_tot_p_t = amp*nint(Ni_tot*(np_top/np_bottom))
                  !write(*,*) 'two.....', Ni_tot_p_t
            Ni_tot = Ni_tot + Ni_tot_p_t 
            !write(*,*) 'Ni_tot for each processor..........', Ni_tot
            !stop
            do l = Ni_tot_2,Ni_tot
               xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
               xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
               xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz-1)-qz(1))
               m_arr(l) = mass
               mrat(l) = mratio
               ! density = nbot ( tanh(x-Lx/2 / Lo) + 1 ) amp/2
               
               beta_p(l) = 1.0
!               beta_p(l) = 1.0/(beta_particle*(tanh((xp(l,3)-qz(nz/2))/Lo)*(amp-1)/2+(amp+1)/2))
!               !                  beta_p(l) = beta_particle
               call get_pindex(i,j,k,l)
               
               
               if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
               if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
               
               
               vp(l,1) = (vbal(k))*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
               vp(l,2) = (vbal(k))*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
               vp(l,3) = (vbal(k))*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())+&
                    0.0*exp(-(xp(l,3)-qz(nz/2))**2/(5*dz_grid(nz/2))**2)*cos(2*PI*(qx(nx-1)-xp(l,1))/(qx(nx-1)-qx(1)))
               
               do m=1,3
                  vp1(l,m) = vp(l,m)
                  input_E = input_E + 0.5*m_arr(l)*(vp(l,m)*km_to_m)**2/(beta * beta_p(l))
                  input_p(m) = input_p(m) + m_arr(l) * vp(l,m) / (beta * beta_p(l))
               enddo
            enddo

            !convert gravity into non-normalized units (km/s)
!            grav = -grav
            call boundary_scalar(grav)
            
!            call MPI_BARRIER(MPI_COMM_WORLD,m)
!            call update_mixed()
!            if (my_rank .eq. 0) then
!            endif
!            stop
            call get_interp_weights()
            call update_np()    
            call update_up(vp)
                  
      end subroutine load_RT_pad
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine loadRT_ppc(pbeta,mass,mratio)
            use dimensions
            use mult_proc, only: procnum, my_rank
            use MPI
            use boundary
            use inputs, only: PI,km_to_m, beta_particle, amp, dx,dy,delz,grad, delz, nf_init, va,Lo,ratio 
            use grid, only: qx,qy,qz,dz_grid
            use gutsp
            use var_arrays, only: vp,vp1,xp,input_p,Ni_tot,input_E,m_arr,mrat,beta,beta_p,grav,mix_ind,ijkp
            implicit none
            integer(4):: Ni_tot_1 
            real, intent(in):: mratio, mass, pbeta
            real:: npb(nz),npt(nz), np_tot(nz)  !need to establish analytical density profile                     
            real:: vth, vx, vy, vz, grav0,vbal(nz)
            integer:: l,m,i,j,k,flg, count1(nx,ny,nz)
            real:: rnd, np_top,np_bottom, vol
            !ratio = 0.75
            !Lo = 4.0*delz
            vth = 50.0!sqrt(pbeta*va**2)
            grav0= .4*vth**2/Lo  !0.4 originally
            beta = 1.0
            vol = (qx(nx-1)-qx(1))*(qy(ny-1)-qy(1))*(qz(nz-1)-qz(1))!dx*dy*delz
            do i=1,nx
               do j=1,ny
                  grav(i,j,:)=-grav0*(tanh((qz(:)-qz(nz/2) - 60*delz)/(10*delz))-tanh((qz(:)-qz(nz/2)+60*delz)/(10*delz)))
               enddo
            enddo
            np_top = nf_init
            np_bottom = nf_init*(amp-1.0)
            npb(:) = np_bottom*((1.0-tanh((qz(:)-qz(nz/2))/(Lo)))/2.0) 
            npt(:) = np_top
            np_tot = npb+npt

            vbal(1) = vth
            do k = 2,nz 
               vbal(k) = sqrt((vbal(1))**2 - (vth**2/np_top)*(npb(k)-npb(1)) + &
                    2*sum(grav(1,1,1:k)*np_tot(1:k))*delz/np_top)
            enddo
            ! Load population 1, where vth is constant, but n varies
            do l = 1, int(real(Ni_tot)*ratio,4)
                  xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
                  xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
                  !do
                        xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz-1)-qz(1))
                  !      if (pad_ranf() .le. ((1.0-tanh((xp(l,3)-qz(nz/2))/(Lo)))/2.0) ) exit
                  !enddo
                  m_arr(l) = mass
                  mrat(l) = mratio
                  beta_p(l) = real(Ni_tot*ratio*procnum)/vol/(np_bottom*((1.0-tanh((xp(l,3)-qz(nz/2))/Lo))/2.0) ) !   beta =  Ni_tot/volume  /  (density)
                  !if (my_rank .eq. 0) write(*,*) xp(l,3), xp(l,3)-qz(nz/2), (1.0-tanh((xp(l,3)-qz(nz/2))/Lo))/2.0
                  call get_pindex(i,j,k,l)
                  if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
                  if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
                  vp(l,1) = (vth)*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
                  vp(l,2) = (vth)*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
                  vp(l,3) = (vth)*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
            enddo
            !call MPI_BARRIER(MPI_COMM_WORLD,l)
            !stop
            ! load population 2, where n is constant, but vth varies
            do l = int(real(Ni_tot)*ratio,4)+1, Ni_tot
                  xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
                  xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
                  !do
                        xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz-1)-qz(1))
                        !rnd = ((1.0-tanh((xp(l,3)-qz(nz/2))/(Lo)))/2.0)
                  !if (pad_ranf() .le. ((1.0+tanh((xp(l,3)-qz(nz/2))/(Lo)))/2.0)) exit
                  !enddo
                  m_arr(l) = mass
                  mrat(l) = mratio
                  beta_p(l) = real(Ni_tot*(1.0-ratio)*procnum)/vol/np_top
                  call get_pindex(i,j,k,l)
                  if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
                  if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
                  vp(l,1) = (vbal(k))*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
                  vp(l,2) = (vbal(k))*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
                  vp(l,3) = (vbal(k))*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
            enddo
            count1 = 0
!            do l = 1, Ni_tot/2
!                  i = ijkp(l,1)
!                  j = ijkp(l,2)
!                  k = ijkp(l,3)
!                  count1(i,j,k) = count1(i,j,k) + 1
!            enddo
            
!            do l = 1,Ni_tot/2
!                  beta_p(l) = real(count1(ijkp(l,1),ijkp(l,2),ijkp(l,3)))/vol/(np_bottom*((1.0+tanh(xp(l,3)-qz(nz/2))/Lo)/2.0) )
!            enddo
!            count1 = 0
!            do l = 1, Ni_tot
!                  i = ijkp(l,1)
!                  j = ijkp(l,2)
!                  k = ijkp(l,3)
!                  count1(i,j,k) = count1(i,j,k) + 1
!            enddo
            !do l = 1,Ni_tot
            !      beta_p(l) = 1.0!real(count1(ijkp(l,1),ijkp(l,2),ijkp(l,3)))/vol/(np_top) 
            !enddo
            
            do l = 1, Ni_tot
                  do m=1,3
                    vp1(l,m) = vp(l,m)
                    input_E = input_E + 0.5*m_arr(l)*(vp(l,m)*km_to_m)**2/(beta * beta_p(l))
                    input_p(m) = input_p(m) + m_arr(l) * vp(l,m) / (beta * beta_p(l))
                  enddo
            enddo
            
            call boundary_scalar(grav)
            
            call get_interp_weights()
            call update_np()    
            call update_up(vp)
                  
      end subroutine loadRT_ppc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      subroutine load_ring_beam(vring,dNi,mass,mratio)
            use dimensions
            use boundary
            use inputs, only: PI, vsw, dx, dy, km_to_m, beta_pu, ion_amu,m_pu,beta_particle, amp, grad
            use grid, only: qx,qy,qz,dz_grid
            use gutsp
            use var_arrays, only: np,vp,vp1,xp,input_p,up,Ni_tot,input_E,ijkp,m_arr,mrat,beta,beta_p,wght
            implicit none
            integer(4), intent(in):: dNi 
            real, intent(in):: vring, mass,mratio
                                  
            integer:: disp, flg, l1
            real:: vth2, vx, vy, vz, rand1, theta2
            integer:: i,j,k,l,m
            
            disp = 0 !Displacement of gradient
!            amp = 100.0
!            grad = 400.0 ! density gradient (larger = more gradual
            
!            v1=1.0
            l1=Ni_tot+1
            
            do l = l1, l1+dni-1
                 ! xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
                  xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
                  flg=0
                        do while (flg .eq. 0)
                            !  xp(l,3) = qz(nz/2-20) + (1.0-pad_ranf())*(qz(nz/2+20)-qz(nz/2-20))
                              xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz-1)-qz(1))
                              rand1=pad_ranf()
                              if (exp(-(xp(l,3)-qz(nz/2))**2/(10*dz_grid(nz/2)**2)) .gt. rand1) then
                                    flg = 1
                              endif
                        enddo
                        
                  flg=0
                        do while (flg .eq. 0)
                              xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
                              rand1=pad_ranf()
                              if (exp(-(xp(l,1)-qx(nx/2))**2/(10*dx**2)) .gt. rand1) then
                                    flg = 1
                              endif
                        enddo
                        
                        beta_p(l) = beta_particle/10.0
                        m_arr(l) = mass
                        mrat(l) = mratio
                  
                  
                  call get_pindex(i,j,k,l)
                  
                  
!                  ii = ijkp(l,1)
!                  kk = ijkp(l,3)
                  
!                  vp(l,1) = -0.0*(exp(-(xp(l,3)-qz(nz/2))**2/(10.*delz)**2)
!               x        *exp(-(xp(l,1)-qx(nx/2))**2/(10.*dx)**2))+vx

            !   Ring beam velocity initializtion
!                  theta2 = pad_ranf()*2*PI
!                  vp(l,1) = vring*cos(theta2)
!                  vp(l,2) = vring*sin(theta2)
!                  vp(l,3) = 0.0
                  
            !   Maxwellian thermal distribution 
            
                  vth2=100.0;
                  
                  vx = vth2*sqrt(-log(pad_ranf()))*cos(2*PI*pad_ranf()) !remember to add in vsw to get the flow velocity
                  vy = vth2*sqrt(-log(pad_ranf()))*cos(2*PI*pad_ranf())
                  vz = vth2*sqrt(-log(pad_ranf()))*cos(2*PI*pad_ranf())
                  
                  
                  do m=1,3
                        vp1(l,m) = vp(l,m)
                        input_E = input_E + 0.5*m_arr(l)*(vp(l,m)*km_to_m)**2/(beta * beta_p(l))
                        input_p(m) = input_p(m) + m_arr(l) * vp(l,m) / (beta * beta_p(l))
                  enddo
                  
            enddo
            Ni_tot=Ni_tot+dNi
            
            call get_interp_weights()
            call update_np()
            call update_up(vp)
            
      end subroutine load_ring_beam
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine load_const_ppc(pbeta,begin,mass,mratio)
            use mult_proc, only: my_rank, procnum
            use dimensions
            use boundary
            use inputs, only: PI, vsw, dx, dy, km_to_m, beta_particle, kboltz, mion, amp, grad, nf_init,b0_init,mu0,ppc
            use grid, only: qx,qy,qz,dz_grid,dy_grid,dx_grid
            use gutsp
            use var_arrays, only: np,vp,vp1,xp,input_p,up,Ni_tot,Ni_tot_sys,Ni_init,input_E,ijkp,m_arr,mrat,beta,beta_p,wght,grav,temp_p
            use misc
            implicit none
            integer(4), intent(in):: begin 
            real, intent(in):: mratio, mass, pbeta
                                  
            integer:: disp,ppcpp,extra_ppc
            integer(4):: Ni_tot_1
            real:: vx, vy, vz, vth, Temp, Tempcalc,vol,new
            integer:: l,m,i,j,k
            
            disp = 0 !Displacement of gradient
            
            ppcpp=int(ppc/procnum)
        !    vth=sqrt(pbeta*va**2)
            
            Ni_tot_1 = begin
            do i=1,nx-2
            do j=1,ny-2
            do k=1,nz-2
            vol = dx_grid(i)*dy_grid(j)*dz_grid(k)  !km^3
            new = 10.0*exp(-((qz(k)-qz(nz/2-disp))/(grad*dz_grid(nz/2-disp)))**2)
            if (new .lt. 1.0) then
                  if (new .gt. pad_ranf()) then
                        new = 1.0
                  else 
                        new = 0.0
                  endif
            endif
            extra_ppc = nint(new)
            
            do l = Ni_tot_1, Ni_tot_1+ppcpp+extra_ppc-1
                  xp(l,1) = qx(i)+(1.0-pad_ranf())*(qx(i+1)-qx(i))
                  xp(l,2) = qy(j)+(1.0-pad_ranf())*(qy(j+1)-qy(j))
                  xp(l,3) = qz(k)+(1.0-pad_ranf())*(qz(k+1)-qz(k))
                  m_arr(l) = mass
                  mrat(l) = mratio

!                  beta_p(l) = 1.0/(beta_particle+beta_particle*amp*exp(-((xp(l,3)-qz(nz/2-disp))/ &
!                        (grad*dz_grid(nz/2-disp)))**2))

!                  beta_p(l) = beta_particle
                  beta_p(l) = (ppcpp+extra_ppc)*procnum/(nf_init+nf_init*(amp-1.0)*exp(-((qz(k)-qz(nz/2-disp))/(grad*dz_grid(nz/2-disp)))**2))/vol
!                   beta_p(l) = ppcpp*procnum/nf_init/vol

!                  call get_pindex(i,j,k,l)
!                  vth2=sqrt(vth*vth*beta_p(l)) !thermal speed dependent on np to set up pressure balance for density gradient

!                  vth2=va*sqrt(pl_beta(ijkp(l,1),ijkp(l,2),ijkp(l,3)))


                  
                  vx = 40.0*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf()) !remember to add in vsw to get the flow velocity
                  vy = 40.0*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
                  vz = 40.0*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
                  
!                  ii = ijkp(l,1)
!                  kk = ijkp(l,3)
                  
!                  vp(l,1) = -0.0*(exp(-(xp(l,3)-qz(nz/2))**2/(10.*delz)**2)
!               x        *exp(-(xp(l,1)-qx(nx/2))**2/(10.*dx)**2))+vx
                  vp(l,1) = vx+57.0!*exp(-(xp(l,3)-qz(nz/2))**2/(30*dz_grid(nz/2))**2) !Gaussian velocity perturbation (20)
                  vp(l,2) = vy 
                  vp(l,3) = vz 
                  
                  
                  
            enddo
            
            Ni_tot_1 = Ni_tot_1+ppcpp+extra_ppc
            Ni_tot=Ni_tot_1-1
            
            enddo
            enddo
            enddo
            
            Ni_tot_sys = Ni_tot*procnum
            beta=1.0;
            
            do l = begin,Ni_tot
                  call get_pindex(i,j,k,l)
                  do m=1,3
                        vp1(l,m) = vp(l,m)
                        input_E = input_E + 0.5*m_arr(l)*(vp(l,m)*km_to_m)**2/(beta * beta_p(l))
                        input_p(m) = input_p(m) + m_arr(l) * vp(l,m) / (beta * beta_p(l))
                  enddo
            enddo
          
            if (my_rank .eq. 0) then
                  write(*,*) 'Particles per processor.....', Ni_tot
                  write(*,*) 'Particles per cell..........', Ni_tot_sys/((nx-2)*(ny-2)*(nz-2)), ppcpp*procnum
                  write(*,*) 'Total particles.............', Ni_tot_sys
                  write(*,*) '***************************************************'
            endif
            
            call get_interp_weights()
            call update_np()
            call update_up(vp)

            
            ! Add a centrifugal gravity term to keep the plasma confined to the torus.  Use T * dn/dz = nmg.  
            ! Depends on the density gradient.  Currently set as a gaussian.
            
!            Temp = vth**2/(3*kboltz)*mion*1.48*10-23!8.61738e-5
!            write(*,*) 'vth.................', vth
!            write(*,*) 'boltzman............', kboltz
!            write(*,*) 'temperature(analytic)..', Temp
            call get_temperature()
            Tempcalc = sum(temp_p(2,2,1:(nz-1)))/1e6/(nz-1) !in kg km^2/s^2
!            write(*,*) 'temperature (2,2,100)..', temp_p(2,2,2:10)/1.6e-19
!            stop
            
            do i=1,nx
            do j=1,ny
            do k=1,nz
                  ! Gravity is based on the analytical expression for the density profile (look at beta_p)
                  ! np = const/(beta*beta_p), and grav = const * (dn/dx) / n.  Units are in km/s^2.
                  
            !      grav = (T*dn/dz)/(n*m)  where n = No*(1+(amp-1)*exp(-(z-z0)^2/dz^2))
                     
            !     grav(i,j,k) = 0.0
                 grav(i,j,k) = -2.0*Tempcalc*(amp-1.0)*(qz(k)-qz(nz/2-disp))/(grad*dz_grid(nz/2-disp))**2*exp(-((qz(k)-qz(nz/2-disp))/ &
                        (grad*dz_grid(nz/2-disp)))**2)/(mion*(1.0+(amp-1.0)*exp(-((qz(k)-qz(nz/2-disp))/(grad*dz_grid(nz/2-disp)))**2)))
                        
                  
            enddo
            enddo
            enddo
            Ni_init = Ni_tot
            call count_ppc()
      
      end subroutine load_const_ppc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine load_den_grad(begin,mass,mratio)
            use mult_proc, only: my_rank, procnum
            use dimensions
            use boundary
            use inputs, only: PI, vsw, dx, dy, km_to_m, beta_particle, kboltz, mion, amp, grad, nf_init,b0_init,mu0,ppc
            use grid, only: qx,qy,qz,dz_grid,dy_grid,dx_grid
            use gutsp
            use var_arrays, only: np,vp,vp1,xp,input_p,up,Ni_tot,Ni_tot_sys,input_E,ijkp,m_arr,mrat,beta,beta_p,wght,grav,temp_p
            use misc
            implicit none
            integer(4), intent(in):: begin 
            real, intent(in):: mratio, mass
                                  
            integer:: disp,ppcpp,extra_ppc
            integer(4):: Ni_tot_1
            real:: vx, vy, vz, vth, Temp, Tempcalc,vol,new
            integer:: l,m,i,j,k
            
            
            disp = 0 !Displacement of gradient
            beta = real(ppc)/(nf_init*dx_grid(1)*dy_grid(1)*dz_grid(1))
        !    vth=sqrt(pbeta*va**2)
            
            write(*,*) 'step 4', beta
            
            Ni_tot_1 = begin
            do i=1,nx-2
            do j=1,ny-2
            do k=1,nz-2
            vol = dx_grid(i)*dy_grid(j)*dz_grid(k)  !km^3
            ppcpp = int(ppc/procnum*beta*(nf_init+nf_init*(amp-1.0)*exp(-((qz(k)-qz(nz/2-disp))/(grad*dz_grid(nz/2-disp)))**2))*vol)
            write(*,*) 'step 5',ppcpp 
            new = 0.0
            if (new .lt. 1.0) then
                  if (new .gt. pad_ranf()) then
                        new = 1.0
                  else 
                        new = 0.0
                  endif
            endif
            extra_ppc = nint(new)
            
            do l = Ni_tot_1, Ni_tot_1+ppcpp+extra_ppc-1
                  xp(l,1) = qx(i)+(1.0-pad_ranf())*(qx(i+1)-qx(i))
                  xp(l,2) = qy(j)+(1.0-pad_ranf())*(qy(j+1)-qy(j))
                  xp(l,3) = qz(k)+(1.0-pad_ranf())*(qz(k+1)-qz(k))
                  m_arr(l) = mass
                  mrat(l) = mratio
                  beta_p(l) = 1.0
                  
                  vx = 40.0*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf()) !remember to add in vsw to get the flow velocity
                  vy = 40.0*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
                  vz = 40.0*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
                  
!                  ii = ijkp(l,1)
!                  kk = ijkp(l,3)
                  
!                  vp(l,1) = -0.0*(exp(-(xp(l,3)-qz(nz/2))**2/(10.*delz)**2)
!               x        *exp(-(xp(l,1)-qx(nx/2))**2/(10.*dx)**2))+vx
                  vp(l,1) = vx+57.0!*exp(-(xp(l,3)-qz(nz/2))**2/(30*dz_grid(nz/2))**2) !Gaussian velocity perturbation (20)
                  vp(l,2) = vy 
                  vp(l,3) = vz 
                  
                  
                  
            enddo
            
            Ni_tot_1 = Ni_tot_1+ppcpp+extra_ppc
            Ni_tot=Ni_tot_1-1
            
            enddo
            enddo
            enddo
            
            
            Ni_tot_sys = Ni_tot*procnum
            
            do l = begin,Ni_tot
                  call get_pindex(i,j,k,l)
                  do m=1,3
                        vp1(l,m) = vp(l,m)
                        input_E = input_E + 0.5*m_arr(l)*(vp(l,m)*km_to_m)**2/(beta * beta_p(l))
                        input_p(m) = input_p(m) + m_arr(l) * vp(l,m) / (beta * beta_p(l))
                  enddo
            enddo
            
            
            call get_interp_weights()
            call update_np()
            call update_up(vp)
            
            call get_temperature()
            Tempcalc = sum(temp_p(2,2,1:(nz-1)))/1e6/(nz-1) !in kg km^2/s^2
            do i=1,nx
                  do j=1,ny
                        do k=1,nz
                        grav(i,j,k) = -2.0*Tempcalc*(amp-1.0)*(qz(k)-qz(nz/2-disp))/(grad*dz_grid(nz/2-disp))**2*exp(-((qz(k)-qz(nz/2-disp))/ &
                              (grad*dz_grid(nz/2-disp)))**2)/(mion*(1.0+(amp-1.0)*exp(-((qz(k)-qz(nz/2-disp))/(grad*dz_grid(nz/2-disp)))**2)))
                        enddo
                  enddo
            enddo
     !       call count_ppc()
            
      end subroutine load_den_grad
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      subroutine load_aniso_Maxwellian(vth,Ni_tot_1)
            use dimensions
            use boundary
            use grid, only: qx,qy,qz,dz_grid
            use inputs, only: mion, dx, dy,delz, vsw, km_to_m, PI, beta_particle
            use gutsp
            use var_arrays, only: np,vp,vp1,xp,input_p,up,Ni_tot,input_E,ijkp,m_arr,mrat,beta,beta_p,wght
            implicit none
            real, intent(in):: vth
            integer(4), intent(in):: Ni_tot_1
            real:: vx,vy,vz
            real:: aniso_frac, vthx, vthy, vthz
            integer:: i,j,k,l,m
            
            aniso_frac = 0.06
            
            
            do l = 1, Ni_tot_1
                  xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
                  xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
                  xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz-1)-qz(1))
                  m_arr(l) = mion
                  mrat(l) = 1.0
                  beta_p(l) = beta_particle
                  
                  call get_pindex(i,j,k,l)
                  
                  vx = vsw+vth*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
                  vy = vth*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
                  vz = vth*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())

!                  ii = ijkp(l,1)
!                  kk = ijkp(l,3)
                  vp(l,1) = -0.0*(exp(-(xp(l,3)-qz(nz/2))**2/(10.*delz)**2) &
                        *exp(-(xp(l,1)-qx(nx/2))**2/(10.*dx)**2))+vx
                  vp(l,2) = vy 
                  vp(l,3) = vz 
                  
                  do m = 1,3
                        vp1(l,m) = vp(l,m)
                        input_E = input_E + 0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 / (beta * beta_p(l))
                        input_p(m) = input_p(m) + m_arr(l)*vp(l,m) / (beta * beta_p(l))
                  enddo
                  
            enddo
            
            vthx = 1200.0
            vthy = 1200.0
            vthz = 500.0
            
            do l = Ni_tot_1+1, Ni_tot_1 + aniso_frac*Ni_tot_1
                  xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
                  xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
                  xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz-1)-qz(1))
                  m_arr(l) = mion
                  mrat(l) = 1.0
                  beta_p(l) = beta_particle
                  
                  call get_pindex(i,j,k,l)
                  
                  vx = vsw+vthx*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
                  vy = vthy*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
                  vz = vthz*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())

!                  ii = ijkp(l,1)
!                  kk = ijkp(l,3)
                  vp(l,1) = -0.0*(exp(-(xp(l,3)-qz(nz/2))**2/(10.*delz)**2) &
                        *exp(-(xp(l,1)-qx(nx/2))**2/(10.*dx)**2))+vx
                  vp(l,2) = vy 
                  vp(l,3) = vz 
                  
                  do m = 1,3
                        vp1(l,m) = vp(l,m)
                        input_E = input_E + 0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 / (beta * beta_p(l))
                        input_p(m) = input_p(m) + m_arr(l)*vp(l,m) / (beta * beta_p(l))
                  enddo
                  
            enddo
            
            Ni_tot = Ni_tot_1 + aniso_frac*Ni_tot_1
            
            do l = Ni_tot + 1, Ni_tot + aniso_frac*Ni_tot
                  xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
                  xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
                  xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz-1)-qz(1))
                  m_arr(l) = 2.0*mion
                  mrat(l) = 0.5
                  beta_p(l) = beta_particle
                  
                  call get_pindex(i,j,k,l)
                  
                  vx = vsw+vth*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
                  vy = vth*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
                  vz = vth*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())

!                  ii = ijkp(l,1)
!                  kk = ijkp(l,3)
                  vp(l,1) = -0.0*(exp(-(xp(l,3)-qz(nz/2))**2/(10.*delz)**2) &
                        *exp(-(xp(l,1)-qx(nx/2))**2/(10.*dx)**2))+vx
                  vp(l,2) = vy 
                  vp(l,3) = vz 
                  
                  do m = 1,3
                        vp1(l,m) = vp(l,m)
                        input_E = input_E + 0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 / (beta * beta_p(l))
                        input_p(m) = input_p(m) + m_arr(l)*vp(l,m) / (beta * beta_p(l))
                  enddo
                  
            enddo
            
            Ni_tot = Ni_tot + aniso_frac*Ni_tot
            
            call get_interp_weights()
            call update_np()
            call update_up(vp)
      
      end subroutine load_aniso_Maxwellian
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine sw_part_setup_maxwl_mix()
            use dimensions
            use inputs, only: mion, vth_bottom
            use gutsp
            use var_arrays, only: np,vp,vp1,xp,input_p,up,np_t_flg,np_b_flg,Ni_tot,input_E,ijkp,m_arr,mrat,beta,beta_p,wght
            implicit none
            integer(4):: Ni_tot_1
            
            
            np_t_flg = 0
            np_b_flg = 0
            
!       Add cold populations first

            Ni_tot_1 = 1
            call load_Maxwellian(vth_bottom,Ni_tot_1,mion,1.0) !mass ratio
! add He++

!      Ni_tot_1 = Ni_tot + 1
!      Ni_tot = 2.0*Ni_tot_0
      
!      call load_Maxwellian(np,vp,vp1,xp,input_p,up,
!     x     vth_bottom, Ni_tot_1, 2.0*mion, 1.0/2.0, 10.0) !inverse ration, this is a mass 2 particle (q/m)
         

! add pickup distribution

!         Ni_tot_1 = Ni_tot + 1
!         Ni_tot = 3*Ni_tot_0

!         do 69 l = Ni_tot_1,Ni_tot
!
!            xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
!            xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
!
!           xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz-1)-qz(1))
!
!            m_arr(l) = mion
!            mrat(l) = 1.0
!            beta_p(l) = beta_pu

!            i=0
! 71         continue
!            i = i + 1
!            if (xp(l,1) .gt. qx(i)) go to 71 !find i on non-uniform 
!            i = i-1
!            ijkp(l,1)= i


!            ijkp(l,2) = floor(xp(l,2)/dy) 
            
!            k=0
! 70         continue
!            k = k + 1
!            if (xp(l,3) .gt. qz(k)) go to 70 !find k on non-uniform 
!            k = k-1
!            ijkp(l,3)= k

!            theta = pad_ranf()*2*PI
            
!            vp(l,1) = vsw+vsw*cos(theta) !+dvx
!            vp(l,2) = vsw*sin(theta) !+dvz 
!            vp(l,3) = 0.0

!            if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
!            if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
            
!            do  m=1,3
!               vp1(l,m) = vp(l,m)
!               input_E = input_E + 
!     x              0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /(beta*beta_p(l))
!               input_p(m) = input_p(m) + m_arr(l)*vp(l,m)/(beta*beta_p(l))
!            enddo
            

! 69      enddo
            call get_interp_weights()
            call update_np()
            call update_up(vp)
            
      end subroutine sw_part_setup_maxwl_mix
      

end module part_init