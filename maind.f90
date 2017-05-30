program hybrid
      
      use mpi
      use dimensions
      use Var_Arrays
      use mult_proc
      use inputs
      use misc
      use initial
      use part_init
      use gutsf
      use gutsp
      use grid_interp
      use chem_rates
      use grid
      use iso_fortran_env
      
      implicit none
      
      real:: Evp,Euf,EB1,EB1x,EB1y,EB1z,EE,EeP,input_EeP
      real:: input_chex, input_bill,dtsub
      real:: pup(3),puf(3),peb(3)
      character(2):: filenum!(16) !max 16 processors
      character(1):: mstart
      integer:: ierr,t1,t2,cnt_rt,m,mstart_n,ndiag,seed
      real(8):: time
      logical:: restart = .false.
      integer(4):: Ni_tot_sw!,Ni_tot_sys
      integer:: i,j,k,n,ntf !looping indicies
      real (real64) :: dp
      integer, parameter :: dp_kind = kind(dp)
      
!      filenum = (/'1 ','2 ','3 ','4 ','5 ','6 ','7 ','8 ','9 ', &
!            '10','11','12','13','14','15','16'/)
                 
      call readInputs()
      call initparameters()
      
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, procnum, ierr)
      
      
      call system_clock(t1, cnt_rt)
      
!Initialize all variables

      write(filenum, '(I2)') my_rank
      
      !if (int(ppc/procnum) .gt. 0) Ni_tot=(nx-2)*(ny-2)*(nz-2)*int(ppc/procnum) !3D
      !if (int(ppc/procnum) .eq. 0) Ni_tot=(nx-2)*(ny-2)*(nz-2)
      Ni_tot = int((real(nx)-2.0)*(real(ny)-2.0)*(real(nz)-2.0)*real(ppc)/real(procnum),4)
      Ni_tot_0 = Ni_tot
      Ni_tot_sw = Ni_tot
      Ni_tot_sys = Ni_tot*procnum
      
      if (my_rank .eq. 0) then
            call check_inputs()
            write(*,*) 'Total particles, PPP, #pc', Ni_tot_sys,Ni_tot,procnum
            write(*,*) 'Partilces per cell... ', Ni_tot_sys/((nz-2)*(ny-2)*(nx-2))
            write(*,*) ' '
      endif
      
      mstart_n = 0 !number of times restarted
      write(mstart, '(I1)') mstart_n
      
      ndiag = 0
      prev_Etot = 1.0
!      nuei = 0.0

      seed = t1 + my_rank*100
!      call random_initialize(seed)
      call seed_mpi(my_rank)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      
      if (.not. restart) then
!            do i=1,nx
!                  do j=1,ny
!                        do k=1,nz
                              input_E = 0.0
                              input_p = 0.0
                              input_chex = 0.0
                              input_bill = 0.0
                             
                              input_Eb = 0.0
                              input_EeP= 0.0
!                        enddo
!                  enddo
!            enddo
      endif
      

      
      call grd7()
!      call grid_gaussian()
      call grd6_setup(b0,bt,b12,b1,b1p2,nu,input_Eb)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call get_beta(Ni_tot_sys,beta)
   
      input_E = 0.0
      bndry_Eflux = 0.0
      
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
      !Initialize particles: use load Maxwellian, or sw_part_setup, etc.
!      call load_Maxwellian(vth,1,mion,1.0)
!      call load_const_ppc(vth,1,mion,1.0)
!      call load_RT_pad(vth,1,mion,1.0)
      call loadRT_ppc(vth,mion,1.0)
  
!      call load_den_grad(1,mion,1.0)
!      call count_ppc()
  
      if (my_rank .eq. 0) then
            call check_inputs()     
      endif
!      call load_ring_beam(57.0,int(Ni_tot*0.1),mion,1.0)
      
!      write(*,*) 'Particles per cell... (Ni_tot_sys/(nz-2)', Ni_tot_sys/(nz-2)
      
      call f_update_tlev(b1,b12,b1p2,bt,b0)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!  Check for restart flag
      if (my_rank .eq. 0) write(*,*) 'restart status...', restart
      if ((restart) .and. (mstart_n .gt. 0)) then
!            if (my_rank .eq. 0) then
            write(*,*) 'opening restar vars....'
            open(210,file=trim(out_dir)//'restart.vars',status='unknown',form='unformatted')
            write(*,*) 'reading restart vars......'
!            read(210) b1,b12,b1p2,bt,btmf,btc,np,np3,vp,vp1,vplus,vminus, &
!                  up,xp,aj,nu,Ep,E,temp_p,mnp,beta,beta_p,Evp,Euf, &
!                  EB1,EB1x,EB1y,EB1z,EE,EeP,input_E,Ni_tot, &
!                  ijkp,input_p,input_EeP,input_Eb,prev_Etot,bndry_Eflux,grav, &
!                  input_chex,input_bill,mrat,m_arr
                  
            read(210) b1,b12,b1p2,bt,btmf,btc,np,np3, &
                        up,aj,nu,E,temp_p,mnp,beta,Evp,Euf, &
                        EB1,EB1x,EB1y,EB1z,EE,EeP, &
                        input_EeP,input_Eb,prev_Etot,bndry_Eflux,grav, &
                        input_chex,input_bill
            write(*,*) 'restarting hybrid ....'
            
!            endif
            
!            if (my_rank .gt. 0) then
                  open(211,file=trim(out_dir)//'restart.part'//trim(filenum),status='unknown',form='unformatted')
                  read(211) vp,vp1,vplus,vminus,xp,Ep,input_E,Ni_tot,ijkp,input_p,mrat,m_arr,beta_p
!            endif
            close(210)
            close(211)
      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write parameter file
            if (my_rank .eq. 0) then

                  open(109, file=trim(out_dir)//'para.dat', &
                        status='unknown',form='unformatted')
         
                  write(109) nx,ny,nz,dx,dy,delz
                  write(109) nt,dtsub_init,ntsub,dt,nout,mion
                  write(109) out_dir
                  write(109) vtop,vbottom
                  write(109) Ni_max
                  write(109) mion,m_pu,m_heavy
                  write(109) np_top,np_bottom
                  write(109) b0_top,b0_bottom
                  write(109) vth_top,vth_bottom
                  write(109) alpha,beta
                  close(109)
                  
! Write fft parameter file
!                  open(401, file = trim(out_dir)//'fft_11400.dat',status='unknown',form='unformatted')
!                  write(401) dt,nt,omega_p
                  
!                  open(402, file = trim(out_dir)//'fft_14000.dat',status='unknown',form='unformatted')
!                  write(402) dt,nt,omega_p
                  
!                  open(403, file = trim(out_dir)//'fft_17000.dat',status='unknown',form='unformatted')
!                  write(403) dt,nt,omega_p

            endif
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Inititalize diagnositc output files

      if ((my_rank .eq. 0) .and. restart) then
            open(110,file=trim(out_dir)//'c.np_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(115,file=trim(out_dir)//'c.np_b_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(120,file=trim(out_dir)//'c.mixed_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(130,file=trim(out_dir)//'c.b1_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(140,file=trim(out_dir)//'c.aj_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(150,file=trim(out_dir)//'c.E_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(160,file=trim(out_dir)//'c_.energy_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            
            !diagnostics chex,bill,satnp
            
            open(180,file=trim(out_dir)//'c.up_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(190,file=trim(out_dir)//'c.momentum_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(192,file=trim(out_dir)//'c.p_conserve_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(300,file=trim(out_dir)//'c.temp_p_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(305,file=trim(out_dir)//'c.xp_0_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(310,file=trim(out_dir)//'c.vp_0_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(315,file=trim(out_dir)//'c.mrat_0_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(317,file=trim(out_dir)//'c.beta_p_0_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(320,file=trim(out_dir)//'c.np_wake_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(330,file=trim(out_dir)//'c.up_t_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
!            open(340,file=trim(out_dir)//'c.up_b_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(342,file=trim(out_dir)//'c.test_part_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(350,file=trim(out_dir)//'c.mnp_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            
       endif
       if ((my_rank .eq. 0) .and. (.not. restart)) then
            open(110,file=trim(out_dir)//'c.np.dat',status='unknown',form='unformatted')
            open(115,file=trim(out_dir)//'c.np_b.dat',status='unknown',form='unformatted')
            open(120,file=trim(out_dir)//'c.mixed.dat',status='unknown',form='unformatted')
            open(130,file=trim(out_dir)//'c.b1.dat',status='unknown',form='unformatted')
            open(140,file=trim(out_dir)//'c.aj.dat',status='unknown',form='unformatted')
            open(150,file=trim(out_dir)//'c.E.dat',status='unknown',form='unformatted')
            open(160,file=trim(out_dir)//'c.energy.dat',status='unknown',form='unformatted')
!            open(170,file=trim(out_dir)//'c.vdist_init.dat',status='unknown',form='unformatted')
!            open(175,file=trim(out_dir)//'c.vdist_add.dat',status='unknown',form='unformatted')
!            open(176,file=trim(out_dir)//'c.vpp_init.dat',status='unknown',form='unformatted')
!            open(177,file=trim(out_dir)//'c.vpp_add.dat',status='unknown',form='unformatted')
            
            !diagnostics chex,bill,satnp
            
            open(180,file=trim(out_dir)//'c.up.dat',status='unknown',form='unformatted')
            open(190,file=trim(out_dir)//'c.momentum.dat',status='unknown',form='unformatted')
            open(192,file=trim(out_dir)//'c.p_conserve.dat',status='unknown',form='unformatted')
            open(300,file=trim(out_dir)//'c.temp_p.dat',status='unknown',form='unformatted')
            open(305,file=trim(out_dir)//'c.xp_0.dat',status='unknown',form='unformatted')
            open(310,file=trim(out_dir)//'c.vp_0.dat',status='unknown',form='unformatted')
            open(315,file=trim(out_dir)//'c.mrat_0.dat',status='unknown',form='unformatted')
            open(317,file=trim(out_dir)//'c.beta_p_0.dat',status='unknown',form='unformatted')
            open(320,file=trim(out_dir)//'c.np_wake.dat',status='unknown',form='unformatted')
            open(330,file=trim(out_dir)//'c.up_t.dat',status='unknown',form='unformatted')
            open(340,file=trim(out_dir)//'c.up_b.dat',status='unknown',form='unformatted')
            open(342,file=trim(out_dir)//'c.test_part.dat',status='unknown',form='unformatted')
            open(350,file=trim(out_dir)//'c.mnp.dat',status='unknown',form='unformatted')
       
            open(501,file=trim(out_dir)//'c.energy_p.dat',status='unknown',form='unformatted')
       endif
       
       if (my_rank .gt. 0) then
!            open(305,file=trim(out_dir)//'c.xp_'//trim(filenum)//'_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
!            open(310,file=trim(out_dir)//'c.vp_'//trim(filenum)//'_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
!            open(315,file=trim(out_dir)//'c.mrat_'//trim(filenum)//'_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
!            open(317,file=trim(out_dir)//'c.beta_p__'//trim(filenum)//'_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
       endif
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       MAIN LOOP
      do m = 1, nt !mstart_n+1, nt
            if (my_rank .eq. 0) then
                  write(*,*) 'time...', m, m*dt,my_rank
            endif
          !  if (m .lt. 600) then
                  !Call ionizing subroutine  (adds ions to the domain)
          !        call Mass_load_Io(m)
          !  endif
            call get_interp_weights()
            call update_np()                  !np at n+1/2
            call update_up(vp)            !up at n+1/2
          
            !energy diagnostics
            call get_bndry_Eflux(b1,E,bndry_Eflux)

            call Energy_diag(Evp,Euf,EB1,EB1x,EB1y,EB1z,EE,EeP)
            call get_gradP()
 
            call curlB(bt,np,aj)
            call edge_to_center(bt,btc)
            call extrapol_up()
            call get_Ep()

            
            call get_vplus_vminus()
            call improve_up()

            call get_Ep()
         
            call get_vplus_vminus()
            call get_vp_final()
      
            call move_ion_half() !1/2 step ion move to n+1/2

            call get_interp_weights()

            call update_np()                  !np at n+1/2

            call update_up(vp)            !up at n+1/2
            
            call get_gradP()
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Subcycling loop

            dtsub = dtsub_init
            ntf=ntsub
            call check_time_step(bt,np,dtsub,ntf)
            
            do n=1,ntf
                  call curlB(bt,np,aj)
                  
                  call predict_B(b12,b1p2,bt,E,aj,up,nu,dtsub)
                  
                  call correct_B(b1,b1p2,E,aj,up,np,nu,dtsub)
                  
                  call f_update_tlev(b1,b12,b1p2,bt,b0)
                  
            enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           
 
            call move_ion_half()       !final ion move to n+1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       diagnositc output
            
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            
            if (my_rank .eq. 0) then
                  write(160) m
                  write(160) input_E,input_EeP,Evp,Euf,EB1,EB1x,EB1y,EB1z,EE, &
                        EeP,input_chex,input_bill
                  write(190) m
                  write(190) pup,puf,peb,input_p
                  
                  !fft output
!                  write(401) b1(2,2,11400,1), b1(2,2,11400,2), b1(2,2,11400,3)
!                  write(402) b1(2,2,14000,1), b1(2,2,14000,2), b1(2,2,14000,3)
!                  write(403) b1(2,2,17000,1), b1(2,2,17000,2), b1(2,2,17000,3)
                  
            endif
            
            ndiag = ndiag+1
            if (ndiag .eq. nout) then
                  call get_temperature()
                  call update_rho()
!                  call get_v_dist()
                  call update_mixed(mixed,mix_cnt)
                  if (my_rank .eq. 0) then
                        write(110) m
                        write(110) np
                        write(115) m
                        write(115) np_b
                        write(120) m
                        write(120) mixed
                        write(130) m
                        write(130) bt
                        write(140) m
                        write(140) aj*alpha*np3
                        write(150) m
                        write(150) E
                        write(180) m
                        write(180) up
                        write(300) m
                        write(300) temp_p/1.6e-19  !output in eV
!                        write(305) m
!                        write(305) xp
!                        write(310) m
!                        write(310) vp
!                        write(315) m
!                        write(315) mrat
!                        write(317) m
!                        write(317) beta_p
                        write(330) m
                        write(330) up_t
!                        write(340) m
!                        write(340) up_b
                        write(350) m
                        write(350) mnp
                        
!                        write(170) m
!                        write(170) vdist_init
!                        write(175) m
!                        write(175) vdist_add
!                        write(176) m
!                        write(176) vpp_init
!                        write(177) m
!                        write(177) vpp_add
                        
!                        write(342) m
!                        write(342) vp(299985:Ni_max,:)
                        
                        ndiag = 0
                   endif
                   
                   if (my_rank .gt. 0) then
!                        write(305) m
!                        write(305) xp
!                        write(310) m
!                        write(310) vp
!                        write(315) m
!                        write(315) mrat
!                        write(317) m
!                        write(317) beta_p
                        
                        ndiag = 0
                  endif
                   
            endif
!            write(*,*) 'minimum density.....', minval(np(:,:,:))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Write restart file

            
            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
      enddo     !End Main Loop
      
!      if (my_rank .eq. 0) then
            if (restart) then
                  if (my_rank .eq. 0) then
                        write(*,*) 'Writing restart file...'
                        
                        !Write grid data
                        
                        open(220,file=trim(out_dir)//'restart.vars',status='unknown',form='unformatted')
                        write(220) b1,b12,b1p2,bt,btmf,btc,np,np3, &
                        up,aj,nu,E,temp_p,mnp,beta,Evp,Euf, &
                        EB1,EB1x,EB1y,EB1z,EE,EeP, &
                        input_EeP,input_Eb,prev_Etot,bndry_Eflux,grav, &
                        input_chex,input_bill
                  endif
                              
                              
                  !write individual processor data
                  open(212,file=trim(out_dir)//'restart.part'//trim(filenum),status='unknown',form='unformatted')
                  write(212) vp,vp1,vplus,vminus,xp,Ep,input_E,Ni_tot,ijkp,input_p,mrat,m_arr,beta_p
                  
                  close(220)
                  close(212)
            endif
!      endif
      
      close(110)
      close(115)
      close(120)
      close(130)
      close(140)
      close(150)
      close(160)
!      close(170)
      close(172)
!      close(175)
!      close(176)
!      close(177)
      close(180)
      close(190)
      close(192)
      close(210)
      close(220)
      close(221)
      close(300)
      close(305)
      close(310)
      close(315)
      close(317)
      close(320)
      close(330)
!      close(340)
      close(342)
      close(350)
!      close(401)
!      close(402)
      close(403)
      close(501)
      
      call system_clock(t2,cnt_rt)
      time=(real(t2,dp_kind) - real(t1,dp_kind))/real(cnt_rt,dp_kind)
      if (my_rank .eq. 0) then
            write(*,*)
            write(*,*) 'Elapsed time .....', time, ' sec'
            write(*,*)
      endif
      
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      
      call MPI_FINALIZE(ierr)

end program hybrid

      
