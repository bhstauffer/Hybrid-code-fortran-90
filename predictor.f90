program predictor
      use MPI
      use dimensions
      use inputs
      use grid
      use initial
      implicit none
      integer:: ierr, my_rank, procnum, ppcpp, Ni_tot, Ni_tot_sys, Ni_tot2
      !integer*8:: Ni_tot_sys, Ni_tot2
      real:: rvar
      call readInputs()
      
      call initparameters()
      
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, procnum, ierr)
      
      if (my_rank .eq. 0) then
                 !open(unit=100, file= 'inputs.dat', status='old')
                  
                 !read(100,*) rvar
                 !read(100,*) rvar
                 !read(100,*) rvar
                 !read(100,*) rvar
                 !read(100,*) rvar
                 !read(100,*) rvar
                 !read(100,*) rvar
                 !read(100,*) rvar
                 !read(100,*) rvar
                 !read(100,*) rvar
                 !read(100,*) rvar
                 !read(100,*) rvar
                 !read(100,*) rvar
                 !read(100,*) rvar
                 !read(100,*) rvar
                 !read(100,*) rvar
                 !read(100,*) ppc
            !write(*,*) 'Particle per cell input.....',ppc 
            !     close(100)
            call grd7()
            !call grd6_setup(b0,bt,b12,b1,b1p2,nu,input_Eb)
            ppcpp=int(ppc/procnum)
            !if (ppcpp .gt. 0) Ni_tot = (nx-2)*(ny-2)*(nz-2)*ppcpp
            !if (ppcpp .eq. 0) Ni_tot = (nx-2)*(ny-2)*(nz-2)
            Ni_tot = int((real(nx)-2.0)*(real(ny)-2.0)*(real(nz)-2.0)*real(ppc)/real(procnum),4)
            Ni_tot_sys = Ni_tot*procnum
            
            Ni_tot2 = nint( Ni_tot / 2 * amp * ( (qz(nz-1) - Lo * log(cosh((qz(nz-1)-qz(nz/2))/Lo))) - &
                  (qz(1) - Lo * log(cosh((qz(1) - qz(nz/2))/Lo))) ) / (qz(nz-1)-qz(1)) )
            write(*,*) 'Particles per processor.....', Ni_tot
            write(*,*) 'Allocated particles.........', Ni_max
            write(*,*) 'Particles per cell..........', Ni_tot_sys/((nx-2)*(ny-2)*(nz-2)), ppcpp*procnum
            write(*,*) 'Total particles.............', Ni_tot_sys
            write(*,*) 'RT_pad particles per proc...', Ni_tot+Ni_tot2
            if (Ni_tot .lt. Ni_max) then
                  write(*,*) 'Sufficient memory allocated.'
            else
                  write(*,*) 'Insufficient memory allocated.'
            endif
            write(*,*) '***************************************************'
      endif
      
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_FINALIZE(ierr)


end program predictor