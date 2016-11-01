program predictor
      use MPI
      use dimensions
      implicit none
      integer:: ierr, my_rank, procnum, ppcpp, Ni_tot, Ni_tot_sys, ppc
      real:: rvar
      
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, procnum, ierr)
      
      if (my_rank .eq. 0) then
                 open(unit=100, file= 'inputs.dat', status='old')
                  
                 read(100,*) rvar
                 read(100,*) rvar
                 read(100,*) rvar
                 read(100,*) rvar
                 read(100,*) rvar
                 read(100,*) rvar
                 read(100,*) rvar
                 read(100,*) rvar
                 read(100,*) rvar
                 read(100,*) rvar
                 read(100,*) rvar
                 read(100,*) rvar
                 read(100,*) rvar
                 read(100,*) rvar
                 read(100,*) rvar
                 read(100,*) rvar
                 read(100,*) ppc
            write(*,*) 'Particle per cell input.....',ppc 
                 close(100)
            ppcpp=int(ppc/procnum)
            Ni_tot = (nx-2)*(ny-2)*(nz-2)*ppcpp
            Ni_tot_sys = Ni_tot*procnum
            write(*,*) 'Particles per processor.....', Ni_tot
            write(*,*) 'Allocated particles.........', Ni_max
            write(*,*) 'Particles per cell..........', Ni_tot_sys/((nx-2)*(ny-2)*(nz-2)), ppcpp*procnum
            write(*,*) 'Total particles.............', Ni_tot_sys
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