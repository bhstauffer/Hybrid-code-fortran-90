program testmpi
      use mpi
      implicit none
      integer:: ierr
      real:: ion_amu = 23
      
      call MPI_INIT(ierr)
      
      write(*,*) 23.0*1.67e-27

      call MPI_FINALIZE(ierr)

end program testmpi

      