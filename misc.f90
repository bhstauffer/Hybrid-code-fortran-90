module misc
      implicit none
      contains

!!!!!!!!!RANDOM NUMBER GENERATOR!!!!!!!!!!!!!!

      real function pad_ranf()
            implicit none
            call random_number(pad_ranf)
      end function pad_ranf


      subroutine random_initialize(seed_input)
!**********************************************************************
!
!! RANDOM_INITIALIZE initializes the FORTRAN90 random number seed.
!
!  Discussion:
!    If you don't initialize the FORTRAN90 random number generator
!    routine RANDOM_NUMBER, which is used by calls like
!
!      call random_number ( harvest = x )
!
!    then its behavior is not specified.  That may be OK for you.  But
!    you may want to be able to force the same sequence to be generated
!    each time, or to force a different sequence.
!
!    To control the sequence of random numbers, you need to set the seed.
!    In FORTRAN90, this is done by calling the RANDOM+SEED routine.
!    You can call it with no arguments, in fact.  But if you call
!    it with no arguments:
!
!      call random_seed ( )
!
!    then its behavior (or more particularly, the behavior of RANDOM_NUMBER)
!    is still not specified.  You might hope that the system will go out
!    and pick a nice random seed for you, but there's no guarantee.
!
!
!    For example, on the DEC ALPHA, if you compile a program that calls
!    RANDOM_NUMBER, then every time you run it, you get the same sequence
!    of "random" values.  If you compile a program that calls RANDOM_SEED
!    with no arguments, and then calls RANDOM_NUMBER, you still get the
!    same sequence each time you run the program.
!
!    In order to actually try to scramble up the random number generator
!    a bit, this routine goes through the tedious process of getting the
!    size of the random number seed, making up values based on the current
!    time, and setting the random number seed.
!
!    Unfortunately, the RANDOM_SEED routine has a very elastic definition.
!    It does not use a single scalar integer SEED.  Instead, it communicates
!    with the user through an integer vector whose size is not specified.
!    You actually have to "ask" the routine to tell you the size of this
!    vector.  Then you can fill up the vector with values to be used to
!    influence the seeding of the random number routine.  The details of
!    how the seed affects the sequence are also unspecified, but the only
!    thing we are reasonably confident about is that specifying the same
!    seed should result in the same sequence, and specifying different
!    seeds should result in different sequences!
!
!    I assume this is far more than you wanted to know.  (It's certainly
!    more than I wanted to know!)
!
!    The argument SEED is an input quantity only, so it is legal
!    to type
!
!      call random_initialize ( 0 )
!
!    or
!
!      call random_initialize ( 18867 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!     Input, integer ( kind = 4 ) SEED_INPUT, the user's "suggestion" for a seed
!     However, if the input value is 0, the routine will come up with
!     its own "suggestion", based on the system clock.

            implicit none
            integer, intent(in):: seed_input
            integer:: seed, count, count_max, count_rate, seed_size
            integer, allocatable:: seed_vector(:)
            logical, parameter:: debug = .false.
            real:: t
            integer, parameter:: warm_up = 100
            integer:: i

            seed = seed_input

            !initialize the random seed routine
            call random_seed()
            !Determine the size of the random number seed vector
            call random_seed(size = seed_size)
            !Allocate a vector of the right size to be used as a random seed
            allocate (seed_vector(seed_size))
            !If user supplied a SEED value, use that

            !Otherwise, use the system clock value to make up a value that is likely to change based on when the routine is called.

            if (seed /= 0 ) then
                  if (debug) then
                        write(*,*) ' '
                        write(*,*) 'RANDOM_INTITIALIZE'
                        write(*,*) 'Initialize RANDOM_NUMBER, user SEED = ', seed
                  endif
            else
                  call system_clock(count,count_rate,count_max)

                  seed = count

                  if (debug) then
                        write(*,*) ' '
                        write(*,*) 'RANDOM_INITIALIZE'
                        write(*,*) 'Initialize RANDOM_NUMBER, arbitrary SEED = ',seed
                  endif
            endif

            !Set the seed vector.  Set all entries to seed

            seed_vector(1:seed_size) = seed

            !Free up the seed space

            deallocate(seed_vector)

            !Call the random number routine a bunch of time to "warm up"

            do i = 1, warm_up
                  call random_number(harvest = t)
            enddo

      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine seed_mpi(my_rank)
            integer, intent(in):: my_rank
            integer::i, time(8)
            integer, dimension(12):: seed
            real:: r

            call date_and_time(values=time)
            
            seed(:) = time(4) * ( 360000*time(5) + 6000*time(6) + 100*time(7) + time(8)) + my_rank*100
            call random_seed(PUT=seed)



            do i=1,100
                  call random_number(r)
            enddo

      end subroutine seed_mpi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_bndry_Eflux(b1,E,bndry_Eflux)
            use dimensions
            use inputs, only: mO,q,dt,dx,dy,km_to_m,mu0
            implicit none
            real, intent(in):: b1(nx,ny,nz,3), E(nx,ny,nz,3)
            real, intent(inout):: bndry_Eflux
            integer:: i,j,k
            real:: exb_flux, mO_q

            mO_q = mO/q
            !k=2 face
            do i = 2, nx
                  do j= 2, ny
!                        m=3
                        k=2
                        exb_flux = (mO_q)**2*(1.0/mu0)*dt*dx*dy* &
                              (E(i,j,k,1)*b1(i,j,k,2) - E(i,j,k,2)*b1(i,j,k,1))* &
                              km_to_m**3
                        bndry_Eflux = bndry_Eflux + exb_flux

                  enddo
            enddo

            !k=nx face

            do i =2,nx
                  do j=2,ny
                        k = nz-1
                        exb_flux = (mO_q)**2*(1.0/mu0)*dt*dx*dy* &
                              (E(i,j,k,1)*b1(i,j,k,2) - E(i,j,k,2)*b1(i,j,k,1))* &
                              km_to_m**3
                        bndry_Eflux = bndry_Eflux - exb_flux

                  enddo
            enddo

      end subroutine get_bndry_Eflux
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_beta(Ni_tot_sys,beta)
            use dimensions
            use grid, only: qx,qy,qz
            use inputs, only: np_top
            implicit none
            integer(4), intent(in):: Ni_tot_sys
            real, intent(out):: beta
            real:: vol

            vol = ((qx(nx-1)-qx(1))*(qy(ny-1)-qy(1))*(qz(nz)-qz(1)))
            beta = (Ni_tot_sys/vol)/np_top

            write(*,*) 'beta....',beta

      end subroutine get_beta
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_np3(nfp,np3)
            use dimensions
            use boundary
            implicit none
            real, intent(in):: nfp(nx,ny,nz)
            real, intent(out):: np3(nx,ny,nz,3)
            integer:: i,j,k

            do i = 2, nx-1
                  do j= 2, ny-1
                        do k = 2, nz-1
                              np3(i,j,k,1) = 0.5*(nfp(i,j,k)+nfp(i+1,j,k))
                              np3(i,j,k,2) = 0.5*(nfp(i,j,k)+nfp(i,j+1,k))
                              np3(i,j,k,3) = 0.5*(nfp(i,j,k)+nfp(i,j,k+1))
                        enddo
                  enddo
            enddo

            call periodic(np3)

      end subroutine get_np3

end module misc
