module misc
      implicit none
      contains
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine seed_mpi(my_rank)
            integer, intent(in):: my_rank
            integer::i, time(8)
            integer, dimension(2):: seed
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
            
            
            vol = ((qx(nx-1)-qx(1))*(qy(ny-1)-qy(1))*(qz(nz-1)-qz(1)))
            beta = (Ni_tot_sys/vol)/np_top
            
            write(*,*) 'beta....',beta
            
      end subroutine get_beta
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_gradP()
            use dimensions
            use grid_interp
            use var_arrays, only: np, gradP
            use grid, only: dx_grid, dy_grid, dz_grid
            use inputs, only: etemp0, mion, kboltz
            implicit none
            real:: np1,gdnp,a0,etemp,gnpf(nx,ny,nz,3)
            integer:: i,j,k
            
            etemp = etemp0*11604.505  !eV to Kelvin
            do i=2,nx-1
                  do j=2,ny-1
                        do k=2,nz-1
                              np1 =  0.5*(np(i+1,j,k)+np(i,j,k))
                              gdnp = (np(i+1,j,k)-np(i,j,k))/dx_grid(i)
                              a0 = kboltz*etemp/(mion*np1)
!                             a(i,j,k,1) = a0*gdnp
                              gnpf(i,j,k,1) = a0*gdnp

                              np1 =  0.5*(np(i,j+1,k)+np(i,j,k))
                              gdnp = (np(i,j+1,k)-np(i,j,k))/dy_grid(j)
                              a0 = kboltz*etemp/(mion*np1)
!                              a(i,j,k,2) = a0*gdnp
                              gnpf(i,j,k,2) = a0*gdnp

                              np1 =  0.5*(np(i,j,k+1)+np(i,j,k))
                              gdnp = (np(i,j,k+1)-np(i,j,k))/dz_grid(k)
                              a0 = kboltz*etemp/(mion*np1)
!                              a(i,j,k,3) = a0*gdnp
                              gnpf(i,j,k,3) = a0*gdnp
                        enddo
                  enddo
            enddo
            
            call face_to_center(gnpf,gradP)
            !gradP(:,:,:,:) = 0
            
      end subroutine get_gradP
            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      subroutine get_np3(nfp,np3)
!            use dimensions
!            use boundary
!            implicit none
!            real, intent(in):: nfp(nx,ny,nz)
!            real, intent(out):: np3(nx,ny,nz,3)
!            integer:: i,j,k
            
!            do i = 2, nx-1
!                  do j= 2, ny-1
!                        do k = 2, nz-1
!                              np3(i,j,k,1) = 0.5*(nfp(i,j,k)+nfp(i+1,j,k))
!                              np3(i,j,k,2) = 0.5*(nfp(i,j,k)+nfp(i,j+1,k))
!                              np3(i,j,k,3) = 0.5*(nfp(i,j,k)+nfp(i,j,k+1))
!                        enddo
!                  enddo
!            enddo
            
!            call boundary_vector(np3)
!            call periodic(np3)
            
!      end subroutine get_np3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            
end module misc