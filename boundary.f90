module boundary
      implicit none
      contains

      subroutine periodic(b)
            use dimensions
            implicit none
            real, intent(inout):: b(nx,ny,nz,3)
            integer:: i,j,k,m
            
!       X direction

            do j=1,ny
                  do k=1,nz
                        do m=1,3
                              b(1,j,k,m) = b(nx-1,j,k,m)
                              b(nx,j,k,m) = b(2,j,k,m)
                        enddo
                  enddo
            enddo

!       Y direction

            do i=1,nx
                  do k=1,nz
                        do m=1,3
                              b(i,1,k,m) = b(i,ny-1,k,m)
                              b(i,ny,k,m) = b(i,2,k,m)
                        enddo
                  enddo
            enddo

!       Z direction

            do i=1,nx
                  do j=1,ny
                        do m=1,3
                              b(i,j,1,m) = b(i,j,nz-1,m)
                              b(i,j,nz,m) = b(i,j,2,m)
                        enddo
                  enddo
            enddo

      end subroutine periodic
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------
      subroutine tangent_B_zero(b) !normal derivative = 0
! The normal derivative of the tangential components is used to 
! determine the tangential boundary values.  ALSO, the normal
! derivative of the normal components temporarily sets the boundary
! values of the normal components.  These values are undated later
! by the requirement that divB=0.  This helps with the values on
! corners and edges.
!---------------------------------------------------------------------
            use dimensions
            implicit none
            real, intent(inout):: b(nx,ny,nz,3)
            integer:: i,j,k
            
!       X surfaces
            do j=2,ny-1
                  do k=2, nz-1
                        b(1,j,k,1) = b(2,j,k,1)
                        b(1,j,k,2) = b(2,j,k,2)
                        b(1,j,k,3) = b(2,j,k,3)
                        
                        b(nx,j,k,2) = b(nx-1,j,k,2)
                        b(nx,j,k,3) = b(nx-1,j,k,3)
                        
                  enddo
            enddo
            
!       Y surfaces
            do i=2,nx-1
                  do k=2,nz-1
                        b(i,1,k,2) = b(i,2,k,2)
                        b(i,1,k,1) = b(i,2,k,1)
                        b(i,1,k,3) = b(i,2,k,3)
                        
                        b(i,ny,k,1) = b(i,ny-1,k,1)
                        b(i,ny,k,3) = b(i,ny-1,k,3)
                        
                  enddo
            enddo
            
!       Z surfaces
            do i=2,nx-1
                  do j=2,ny-1
                        b(i,j,1,3) = b(i,j,2,3)
                        b(i,j,1,1) = b(i,j,2,1)
                        b(i,j,1,2) = b(i,j,2,2)
                       
                        b(i,j,nz,1) = b(i,j,nz-1,1)
                        b(i,j,nz,2) = b(i,j,nz-1,2)
                        
                  enddo
            enddo
            
      end subroutine tangent_B_zero
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine copy_to_boundary(b)
            use dimensions
            implicit none
            real, intent(inout):: b(nx,ny,nz,3)
            integer:: i,j,k,m
            
!       X surfaces, periodic
            do j=1,ny
                  do k=1,nz
                        do m = 1,3
                              b(1,j,k,m) = b(nx-1,j,k,m)
                              b(nx,j,k,m) = b(2,j,k,m)
                        enddo
                  enddo
            enddo
            
!       Y surfaces, periodic
            do i=1,nx
                  do k=1,nz
                        do m = 1,3
                              b(i,1,k,m) = b(i,2,k,m)
                              b(i,ny,k,m) = b(i,ny-1,k,m)
                        enddo
                  enddo
            enddo          

!       Z surfaces, periodic
            do i=1,nx
                  do j=1,ny
                        do m = 1,3
                              b(i,j,1,m) = b(i,j,2,m)
                              b(i,j,nz,m) = b(i,j,nz-2,m)
                        enddo
                  enddo
            enddo    
            
      end subroutine copy_to_boundary
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine periodic_scalar(b)
            use dimensions
            implicit none
            real, intent(inout):: b(nx,ny,nz)
            integer:: i,j,k
            
!       X surfaces
            do j=1,ny
                  do k=1,nz
                        b(1,j,k) = b(nx-1,j,k)
                        b(nx,j,k) = b(2,j,k)
                  enddo
            enddo

!       Y surfaces
            do i=1,nx
                  do k=1,nz
                        b(i,1,k) = b(i,ny-1,k)
                        b(i,ny,k) = b(i,2,k)
                  enddo
            enddo
            
!       Z surfaces
            do i=1,nx
                  do j=1,ny
                        b(i,j,1) = b(i,j,nz-1)
                        b(i,j,nz) = b(i,j,2)
                  enddo
            enddo    
            
      end subroutine periodic_scalar
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine fix_normal_b(b)
            use dimensions
            use grid, only: dx_grid,dy_grid,dz_grid
            use inputs, only: dx,dy
            implicit none
            real, intent(inout):: b(nx,ny,nz,3)
            integer:: i,j,k
            
!       Normal X components
            do j=2,ny-1
                  do k = 2,nz-1
                        b(2,j,k,1) = dx*(b(2,j+1,k,2) - b(2,j,k,2))/dy + &
                              dx*(b(2,j,k+1,3) - b(2,j,k,3))/dz_grid(k) + &
                              b(3,j,k,1)
                              
                        b(nx-1,j,k,1) = b(nx-2,j,k,1) - &
                              dx*(b(nx-1,j+1,k,2) - b(nx-2,j,k,2))/dy - &
                              dx*(b(nx-2,j,k+1,3) - b(nx-2,j,k,3))/dz_grid(k)
                  enddo
            enddo
            
!       normal y components
!            do i=2,nx-1
!                  do k=2,nz-1
!                        b(i,2,k,2) = dy*(b(i+1,2,k,1) - b(i,2,k,1))/dx + &
!                              dy*(b(i,2,k+1,3) - b(i,2,k,3))/dz_grid(k) + &
!                              b(i,3,k,2)

!                        b(i,ny,k,2) = b(i,ny-1,k,2) - &
!                              dy*(b(i+1,ny-1,k,1) - b(i,ny-1,k,1))/dx - &
!                              dy*(b(i,ny-1,k+1,3) - b(i,ny-1,k,3))/dz_grid(k)
!                  enddo
!            enddo

!       normal z components
!            do i=2,nx-1
!                  do j=2,ny-1
!                        b(i,j,2,3) = dz_grid(2)*(b(i+1,j,2,1) - b(i,j,2,1))/dx + &
!                              dz_grid(2)*(b(i,j+1,2,2) - b(i,j,2,2))/dy + &
!                              b(i,j,3,3)

!                        b(i,j,nz,3) = b(i,j,nz-1,3) - &
!                              dz_grid(nz)*(b(i+1,j,nz-1,1) - b(i,j,nz-1,1))/dx + &
!                              dz_grid(nz)*(b(i,j+1,nz-1,2) - b(i,j,nz-1,2))/dy

!                  enddo
!            enddo
      
      end subroutine fix_normal_b
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      subroutine smooth_boundary(b)
            use dimensions
            implicit none
            real, intent(inout):: b(nx,ny,nz,3)
            integer:: i,j,k,m
            
!       X surfaces
            do j=2,ny-1
                  do k=2,nz-1
                        do m=1,3
                              b(1,j,k,m) = b(2,j,k,m)
                              b(nx,j,k,m) = b(nx-1,j,k,m)
                        enddo
                  enddo
            enddo
            
!       Y surfaces
            do i=1,nz
                  do k=1,nz
                        do m = 1,3
                              b(i,1,k,m) = b(i,2,k,m)
                              b(i,ny,k,m) = b(i,ny-1,k,m)
                        enddo
                  enddo
            enddo
            
!       Z surfaces
            do i =1,nx
                  do j= 1,ny
                        do m=1,3
                              b(i,j,1,m) = b(i,j,2,m)
                              b(i,j,nz,m) = b(i,j,nz-1,m)
                        enddo
                  enddo
            enddo
            
      end subroutine smooth_boundary
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine fix_tangential_E(E)
            use dimensions
            implicit none
            real, intent(inout):: E(nx,ny,nz,3)
            integer:: j,k
            
            do j=2,ny   !periodic boundary conditions
                  do k=2,nz
                        E(nx,j,k,2) = -2.3
                        E(nx,j,k,3) = 0.0
                  enddo
            enddo
            
      end subroutine fix_tangential_E
            
end module boundary

