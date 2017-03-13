module initial
      use dimensions
      implicit none
      save
      contains
      
      
      subroutine grd6_setup(b0,bt,b12,b1,b1p2,nu,input_Eb)
            use inputs, only: q, mO, PI, b0_top, b0_bottom, b0_init, nu_init, km_to_m, mu0
            use grid, only: dx_cell, dy_cell, dz_cell
            implicit none
            real, intent(out):: b0(nx,ny,nz,3), &
                                bt(nx,ny,nz,3), &
                                b12(nx,ny,nz,3), &
                                b1(nx,ny,nz,3), &
                                b1p2(nx,ny,nz,3), &
                                nu(nx,ny,nz)
            real, intent(inout):: input_Eb                    
                                
            real:: eoverm, mO_q, vol
            real:: b0_1x, b0_2x, b0_1y, b0_2y, phi
            integer:: i,j,k,m
            
            eoverm = q/mO
            mO_q = mO/q
            
            phi = 2.0*PI/180.0
            
            b0_1x = b0_top*eoverm*sin(phi)
            b0_2x = b0_bottom*eoverm*sin(phi)
            
            b0_1y = -b0_top*eoverm*cos(phi)
            b0_2y = -b0_bottom*eoverm*cos(phi)
            
            do i=1,nx
                  do j=1,ny
                        do k=1,nz
                              b0(i,j,k,1) = 0.0
                              b0(i,j,k,2) = 0.0
                              b0(i,j,k,3) = b0_init*eoverm
                        enddo
                  enddo
            enddo
            
!            input_Eb = 0.0
            do i=1,nx
                  do j=1,ny
                        do k= 1,nz
                              nu(i,j,k) = nu_init !+ 0.05*(tanh(real(k-570))-tanh(real(k-30))+2)
                              do m = 1,3
                                    bt(i,j,k,m) = b0(i,j,k,m)
                                    b12(i,j,k,m) = b0(i,j,k,m)
                                    b1(i,j,k,m) = b0(i,j,k,m)
                                    b1p2(i,j,k,m) = b0(i,j,k,m)
                                    vol = dx_cell(i)*dy_cell(j)*dz_cell(k)*km_to_m**3
                                    input_Eb = input_Eb + (vol/(2.0*mu0))*(mO_q*b0(i,j,k,m))**2
                              enddo
                        enddo
                  enddo
            enddo
            
            
!            open(40,file='b0.dat',status='unknown',form='unformatted')
!            write(40) nz
!            write(40) b0
!            close(40)
            
      end subroutine grd6_setup
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 
      subroutine grd7()
            use grid
            use mult_proc, only: my_rank
            use inputs, only: dx,dy,delz,out_dir,zsf!,nrgrd
            implicit none
            integer, parameter:: nrgrd = 0
            integer:: i,j,k,ind
            real:: xsf,zplus,zminus,xplus,xminus,yplus,yminus
            
            rk = nz/2
            rj= ny/2
            ri = nx/2
!!!!!!!!!Unstretched grids!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
            do j=1,ny
                  qy(j) = j*dy
                  dy_grid(j) = dy
            enddo
            
            do i = 1,nx
                  qx(i) = i*dx
                  dx_grid(i) = dx
            enddo

            
!!!!!!!Stretch X direction!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!            xsf = 0.0   !x stretch factor
!            !up from center
!            do i = ri, ri+nrgrd
!                  dx_grid(i) = dx
!            enddo
!            do i = ri + nrgrd + 1, nx
!                  dx_grid(i) = dx + xsf*dx*(i-(ri+nrgrd+1))/(nx-(ri+nrgrd+1))
!            enddo
!            !down from center
!            do i=ri-nrgrd, ri-1
!                  dx_grid(i) = dx
!            enddo
!            do i = 1, ri-nrgrd-1
!                  ind = ri-nrgrd-i
!                  dx_grid(ind) = dx + xsf*dx*(ri-nrgrd-1-ind)/(ri-nrgrd-1)
!            enddo
            
!            qx(1) = dx
            
!            do i = 2,nx
!                  qx(i) = qx(i-1)+dx_grid(i)
!            enddo
            
!            do i=1, nx-1
!                  dx_grid(i) = qx(i+1)-qx(i)
!            enddo
            
!            dx_grid(nx) = dx_grid(nx-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
            
            
!!!!!!!Stretch Z direction!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
!            zsf = 0.0   !z stretch factor
            !up from center
            
            do k = rk, rk + nrgrd
                  dz_grid(k) = delz
            enddo
                      
            
            do k = rk + nrgrd + 1, nz
                  dz_grid(k) = delz + zsf*delz*(k-(rk+nrgrd+1))/(nz-(rk+nrgrd+1))
            enddo
            
            !down from center
            do k = rk - nrgrd, rk - 1
                  dz_grid(k) = delz
            enddo
            do k = 1, rk - nrgrd - 1
                  ind = rk-nrgrd-k
                  dz_grid(ind) =  delz + zsf*delz*(rk-nrgrd-1-ind)/(rk-nrgrd-1)
            enddo
            
            qz(1) = delz
            
            do k=2,nz
                  qz(k) = qz(k-1)+dz_grid(k)
            enddo
            
            do k=1, nz-1
                  dz_grid(k) = qz(k+1)-qz(k)
            enddo
            
            dz_grid(nz) = dz_grid(nz-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            dz_cell(1) = dz_grid(1)
            dz_cell(nz) = dz_grid(nz)
            zrat(1) = 0.5
            zrat(nz) = 0.5
            do k=2, nz-1
                  dz_cell(k) = ((qz(k+1) + qz(k))/2.0) - ((qz(k) + qz(k-1))/2.0)
                  zplus = (qz(k+1) + qz(k))/2.0
                  zminus = (qz(k) + qz(k-1))/2.0
                  zrat(k) = (qz(k) - zminus)/(zplus - zminus)
            enddo
            
            dx_cell(1) = dx_grid(1)
            dx_cell(nx) = dx_grid(nx)
            xrat(1) = 0.5
            xrat(nx) = 0.5
            do i=2, nx-1
                  dx_cell(i) = ((qx(i+1) + qx(i))/2.0) - ((qx(i) + qx(i-1))/2.0)
                  xplus = (qx(i+1) + qx(i))/2.0
                  xminus = (qx(i) + qx(i-1))/2.0
                  xrat(i) = (qx(i) - xminus)/(xplus - xminus)
            enddo
            
            dy_cell(1) = dy_grid(1)
            dy_cell(ny) = dy_grid(ny)
            yrat(1) = 0.5
            yrat(ny) = 0.5
            do j=2, ny-1
                  dy_cell(j) = ((qy(j+1) + qy(j))/2.0) - ((qy(j) + qy(j-1))/2.0)
                  yplus = (qy(j+1) + qy(j))/2.0
                  yminus = (qy(j) + qy(j-1))/2.0
                  yrat(j) = (qy(j) - yminus)/(yplus - yminus)
            enddo
            
            if (my_rank .eq. 0) then
                  open(40,file=trim(out_dir)//'c.coord.dat',status='unknown',form='unformatted')
                  
                  write(40) nx
                  write(40) ny
                  write(40) nz
                  write(40) qx
                  write(40) qy
                  write(40) qz
                  write(40) dz_grid
                  write(40) dz_cell
                  close(40)
                  
            endif
            
      end subroutine grd7
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      subroutine grid_gaussian()
! Generates a grid that scales the z axis so that the spacing between grid points is a multiple of lambda_i
! using an analytical gaussian.
            use grid
            use mult_proc, only: my_rank
            use inputs, only: dx,dy,delz,out_dir,zsf,nrgrd,grad,amp,q,mion,nf_init
            implicit none
            integer:: i,j,k,ind
            real:: xsf,zplus,zminus,xplus,xminus,yplus,yminus
            
            rk = nz/2
            rj= ny/2
            ri = nx/2
!!!!!!!!!Unstretched grids!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
            do j=1,ny
                  qy(j) = j*dy
                  dy_grid(j) = dy
            enddo
            
            do i = 1,nx
                  qx(i) = i*dx
                  dx_grid(i) = dx
            enddo
            
!!!!!!!!!Stretch z!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            qz(nz/2) = 0
            do k=rk+1, nz
                  qz(k) = qz(k-1) + zsf*3e8/1e3*sqrt(8.85e-12*mion/(q*q* &
                  (nf_init/1e9+nf_init/1e9*(amp-1.0)*exp(-(qz(k-1)/(grad*delz))**2))))
            enddo
            
!            do k = 1, rk-1
!                  ind = rk-k
!                  qz(ind) = qz(ind+1) - zsf*3e8/1e3*sqrt(8.85e-12*mion/(q*q* &
!                  (nf_init/1e9+nf_init/1e9*(amp-1.0)*exp(-(qz(ind+1)/(grad*delz))**2))))
!            enddo

            do k = 1 , rk-1
                  ind = rk - k
                  qz(ind) = -qz(rk+k)
            enddo
            zplus = qz(nz-1)
            do k=1,nz
                  qz(k) = qz(k) + zplus
            enddo
            
            do k=1, nz-1
                  dz_grid(k) = qz(k+1)-qz(k)
            enddo
            
            dz_grid(nz) = dz_grid(nz-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            dz_cell(1) = dz_grid(1)
            dz_cell(nz) = dz_grid(nz)
            zrat(1) = 0.5
            zrat(nz) = 0.5
            do k=2, nz-1
                  dz_cell(k) = ((qz(k+1) + qz(k))/2.0) - ((qz(k) + qz(k-1))/2.0)
                  zplus = (qz(k+1) + qz(k))/2.0
                  zminus = (qz(k) + qz(k-1))/2.0
                  zrat(k) = (qz(k) - zminus)/(zplus - zminus)
            enddo
            
            dx_cell(1) = dx_grid(1)
            dx_cell(nx) = dx_grid(nx)
            xrat(1) = 0.5
            xrat(nx) = 0.5
            do i=2, nx-1
                  dx_cell(i) = ((qx(i+1) + qx(i))/2.0) - ((qx(i) + qx(i-1))/2.0)
                  xplus = (qx(i+1) + qx(i))/2.0
                  xminus = (qx(i) + qx(i-1))/2.0
                  xrat(i) = (qx(i) - xminus)/(xplus - xminus)
            enddo
            
            dy_cell(1) = dy_grid(1)
            dy_cell(ny) = dy_grid(ny)
            yrat(1) = 0.5
            yrat(ny) = 0.5
            do j=2, ny-1
                  dy_cell(j) = ((qy(j+1) + qy(j))/2.0) - ((qy(j) + qy(j-1))/2.0)
                  yplus = (qy(j+1) + qy(j))/2.0
                  yminus = (qy(j) + qy(j-1))/2.0
                  yrat(j) = (qy(j) - yminus)/(yplus - yminus)
            enddo
            
            if (my_rank .eq. 0) then
                  open(40,file=trim(out_dir)//'c.coord.dat',status='unknown',form='unformatted')
                  
                  write(40) nx
                  write(40) ny
                  write(40) nz
                  write(40) qx
                  write(40) qy
                  write(40) qz
                  write(40) dz_grid
                  write(40) dz_cell
                  close(40)
                  
            endif

      end subroutine grid_gaussian

      
end module initial
            