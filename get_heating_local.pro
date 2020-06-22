;---------------------------------------------------------------------------
pro plot_image,img,x,y,sclf,loc,ii,kk,nxz,plt_tit,l_tit
;---------------------------------------------------------------------------

   ctbl = colortable(15,/reverse)
   im = image(img,x,y,/current,rgb_table=ctbl,layout=[1,1,loc],font_size=16,$
              aspect_ratio=1.0)

   xax = axis('x',axis_range=[min(x),max(x)],location=[min(y)],thick=1,tickdir=1,target=im,tickfont_size=16)
   yax = axis('y',axis_range=[min(y),max(y)],location=[min(x)],thick=1,tickdir=1,target=im,tickfont_size=16)

   im.xtitle='$x (c/\omega_{pi})$'
   im.ytitle='$z (c/\omega_{pi})$'
   whx = where(x gt 0)
   why = where(y gt 0)
   im.xrange=[min(x(whx)),max(x)]
   im.yrange=[min(y(why)),max(y)]
   im.title = plt_tit
   im.position=[0.2,0.2,0.8,0.8]
   im.font_size=16
;   ct = colorbar(target = im,title=tit,orientation=1,textpos=1,font_size=12);,$                                                                      
;                 position=[max(x), min(y),                                                                                                           
;                 0.05*(max(x)-min(x)),max(y)],/data)                                                                                                 
   im.scale,sclf,sclf
   ct = colorbar(target = im,title=l_tit,textpos=1,font_size=14,orientation=1,$
                 position=[1.01,0.0,1.05,1.0],/relative)

;   im = plot([x(ii-nxz/2),x(ii+nxz/2)],[y(kk-nxz/2),y(kk-nxz/2)],'4w',/overplot)
;   im = plot([x(ii-nxz/2),x(ii+nxz/2)],[y(kk+nxz/2),y(kk+nxz/2)],'4w',/overplot)
;   im = plot([x(ii-nxz/2),x(ii-nxz/2)],[y(kk-nxz/2),y(kk+nxz/2)],'4w',/overplot)
;   im = plot([x(ii+nxz/2),x(ii+nxz/2)],[y(kk-nxz/2),y(kk+nxz/2)],'4w',/overplot)
   
end
;---------------------------------------------------------------------------

;---------------------------------------------------------------------------
pro plot_psd,kx,s,psd_arr_sum,psd_k_mhd,psd_k_kaw
;---------------------------------------------------------------------------
  
  wh = where((kx gt kx(0)) and (kx le s.kp_rhoi))
  fkx_mhd = poly_fit(alog(kx(wh)),alog(psd_arr_sum(wh)),1)
  psd_k_mhd = [psd_k_mhd,fkx_mhd(1)]
  
;  plot,alog(kx(1:*)),alog(psd_arr_sum(1:*)),charsize=2 ;,/ylog,/xlog
;  oplot,[(1/sqrt(beta))*2*!pi/cwpi,(1/sqrt(beta))*2*!pi/cwpi],[min(psd_arr_sum),max(psd_arr_sum)],linestyle=1
;  oplot,alog(kx(wh)),fkx_mhd(0) + fkx_mhd(1)*alog(kx(wh))
;  oplot,alog(kx(wh)),fkx_mhd(0) - (5./3.)*alog(kx(wh))
  
  wh = where((kx ge s.kp_rhoi))
  fkx_kaw = poly_fit(alog(kx(wh)),alog(psd_arr_sum(wh)),1)
  psd_k_kaw = [psd_k_kaw,fkx_kaw(1)]
  
;  oplot,alog(kx(wh)),fkx_kaw(0) + fkx_kaw(1)*alog(kx(wh))

end
;---------------------------------------------------------------------------      

;---------------------------------------------------------------------------      
pro plot_q,kx,beta,cwpi,pwr_arr_sum,pwr_kaw_arr_sum
;---------------------------------------------------------------------------      
  
  plot,kx/((1/sqrt(beta))*2*!pi/cwpi),pwr_arr_sum/1e-15,psym=10,charsize=2,/ylog,yrange=[0.01,100]
  oplot,kx/((1/sqrt(beta))*2*!pi/cwpi),pwr_kaw_arr_sum/1e-15,psym=10,thick=2

end
;---------------------------------------------------------------------------      


;---------------------------------------------------------------------------      
pro get_particle_heating,info,s,poft,dEdt_1,p1,tm,w4
;---------------------------------------------------------------------------      
  
  for i = info.nfr0,info.nfr1 do begin
     nfrm = i
     
     c_read_3d_m_32,info.dir,'c.temp_p',nfrm,tp
     c_read_3d_m_32,info.dir,'c.np',nfrm,np
     
     p = (3./2.)*np*tp
     
     parr = reform(p(*,*,*))
     parr1 = reform(p(*,*,*))
     
     sz = size(parr)
     print,total(parr1)/(float(sz(1))*float(sz(3))*float(sz(2)))
     poft = [poft,total(parr1)/(float(sz(1))*float(sz(3))*float(sz(2)))]
     
  endfor
  
  poft_arr = poft*1.6e-19/1e9
  tm = info.dt*info.t0*findgen(n_elements(poft_arr))
  w4=window()
  dEdt_1 = (shift(poft_arr,-1)-shift(poft_arr,1))/(2*info.dt*info.t0)/1e-15
;  dEdt_1 = (shift(poft_arr,-1)-poft_arr)/(info.dt*info.t0)/1e-15
  p1 = plot(tm(2:*)*s.Omega_i,(dEdt_1(2:*))>0,'2r-D',/current,NAME='q')
  p1.sym_filled=1
  p1.ytitle='Heating rate density [$10^{-15}$ W/m$^{3}$]'
  p1.xtitle='time ($\Omega_i^{-1}$)'
  p1.save,'heat_rate.png'
  p1.font_size=18
  
end
;---------------------------------------------------------------------------


;---------------------------------------------------------------------------
pro get_curl,bx,bz,dx,curlb
;---------------------------------------------------------------------------
  sz = size(bx)
  curlb = fltarr(sz(1),sz(2))
  for i = 1,sz(1)-2 do begin
     for k = 1,sz(2)-2 do begin
        curlb(i,k) = -(bz(i+1,k)-bz(i-1,k))/dx + (bx(i,k+1) - bx(i,k-1))/dx
     endfor
  endfor
  curlb(*,sz(2)-1) = curlb(*,sz(2)-3)
  curlb(*,sz(2)-2) = curlb(*,sz(2)-3)
end
;---------------------------------------------------------------------------


;---------------------------------------------------------------------------
function geo_mean,a
  return, 10^(mean(alog10(a)))
end
;---------------------------------------------------------------------------

;function dotp,a,b
;  return, a(0)*b(0) + a(1)*b(1) + a(2)*b(2)
;end

;;---------------------------------------------------------------------------
;pro get_bperp,b1,b2,b3
;;---------------------------------------------------------------------------
;  b2 = fltarr(3)
;  b2(1) = 0                     ;choose y compondent to be zero for b2
;  b2(2) = 1.0
;  b2(0) = -b1(2)/b1(0)
;  b2hat = b2/sqrt(dotp(b2,b2))
;  b3 = crossp(b1,b2)
;  b3hat = b3/sqrt(dotp(b3,b3))
;  b2 = dotp(b1,b2hat)*b2hat
;  b3 = dotp(b1,b3hat)*b3hat
;  print,'b1...',b1
;  print,'b2...',b2
;  print,'b3...',b3
;;  if (dotp(b2,b3) gt 1e-10) then print,'err'
;end
;;---------------------------------------------------------------------------      

;w = window()
@clr_win
@get_const

;dir = './run_va_0.8_beta_1/'
dir = './run_va_0.8_beta_3/'

nfr0 = 18
nfr1 = 19                        ;number of frames.
nxz = 6   ;fft domain

;initialize
read_para,dir
restore,filename=dir+'para.sav'
read_coords,dir,x,y,z

beta = 3.0
Omega_i = q*b0_top/mproton
wpi = sqrt(q*q*np_top/1e9/(epo*mproton))
cwpi = 3e8/wpi
cwpi = cwpi
rhoi = sqrt(beta)*cwpi

info = {info, nfr0:nfr0, nfr1:nfr1, dt:dt, t0:200., dir:dir, nxz:nxz}
s = {plasma_params, Omega_i:Omega_i, wpi:wpi, $
              cwpi:cwpi, beta:beta, rhoi:rhoi, kp_rhoi:(1/sqrt(beta))*2*!pi/cwpi }

dx = x(1)-x(0)

poft = 0.0
get_particle_heating,info,s,poft,dEdt_1,p1,tm,w4

barr = fltarr(nx)
pmhd = 0.0
pkaw = 0.0
pmhd_sd = 0.0
pkaw_sd = 0.0
psd_k = 0.


psd_av_mhd = 0.0
psd_sd_mhd = 0.0
psd_av_kaw = 0.0
psd_sd_kaw = 0.0

window,1
device,decomposed=0
w1 = window(dimensions=[900,600])
w2 = window(dimensions=[900,600])
w3 = window(dimensions=[900,600])
w5 = window(dimensions=[900,600])
w6 = window(dimensions=[900,600])
qmhd = 0.
qkaw = 0.

for j = info.nfr0,info.nfr1 do begin
   qmhd_1 = 0.
   qkaw_1 = 0.
;   for jj = 1,ny-2,10 do begin
   jj = ny/2
      nfrm = j
      print,'nfrm....',nfrm
      c_read_3d_vec_m_32,dir,'c.b1',nfrm,b1
      c_read_3d_vec_m_32,dir,'c.up',nfrm,up
      c_read_3d_m_32,dir,'c.np',nfrm,np
      
      b1 = b1*mproton/q
      cnt = 0
      pmhd_arr = 0
      pkaw_arr = 0
      
;b2p = b1
;b3p = b1
;for i = 0,nx-1 do begin
;   for j = 0,ny-1 do begin
;      for k = 0,nz-1 do begin
;         get_bperp,reform(b1(i,j,k,*)),b2,b3
;         b2p(i,j,k) = sqrt(dotp(b2,b2))
;         b3p(i,j,k) = sqrt(dotp(b3,b3))
;      endfor
;   endfor
;endfor
;
;!p.multi=[0,1,2]
;bperp = reform(sqrt(b1(*,jj,*,0)^2 + b1(*,jj,*,2)^2))
;surface,bperp,charsize=5
;surface,reform(sqrt(b2p(*,jj,*)^2 + b3p(*,jj,*)^2)),charsize=5
      
      cnt = 0
      pwr = 0.
      pwrk = 0.
      psd_k_mhd = 0.
      psd_k_kaw = 0.
      
;loadct,39
;tvscl,b1(*,ny/2,*,0)^2 + b1(*,ny/2,*,2)^2
;plot_image,reform(b1(*,ny/2,*,0)^2 + b1(*,ny/2,*,2)^2),indgen(n_elements(x)),indgen(n_elements(z)),1,1,'string(nxz)',nxz
      
      qmhd_arr = fltarr(nx,nz)
      qkaw_arr = fltarr(nx,nz)
      kk = 0
      cnt = 0
;   jj = ny/2+60
      for ii = nxz/2+1, nx-nxz/2-2 do begin
         for kk = nxz/2+1, nz-nxz/2-2 do begin
;while (kk lt nz-1) do begin
;   cursor,ii,kk,/device,/down
;   print,ii,kk
;   plot_image,reform(sqrt(b1(*,ny/2,*,0)^2 + b1(*,ny/2,*,2)^2)/1e-9),indgen(n_elements(x)),$
;              indgen(n_elements(z)),0.8,1,ii,kk,nxz
            
            
;stop
            
;   for jj = 2, ny-2 do begin
;   jj = ny/2
;   x0 = nxz/2 + 1
;   while((x0 + nxz/2) lt nx-2) do begin
            
            x1 = ii-nxz/2
            x2 = ii+nxz/2-1
            z1 = kk-nxz/2
            z2 = kk+nxz/2-1
            
            fx = fft(reform(b1(x1:x2,jj,z1:z2,0)))
            fz = fft(reform(b1(x1:x2,jj,z1:z2,2)))
;      fx = fft(b2p)
;      fz = fft(b3p)
            fperp = fft(sqrt(reform(b1(x1:x2,jj,z1:z2,0)^2) + reform(b1(x1:x2,jj,z1:z2,2)^2)))
            
            density = mean(reform(np(x1:x2,jj,z1:z2)))
            sz = size(fx)
            N = sz(1)
            dx = (x(1)-x(0))*1e3
            kx = indgen(N)
            N21x = N/2+1
            kx(N21x)= N21x - N + findgen(N21x-2)
            kx = 2*!pi*kx/(N*dx)
            
            N = sz(2)
            dx = (x(1)-x(0))*1e3
            kz = indgen(N)
            N21z = N/2+1
            kz(N21z)= N21z - N + findgen(N21z-2)
            kz = 2*!pi*kz/(N*dx)
            
            psd = fltarr(n_elements(kx),n_elements(kz))
            fftpwr = psd
            pwr_mhd = fltarr(n_elements(kx),n_elements(kz))
            pwr_kaw = pwr_mhd
            
            for i = 0,n_elements(kx)-1 do begin
               for k = 0,n_elements(kz)-1 do begin
                  k_perp = sqrt(kx(i)^2 + kz(k)^2)
                  rho = mp*density/1e9
                  psd(i,k) = (abs(fx(i,k))^(2.) + abs(fz(i,k))^(2.))^(3./2.)
;               psd(i,k) = (abs(fx(i,k))^(3.) + abs(fz(i,k))^(3.)) 
                  fftpwr(i,k) = (abs(fx(i,k))^(2.) + abs(fz(i,k))^(2.)) 
                  
                  if ((k_perp lt s.kp_rhoi) and (k_perp gt kx(0))) then pwr_mhd(i,k) = $
                     (psd(i,k)*k_perp)/sqrt(muo^3*rho)
;               if ((k_perp ge s.kp_rhoi) and (k_perp gt kx(0))) then pwr_kaw(i,k) = (k_perp*rhoi)*(psd(i,k)*k_perp)/sqrt(muo^3*rho)
                  ;if (k_perp ge s.kp_rhoi) then pwr_kaw(i,k) = (sqrt(0.75)*k_perp*rhoi)*(psd(i,k)*k_perp)/sqrt(muo^3*rho)
                  if (k_perp gt kx(0)) then begin
                     b = double(k_perp^2*rhoi^2)
                     gamma0 = beseli(b,0)*exp(-b)
                     gamma1 = beseli(b,1)*exp(-b)
                                ;print,'gamma...',b,gamma0,gamma1
                     ;original formulation
                     ;pwr_kaw(i,k) = sqrt(1+(0.75*k_perp^2*rhoi^2))*(psd(i,k)*k_perp)/sqrt(muo^3*rho)
                                ; modified for no flow energy in wave
                                ;factor to reduce wave energy for KAW
                     KAW_Walen = 0.5+0.5*(1./(1. + k_perp^2*rhoi^2))*(1./(1. + 1.25*k_perp^2*rhoi^2))^2
                     pwr_kaw(i,k) = KAW_Walen*sqrt(1+(0.75*k_perp^2*rhoi^2))*(psd(i,k)*k_perp)/sqrt(muo^3*rho)
                     ;pwr_kaw(i,k) = sqrt((1 - gamma0)/b)*((gamma0-gamma1)*b/(1 - gamma0))*(psd(i,k)*k_perp)/sqrt(muo^3*rho)
                     ;pwr_kaw(i,k) = (psd(i,k)*k_perp)/sqrt(muo^3*rho)/(1 + 1.25*b)
                  endif
               endfor
            endfor
            
            pwr_arr = pwr_mhd(0:N21x-1,0:N21z-1)
            pwr_arr_kaw = pwr_kaw(0:N21x-1,0:N21z-1)
            
            kx = kx(0:N21x-1)
            kz = kz(0:N21z-1)
            
            psd_arr_sum = fltarr(n_elements(kx))
            pwr_arr_sum = fltarr(n_elements(kx))
            pwr_kaw_arr_sum = fltarr(n_elements(kx))
            dk = kx(1)-kx(0)
            for i = 0,n_elements(kx)-1 do begin
               for k = 0,n_elements(kz)-1 do begin
                  k_perp = sqrt(kx(i)^2 + kz(k)^2)
                  wh = where(abs(k_perp+0.0*dk - kx) eq min(abs(k_perp+0.0*dk-kx)))
                  pwr_arr_sum(wh(0)) = pwr_arr_sum(wh(0)) + pwr_arr(i,k)
                  pwr_kaw_arr_sum(wh(0)) = pwr_kaw_arr_sum(wh(0)) + pwr_arr_kaw(i,k)
                  psd_arr_sum(wh(0)) = psd_arr_sum(wh(0)) + fftpwr(i,k)
               endfor
            endfor
            
            wh = where((pwr_arr_sum gt 0) and (kx/s.kp_rhoi le 1.0))
            pmhd = mean(pwr_arr_sum(wh))
                                ;print,pmhd
            
            wh = where((pwr_arr_kaw gt 0) and (kx/s.kp_rhoi ge 1.0))
            pkaw = mean(pwr_kaw_arr_sum(wh))
                                ;print,pkaw
            
            cnt = cnt+1
            
            qmhd_arr(ii,kk) = pmhd
            qkaw_arr(ii,kk) = pkaw
            
         endfor
      endfor
      
      bperp = reform(sqrt(b1(*,jj,*,0)^2 + b1(*,jj,*,2)^2))
      tvscl,bperp

      bx = reform(b1(*,jj,*,0))
      bz = reform(b1(*,jj,*,2))
      ux = reform(up(*,jj,*,0))*1e3
      uz = reform(up(*,jj,*,2))*1e3
      rho = mp*reform(np(*,jj,*))/1e9
      
;   wh = where(qmhd_arr lt 1e-17)
;   qmhd_arr(wh)  = 0.0

;      qmhd_arr = median(qmhd_arr,8)
;      qkaw_arr = median(qkaw_arr,8)
      
      ctbl = colortable(3)
      cclr = [10,50,100,150,250]
      w1.erase
      w1.SetCurrent
      plt_tit = '$B_{in-plane}$'
      l_tit = '$|B|$ (nT)'
      plot_image,bperp/1e-9,x*1e3/s.cwpi,z*1e3/s.cwpi,$
                 1,1,ii,kk,nxz,plt_tit,l_tit
      c = contour(bperp,x*1e3/s.cwpi,z*1e3/s.cwpi,/overplot,n_levels=5,c_value=0,c_label_show=0,c_thick=1,$
                  rgb_table=ctbl,c_color=cclr)
      
;      ctbl = colortable(3,/reverse)      
      w5.erase
      w5.SetCurrent
      plt_tit = 'Els'+string(228b)+'sser vortices'
      l_tit = '$ln( |\omega^+||\omega^-|)$'
      get_curl,bx,bz,dx,curlb
      get_curl,ux,uz,dx,curlu
      muo=!pi*4e-7
      apx = ux + bx/sqrt(muo*rho)
      amx = ux - bx/sqrt(muo*rho)
      apz = uz + bz/sqrt(muo*rho)
      amz = uz - bz/sqrt(muo*rho)
      get_curl,amx,amz,dx,omegam
      get_curl,apx,apz,dx,omegap
      awx = -amz*omegap
      awz = amx*omegap
      get_curl,awx,awz,dx,curlaw
      ;get_curl,bx,bz,dx,curlb
      ;plot_image,smooth(abs(curlaw),2)>0.4,x*1e3/s.cwpi,z*1e3/s.cwpi,$
      ;           1,1,ii,kk,nxz,plt_tit,l_tit
      plot_image,smooth(alog((abs(omegam)*abs(omegap))>0.06),3),x*1e3/s.cwpi,z*1e3/s.cwpi,$
                 1,1,ii,kk,nxz,plt_tit,l_tit
      c = contour(bperp,x*1e3/s.cwpi,z*1e3/s.cwpi,/overplot,n_levels=5,c_value=0,c_label_show=0,c_thick=1,$
                 rgb_table=ctbl,c_color=cclr)
      
      w2.erase
      w2.SetCurrent
;   w2 = window(dimensions=[900,600])
      plt_tit = '$q_{MHD}$'
      l_tit = 'q ($10^{-15}$ W/m$^3$)'
      plot_image,qmhd_arr/1e-15,x*1e3/s.cwpi,z*1e3/s.cwpi,$
                 1,1,ii,kk,nxz,plt_tit,l_tit
      c = contour(bperp,x*1e3/s.cwpi,z*1e3/s.cwpi,/overplot,n_levels=5,c_value=0,c_label_show=0,c_thick=1,$
                 rgb_table=ctbl,c_color=cclr)
      
      w3.erase
      w3.SetCurrent
;   w3 = window(dimensions=[900,600])
      plt_tit = 'Heating rate density'
      plot_image,qkaw_arr/1e-15<200,x*1e3/s.cwpi,z*1e3/s.cwpi,$
                 1,1,ii,kk,nxz,plt_tit,l_tit
      c = contour(bperp,x*1e3/s.cwpi,z*1e3/s.cwpi,/overplot,n_levels=5,c_value=0,c_label_show=0,c_thick=1,$
                  rgb_table=ctbl,c_color=cclr)
; draw box
      x1 = x(10)*1e3/s.cwpi
      x2 = x(16)*1e3/s.cwpi
      z1 = z(10)*1e3/s.cwpi
      z2 = z(16)*1e3/s.cwpi
      p = plot([x1,x2],[z1,z1],/overplot)
      p = plot([x1,x2],[z2,z2],/overplot)
      p = plot([x1,x1],[z1,z2],/overplot)
      p = plot([x2,x2],[z1,z2],/overplot)

      w6.erase
      w6.SetCurrent
      wh = qkaw_arr(*,fix(nz/2))/1e-15 > 0.0
      p = plot(x*1e3/s.cwpi,qkaw_arr(*,fix(nz/2))/1e-15, ylog=1, /current,/xsty,/ysty,layout=[1,2,1])
      p.ytitle='q (10$^{-15}$ W/m$^3$)'
      p.position=[0.2,0.51,0.9,0.81]
      p.xshowtext=0
      p.font_size=14

      ;print,'z coord....',z(fix(nz/2))*1e3/s.cwpi
      ;stop
      p11 = plot(x(4:nx-5)*1e3/s.cwpi,b1(4:nx-5,fix(ny/2),fix(nz/2),0)/1e-9,'r', $
                 /current, /xsty,/ysty, layout=[1,2,2],NAME = '$B_x$')
      
      p12 = plot(x(4:nx-5)*1e3/s.cwpi,b1(4:nx-5,fix(ny/2),fix(nz/2),2)/1e-9,'g', $
                 /current, /xsty,/ysty, /overplot, layout=[1,2,2],NAME = '$B_z$')
      p12.position=[0.2,0.2,0.9,0.5]
      ;p13 = plot(x*1e3/s.cwpi,b1(*,fix(ny/2),fix(nz/2),1)/1e-9,'g', $
                                ;          /current, /xsty,/ysty,
                                ;          /overplot,
                                ;          layout=[1,2,2],NAME =
                                ;          '$B_y$')

      p11.ytitle='B (nT)'
     
      p11.xtitle='x (c/$\omega_{pi}$)'
      p11.font_size=14
      p14 = plot([x(4),x(nx-5)]*1e3/s.cwpi,[0,0],linestyle=1,/current,/overplot,layout=[1,2,2])
      l = legend(target=[p11,p12])
    
      ;p.xrange=[0.1,200]
      
      qkaw_back = geo_mean(qkaw_arr(*,0:20))/1e-15
      print,'qkaw_back...',qkaw_back
      
      qmhd_back = geo_mean(qmhd_arr(*,0:20))/1e-15
      print,'qmhd_back...',qmhd_back
      
      whkaw = where(qkaw_arr/1e-15 gt qkaw_back)
      print,'mean qkaw...',jj,mean(qkaw_arr(whkaw)/1e-15)
      
      whmhd = where(qmhd_arr/1e-15 gt qmhd_back)
      print,'mean qmhd...',jj,mean(qmhd_arr(whmhd)/1e-15)
      
;   w5.erase
;   w5.SetCurrent
;   h = histogram(qkaw_arr/1e-15,locations = bins,nbins=1000,min=0.1,max=100)
;   x0 = 0.1
;   x1= x0*2
;   xarr = x0
;   harr = 0
;   for i = 0,10 do begin
;      wh = where((bins gt x0) and (bins le x1))
;      harr = [harr,total(h(wh))]
;      xarr = [xarr,x1*2]
;      x0 = x1
;      x1 = x0*2
;   endfor
;   print,xarr,harr
;   b = plot(xarr(1:*),harr(1:*),'ro:',/xlog,/ylog,/current)
;;   b.xrange=[0.02,10]
;   b.xtitle='q ($10^{-15}$ W/m$^{-3}$)'

      ;background subtraction
;      qmhd_1 = [qmhd_1,mean(qmhd_arr(whmhd)/1e-15)]
;      qkaw_1 = [qkaw_1,mean(qkaw_arr(whkaw)/1e-15)]

      ;no background subtraction
      qmhd_1 = [qmhd_1,mean(qmhd_arr/1e-15)]
      qkaw_1 = [qkaw_1,mean(qkaw_arr/1e-15)]
      
;   endfor
   qmhd = [qmhd, mean(qmhd_1(1:*))]
   qkaw = [qkaw, mean(qkaw_1(1:*))]
;print,'mean qkaw 1...',qkaw
;print,'mean qmhd 1...',qmhd
endfor

qmhd = qmhd(1:*)
qkaw = qkaw(1:*)

;tm = info.dt*info.t0*findgen(n_elements(qmhd))
w4.SetCurrent
p1a = plot(tm*s.Omega_i,qmhd,'2b-tu',/overplot,NAME='$q_{MHD}$')
p1a.sym_filled=1
p1b = plot(tm*s.Omega_i,qkaw,'2g-s',/overplot,NAME='$q_{KAW}$')
p1b.sym_filled=1
p1.ylog=1
p1.yrange=[0.1,30]
p1.xrange=[90,250]
p1.title='$\beta$ = 3'
p1.position=[0.2,0.15,0.95,0.9]
p1.font_size=18
l1 = legend(target=[p1,p1a,p1b])
l1.font_size=18
;save,tm,s.Omega_i,qmhd,qkaw,filename='beta_3_8x8.sav'


end
