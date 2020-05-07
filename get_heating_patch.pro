;---------------------------------------------------------------------------
pro plot_image,img,x,y,sclf,loc,ii,kk,nxz,plt_tit
;---------------------------------------------------------------------------

   im = image(img,x,y,/current,rgb_table=33,layout=[1,1,loc],font_size=12)

   xax = axis('x',axis_range=[min(x),max(x)],location=[min(y)],thick=2,tickdir=1,target=im,tickfont_size=12)
   yax = axis('y',axis_range=[min(y),max(y)],location=[min(x)],thick=2,tickdir=1,target=im,tickfont_size=12)

   im.xtitle='$x (c/\omega_{pi})$'
   im.ytitle='$z (c/\omega_{pi})$'
   whx = where(x gt 0)
   why = where(y gt 0)
   im.xrange=[min(x(whx)),max(x)]
   im.yrange=[min(y(why)),max(y)]
   im.title = plt_tit
   
;   ct = colorbar(target = im,title=tit,orientation=1,textpos=1,font_size=12);,$                                                                      
;                 position=[max(x), min(y),                                                                                                           
;                 0.05*(max(x)-min(x)),max(y)],/data)                                                                                                 
   im.scale,sclf,sclf
   ct = colorbar(target = im,title='$B_\perp$ (nT$^2$)',textpos=1,font_size=12,orientation=1,$
                 position=[1.01,0,1.06,1.0],/relative)

   im = plot([x(ii-nxz/2),x(ii+nxz/2)],[y(kk-nxz/2),y(kk-nxz/2)],'4w',/overplot)
   im = plot([x(ii-nxz/2),x(ii+nxz/2)],[y(kk+nxz/2),y(kk+nxz/2)],'4w',/overplot)
   im = plot([x(ii-nxz/2),x(ii-nxz/2)],[y(kk-nxz/2),y(kk+nxz/2)],'4w',/overplot)
   im = plot([x(ii+nxz/2),x(ii+nxz/2)],[y(kk-nxz/2),y(kk+nxz/2)],'4w',/overplot)
   
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
pro get_particle_heating,info,s,poft,dEdt_1,p1,tm
;---------------------------------------------------------------------------      
  
  for i = 1,info.nfr do begin
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
  w1=window()
  dEdt_1 = (shift(poft_arr,-1)-shift(poft_arr,1))/(2*info.dt*info.t0)/1e-15
  p1 = plot(tm(2:*)*s.Omega_i,(dEdt_1(2:*))>0,'2r-D',/current,NAME='q')
  p1.sym_filled=1
  p1.ytitle='Heating rate density [$10^{-15}$ W/m$^{3}$]'
  p1.xtitle='time ($\Omega_i^{-1}$)'
  p1.save,'heat_rate.png'
  p1.font_size=18
  
end
;---------------------------------------------------------------------------      

;w = window()
@clr_win
@get_const

;dir = './run_va_0.8_beta_1/'
dir = './run_va_0.8_beta_3/'

nfr = 16   ;number of frames.
nxz = 12   ;fft domain

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

info = {info, nfr:nfr, dt:dt, t0:200., dir:dir, nxz:nxz}
s = {plasma_params, Omega_i:Omega_i, wpi:wpi, $
              cwpi:cwpi, beta:beta, rhoi:rhoi, kp_rhoi:(1/sqrt(beta))*2*!pi/cwpi }

dx = x(1)-x(0)

poft = 0.0
;get_particle_heating,info,s,poft,dEdt_1,p1,tm

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
for j = 1,nfr do begin
   nfrm = j
   print,'nfrm....',nfrm
   c_read_3d_vec_m_32,dir,'c.b1',nfrm,b1
   c_read_3d_m_32,dir,'c.np',nfrm,np
endfor
b1 = b1*mproton/q
cnt = 0
pmhd_arr = 0
pkaw_arr = 0

cnt = 0
pwr = 0.
pwrk = 0.
psd_k_mhd = 0.
psd_k_kaw = 0.

loadct,39
tvscl,b1(*,ny/2,*,0)^2 + b1(*,ny/2,*,2)^2
;plot_image,reform(b1(*,ny/2,*,0)^2 + b1(*,ny/2,*,2)^2),indgen(n_elements(x)),indgen(n_elements(z)),1,1,'string(nxz)',nxz

kk = 0
cnt = 0
w2 = window(dimensions=[900,600])
while (kk lt nz-1) do begin
   cursor,ii,kk,/device,/down
   print,ii,kk
;   plot_image,reform(sqrt(b1(*,ny/2,*,0)^2 + b1(*,ny/2,*,2)^2)/1e-9),indgen(n_elements(x)),$
;              indgen(n_elements(z)),0.8,1,ii,kk,nxz


;stop
   
;   for jj = 2, ny-2 do begin
   jj = ny/2
;   x0 = nxz/2 + 1
;   while((x0 + nxz/2) lt nx-2) do begin
   
   x1 = ii-nxz/2
   x2 = ii+nxz/2-1
   z1 = kk-nxz/2
   z2 = kk+nxz/2-1
   
   fx = fft(reform(b1(x1:x2,jj,z1:z2,0)))
   fz = fft(reform(b1(x1:x2,jj,z1:z2,2)))
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
         fftpwr(i,k) = (abs(fx(i,k))^(2.) + abs(fz(i,k))^(2.)) 
         
         if ((k_perp lt s.kp_rhoi) and (k_perp gt kx(0))) then pwr_mhd(i,k) = $
            (psd(i,k)*k_perp)/sqrt(muo^3*rho)
         if (k_perp ge s.kp_rhoi) then pwr_kaw(i,k) = (k_perp*rhoi)*(psd(i,k)*k_perp)/sqrt(muo^3*rho)
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
   print,pmhd
    
   wh = where((pwr_arr_kaw gt 0) and (kx/s.kp_rhoi ge 1.0))
   pkaw = mean(pwr_kaw_arr_sum(wh))
   print,pkaw

   cnt = cnt+1

   w2.erase
   plt_tit = '$q_{MHD}$: '+string(pmhd/1e-15)+' ($10^{-15}$ W/m$^3$);    $q_{KAW}$: '+string(pkaw/1e-15)+' ($10^{-15}$ W/m$^3$)
   plot_image,reform(sqrt(b1(*,ny/2,*,0)^2 + b1(*,ny/2,*,2)^2)/1e-9),x*1e3/s.cwpi,z*1e3/s.cwpi,$
              1,1,ii,kk,nxz,plt_tit
   
endwhile

end
