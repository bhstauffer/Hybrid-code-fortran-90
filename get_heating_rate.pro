;--------------------------------------------------------------------------
pro plot_image,img,x,y,sclf,loc,tit
;--------------------------------------------------------------------------
  
   im = image(img,x,y,/current,rgb_table=33,layout=[1,1,loc],font_size=12)

   xax = axis('x',axis_range=[min(x),max(x)],location=[0,min(y)],thick=2,tickdir=1,target=im,tickfont_size=12)
   yax = axis('y',axis_range=[min(y),max(y)],location=[0,0],thick=2,tickdir=1,target=im,tickfont_size=12)

   im.xtitle='$x (c/\omega_{pi})$'
   im.ytitle='$z (c/\omega_{pi})$'

;   ct = colorbar(target = im,title=tit,orientation=1,textpos=1,font_size=12);,$                                                                      
;                 position=[max(x), min(y),                                                                                                           
;                 0.05*(max(x)-min(x)),max(y)],/data)                                                                                                 
   im.scale,sclf,sclf
   ct = colorbar(target = im,title=tit,textpos=1,font_size=12,orientation=1,$
                 position=[1.01,0,1.06,1.0],/relative)

end
;--------------------------------------------------------------------------


;--------------------------------------------------------------------------
function geomean,arr
;--------------------------------------------------------------------------
  sz = size(arr)

  prod= 1.
  for i = 0,sz(1)-1 do begin
     prod = prod*arr(i)
  endfor
  geomean = prod^(1/float(sz(1)))
  return,geomean
end

;--------------------------------------------------------------------------

@clr_win
@get_const

dir = './run_va_0.8_beta_3/'
read_para,dir
restore,filename=dir+'para.sav'
read_coords,dir,x,y,z


wpi = sqrt(q*q*np_top/1e9/(epo*mproton))
cwpi = 3e8/wpi
;cwpi = cwpi/1e3

dx = (x(1)-x(0))*1e3

nfr = 25
poft = 0.0

w = window()
for i = 1,nfr do begin
   nfrm = i
   c_read_3d_vec_m_32,dir,'c.b1',nfrm,b1
;   c_read_3d_vec_m_32,dir,'c.up',nfrm,up
   c_read_3d_m_32,dir,'c.temp_p',nfrm,tp
   c_read_3d_m_32,dir,'c.np',nfrm,np

   b1 = b1*mproton/q
   
   p = np*tp
   
;   parr = reform(p(*,*,nz/2))
;   parr1 = reform(p(*,*,nz/2))

   parr = reform(p(*,*,*))
   parr1 = reform(p(*,*,*))

;   plot_image,reform(p(*,ny/2,*)),x,z,0.8,1,'p'
   
   sz = size(parr)
;   print,total(parr1)/(float(sz(1))*11.)
   poft = [poft,total(parr1)/(float(sz(1))*float(sz(3))*float(sz(2)))]
   
endfor

poft_arr = poft(4:*)*1.6e-19/1e9

tm = dt*200.*findgen(n_elements(poft_arr))

;p=plot(tm,poft_arr,yrange=[min(poft_arr),max(poft_arr)],/current)
;p.ytitle='np*Tp ($10^{-11}$ J/m$^{3}$)'
;p.xtitle='time (s)'

;fit = poly_fit(tm,poft_arr,1,sigma=sig)

;p = plot(tm,fit(0)+fit(1)*tm,'b',/current,/overplot)

w1=window()
dEdt_1 = (shift(poft_arr,-1)-shift(poft_arr,1))/(2*dt*200.)/1e-15
p1 = plot(tm,dEdt_1>0,'4r-D',/current)
p1.sym_filled=1
p1.ytitle='Heating rate density [$10^{-15}$ W/m$^{3}$]'
p1.xtitle='time (s)'
;dEdt = (1e17-7.5e16)/1000
p1.save,'heat_rate.png'

;dEdt = fit(1)
;print,dEdt,sig

;p.title='Heating rate density: '+strmid(strtrim(string(dEdt/1e-15),2),0,4)+'$\pm$'+strmid(strtrim(string(sig(1)/1e-15),2),0,4)+' [$10^{-15}$ W/m$^{3}$]'
;p.save,'heat_rate_run7.png'

;bx = b1(*,ny/2,nz/2,0)
;bz = b1(*,ny/2,nz/2,2)
;ux = up(*,ny/2,nz/2,0)*1e3
;uz = up(*,ny/2,nz/2,2)*1e3
;dens = np(*,ny/2,nz/2)/1e9
;temp_p = tp(*,ny/2,nz/2)

;print,'bx...',b1(*,ny/2,nz/2,0)

;x = x*1e3
;save,filename='KH_profiles.sav',x,bx,bz,ux,uz,dens,temp_p

barr = fltarr(nx)
;for i = 0,nx-1 do begin
pmhd = 0
pkaw = 0
w2 = window()
for i = 1,nfr do begin
   nfrm = i
   c_read_3d_vec_m_32,dir,'c.b1',nfrm,b1
   c_read_3d_m_32,dir,'c.np',nfrm,np
   b1 = b1*mproton/q
   cnt = 0
   pmhd_arr = 0
   pkaw_arr = 0

   
   
   for k = nz/2-15,nz/2+15 do begin
      for j = 1,ny-2 do begin
;k = nz/2
;      j= ny/2
      bx = reform(b1(*,j,k,0))
      bz = reform(b1(*,j,k,2))
      density = reform(np(*,j,k))
      barr = sqrt(bx^2 + bz^2) ;reform(b1(*,k,nz/2-1,2))
      ;barr = bx
;      plot,x,barr
      cnt = cnt+1
      fx = FFT_PowerSpectrum(bx, dx, FREQ=freq)
      fz = FFT_PowerSpectrum(bz, dx, FREQ=freq)
      
      freq = 2*!pi*freq(1:*)
      fx = fx(1:*)
      fz = fz(1:*);*n_elements(freq)
      
      psd = (2*(fz+fx))^(3./2); + (2*fx)^(3./2)
      rho = mp*density/1e9
      wh_mhd = where(freq lt 0.5*2*!pi/cwpi)
      wh_kaw = where(freq ge 0.1*2*!pi/cwpi)
      
      pwr_mhd = (psd(wh_mhd)*freq(wh_mhd))/sqrt(muo^3*rho(wh_mhd))
      pwr_kaw = (freq(wh_kaw)*cwpi)*(psd(wh_kaw)*freq(wh_kaw))/sqrt(muo^3*rho(wh_mhd))
;           pwr_kaw = (psd(wh_kaw)*freq(wh_kaw))/sqrt(muo^3*rho)

; Plot the results
;      w2.erase
;      p2 = plot(freq(wh_mhd), pwr_mhd, /YLOG, /xlog,XTITLE='$k_\perp$',/xsty,/ysty,/current)
;      p2 = plot(freq,fx,/ylog,/xlog,XTITLE='$k_\perp$',/xsty,/ysty,/current)
;      xrange=[min(freq),max(freq)]
;      yrange=[min(pwr_mhd),max(pwr_mhd)]
;      P2.ytitle='heating rate density [W/m$^3$]'
;      p2 = plot([2*!pi/cwpi,2*!pi/cwpi],[min(pwr_mhd),max(pwr_mhd)],':',/overplot,/current)
;      p2 = plot([freq(4),freq(4)],[min(pwr_mhd),max(pwr_mhd)],':',/overplot,/current)
;     fkx = 5e-32*freq^(-5./3.)
;     wh = where(freq lt 0.5*2*!pi/cwpi)
;     p = plot(freq(wh),fkx(wh),/overplot,/current,'2r')
;     fkx = 5e-37*freq^(-8./3.)
;     wh = where(freq gt 0.5*2*!pi/cwpi)
;     p = plot(freq(wh),fkx(wh),/overplot,/current,'2b')
      
;      wait,0.2

;     p2 = plot(freq(wh_kaw),pwr_kaw,/overplot,'2b')
      
      
      wh = where((freq lt 0.25*2*!pi/cwpi) and (freq gt freq(0)))
;      wh = where((freq lt max(freq)/2) and (freq gt f(2)))
;      wh = where(freq lt 2*!pi/cwpi)
;      print,'wh....',n_elements(wh)
;print,'pwr_mhd...',total(pwr_mhd(wh))/n_elements(wh)
;print,'pwr_kaw...',total(pwr_kaw)/n_elements(wh_kaw)
      
      pmhd_arr = pmhd_arr + mean(pwr_mhd(wh));total(pwr_mhd(wh))/n_elements(wh)
      pkaw_arr = pkaw_arr + total(pwr_kaw)/n_elements(wh_kaw)
      
;      p2 = plot(freq(wh),pwr_mhd(wh),/overplot,/current,'2r')
      
   endfor
   endfor
   print,'pmhd...',pmhd_arr/cnt
   print,'pkaw...',pkaw_arr/cnt
   
   pmhd = [pmhd,(pmhd_arr)/cnt]
   pkaw = [pkaw,(pkaw_arr)/cnt]

endfor

p1.window.SetCurrent
p1 = plot(tm,pmhd(4:*)/1e-15,'4b-D',/overplot)
p1.sym_filled=1
p1 = plot(tm,pkaw(4:*)/1e-15,'4g-D',/overplot)
p1.sym_filled=1
p1.ylog=1
p1.yrange=[1e-2,10]

;w4 = window()

;p = plot(freq,f,/current,/xlog,/ylog)

fkx = 5e-20*freq^(-5./3.)

wh = where(freq lt 2*!pi/cwpi)
p = plot(freq(wh),fkx(wh),/overplot,/current,'2r')

;fkx = 4e-11*freq^(-8./3.)

;wh = where(freq gt 2*!pi/cwpi)
;p = plot(freq(wh),fkx(wh),/overplot,/current,'2b')

;fkx = 3e-10*freq^(-7./3.)

;wh = where(freq gt 2*!pi/cwpi)
;p = plot(freq(wh),fkx(wh),/overplot,/current,'2g')

;p.save,'ps_kaw.png'

end
