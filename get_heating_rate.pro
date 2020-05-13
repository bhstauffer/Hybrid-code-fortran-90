;--------------------------------------------------------------------------
pro plot_image,img,x,y,sclf,loc,tit
;--------------------------------------------------------------------------
  
   im = image(img,x,y,/current,rgb_table=33,layout=[1,1,loc],font_size=12)

   xax = axis('x',axis_range=[min(x),max(x)],location=[0,min(y)],thick=2,$
              tickdir=1,target=im,tickfont_size=12)
   yax = axis('y',axis_range=[min(y),max(y)],location=[0,0],thick=2,$
              tickdir=1,target=im,tickfont_size=12)

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

w = window()
@clr_win
@get_const

dir = './run_va_0.8_beta_3/'
read_para,dir
restore,filename=dir+'para.sav'
read_coords,dir,x,y,z

beta = 3
wpi = sqrt(q*q*np_top/1e9/(epo*mproton))
cwpi = 3e8/wpi
rhoi = sqrt(beta)*cwpi
;cwpi = cwpi/1e3

dx = (x(1)-x(0))*1e3

nfr = 16
poft = 0.0

;w = window()
;for i = 1,nfr do begin
;   nfrm = i
;   c_read_3d_vec_m_32,dir,'c.b1',nfrm,b1
;;   c_read_3d_vec_m_32,dir,'c.up',nfrm,up
;   c_read_3d_m_32,dir,'c.temp_p',nfrm,tp
;   c_read_3d_m_32,dir,'c.np',nfrm,np

;   b1 = b1*mproton/q
   
;   p = np*tp
   
;   parr = reform(p(*,*,*))
;   parr1 = reform(p(*,*,*))

;   sz = size(parr)

;   poft = [poft,total(parr1)/(float(sz(1))*float(sz(3))*float(sz(2)))]
   
;endfor

;poft_arr = poft(2:*)*1.6e-19/1e9

;tm = dt*200.*findgen(n_elements(poft_arr))

;w1=window()
;dEdt_1 = (shift(poft_arr,-1)-shift(poft_arr,1))/(2*dt*200.)/1e-15
;p1 = plot(tm,dEdt_1>0,'4r-D',/current)
;p1.sym_filled=1
;p1.ytitle='Heating rate density [$10^{-15}$ W/m$^{3}$]'
;p1.xtitle='time (s)'
;p1.save,'heat_rate.png'

barr = fltarr(nx)
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
   
   bx = reform(b1(*,1,1,0))
   bz = reform(b1(*,1,1,2))
   barr = sqrt(bx^2 + bz^2) 
   fx = FFT_PowerSpectrum(barr, dx, FREQ=freq)
   parr = fltarr(n_elements(fx))
   
   for k = nz/2-5,nz/2+5 do begin
      for j = ny/2-20,ny/2+20 do begin
      bx = reform(b1(*,j,k,0))
      bz = reform(b1(*,j,k,2))
      density = reform(np(*,j,k))
      rho = mean(mp*density/1e9)

;      tvscl,reform(p(*,ny/2,*))

      barr = sqrt(bx^2 + bz^2)  ;reform(b1(*,k,nz/2-1,2))
      cnt = cnt+1
      fx = FFT_PowerSpectrum(bx, dx, FREQ=freq)
      fz = FFT_PowerSpectrum(bz, dx, FREQ=freq)
      fperp = FFT_PowerSpectrum(barr, dx, FREQ= freq)
      
      freq = 2*!pi*freq(1:*)
      kperp = freq
      fx = fx(1:*)
      fz = fz(1:*)
      fperp = fperp(1:*)
      
;      psd = ((abs(fz)^2+abs(fx)^2))  ; + (2*fx)^(3./2)
      parr = parr + fx + fz
;      wh_mhd = where(freq lt sqrt(beta)*2*!pi/cwpi)
;      wh_kaw = where(freq ge sqrt(beta)*2*!pi/cwpi)
      
;      pwr_mhd = (psd(wh_mhd)*kperp(wh_mhd))/sqrt(muo^3*rho)
;      pwr_kaw = (kperp(wh_kaw)*rhoi)*(psd(wh_kaw)*kperp(wh_kaw))/sqrt(muo^3*rho)
      
;      wh = where((kperp lt 0.25*2*!pi/cwpi) and (freq gt freq(0)))
      
;      pmhd_arr = pmhd_arr + mean(pwr_mhd(wh))
;      pkaw_arr = pkaw_arr + total(pwr_kaw)/n_elements(wh_kaw)
      
   endfor
   endfor
;   print,'pmhd...',pmhd_arr/cnt
;   print,'pkaw...',pkaw_arr/cnt
   
;   pmhd = [pmhd,(pmhd_arr)/cnt]
;   pkaw = [pkaw,(pkaw_arr)/cnt]

endfor

parr = parr/cnt
k = freq/((1/sqrt(beta))*2*!pi/cwpi)
psd = parr/(max(parr));>1e-8

wh = where((k ge k(1)) and (k le 0.3))
fkx_mhd = poly_fit(alog(k(wh)),alog(psd(wh)),1)

p = plot(alog(k(1:*)),alog(psd(1:*)))
p.xtitle='$k_\perp \rho_i$'
p.ytitle='Power'
;p.xrange=[min(k),max(k)]
p1 = plot(alog(k(wh)), fkx_mhd(0) + fkx_mhd(1)*alog(k(wh)),'r',/overplot,name=string(fkx_mhd(1)))

wh = where((k gt 0.3) and (k le 1.0))
fkx_mid = poly_fit(alog(k(wh)),alog(psd(wh)),1)

p3 = plot(alog(k(wh)), fkx_mid(0) + fkx_mid(1)*alog(k(wh)),'b',/overplot,name=string(fkx_mid(1)))

wh = where((k gt 1.0) and (k le 2.0))
fkx_kaw = poly_fit(alog(k(wh)),alog(psd(wh)),1)
p2 = plot(alog(k(wh)), fkx_kaw(0) + fkx_kaw(1)*alog(k(wh)),'g',/overplot,name=string(fkx_kaw(1)))

wh = where((k gt k(0)) and (k le 2.0))
fkx_tot = poly_fit(alog(k(wh)),alog(psd(wh)),1)
p4 = plot(alog(k(wh)), fkx_tot(0) + fkx_tot(1)*alog(k(wh)),'c',/overplot,name=string(fkx_tot(1)))
l = legend(target=[p1,p3,p2,p4])


;p1.window.SetCurrent
;p1 = plot(tm,pmhd(4:*)/1e-15,'4b-D',/overplot)
;p1.sym_filled=1
;p1 = plot(tm,pkaw(4:*)/1e-15,'4g-D',/overplot)
;p1.sym_filled=1
;p1.ylog=1
;p1.yrange=[1e-2,10]

;fkx = 5e-20*freq^(-5./3.)

;wh = where(freq lt 2*!pi/cwpi)
;p = plot(freq(wh),fkx(wh),/overplot,/current,'2r')

end
