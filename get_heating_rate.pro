
pro plot_image,img,x,y,sclf,loc,tit

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

@get_const

dir = './run_va_0.8_159/'
read_para,dir
restore,filename=dir+'para.sav'
read_coords,dir,x,y,z


wpi = sqrt(q*q*np_top/1e9/(epo*mproton))
cwpi = 3e8/wpi
cwpi = cwpi/1e3

dx = x(1)-x(0)

nfr = 16
poft = 0.0
for i = 1,nfr do begin
   nfrm = i
   c_read_3d_vec_m_32,dir,'c.b1',nfrm,b1
   c_read_3d_vec_m_32,dir,'c.up',nfrm,up
   c_read_3d_m_32,dir,'c.temp_p',nfrm,tp
   c_read_3d_m_32,dir,'c.np',nfrm,np

   b1 = b1*mproton/q
   
   p = np*tp
   
;   parr = reform(p(*,ny/2,*))
;   parr1 = reform(p(*,ny/2,*))

   parr = reform(p(*,*,nz/2-60:nz/2+60))
   parr1 = reform(p(*,*,nz/2-60:nz/2+60))

   
;   plot_image,reform(p(*,ny/2,*)),x,z,0.8,1,'p'
   
   sz = size(parr)
;   print,total(parr1)/(float(sz(1))*11.)
   poft = [poft,total(parr1)/(float(sz(1))*float(sz(3))*float(sz(2)))]
   
endfor

poft_arr = poft(4:*)*1.6e-19/1e9/1e-11

tm = dt*200.*findgen(n_elements(poft_arr))

p=plot(tm,poft_arr,yrange=[min(poft_arr),max(poft_arr)])
p.ytitle='np*Tp ($10^{-11}$ J/m$^{3}$)'
p.xtitle='time (s)'

fit = poly_fit(tm,poft_arr,1,sigma=sig)

p = plot(tm,fit(0)+fit(1)*tm,'b',/current,/overplot)

dEdt_1 = (poft_arr-shift(poft_arr,1))/(dt*200.)*1e-11/1e-15
p1 = plot(tm,smooth(dEdt_1>0,2),'4rD')
p1.sym_filled=1
p1.ytitle='Heating rate density [$10^{-15}$ W/m$^{3}$]'
p1.xtitle='time (s)'
;dEdt = (1e17-7.5e16)/1000
p1.save,'heat_rate.png'

dEdt = fit(1)
print,dEdt,sig

p.title='Heating rate density: '+strmid(strtrim(string(dEdt*1e-11/1e-15),2),0,4)+'$\pm$'+strmid(strtrim(string(sig(1)*1e-11/1e-15),2),0,4)+' [$10^{-15}$ W/m$^{3}$]'
p.save,'heat_rate_run7.png'

bx = b1(*,ny/2,nz/2,0)
bz = b1(*,ny/2,nz/2,2)
ux = up(*,ny/2,nz/2,0)*1e3
uz = up(*,ny/2,nz/2,2)*1e3
dens = np(*,ny/2,nz/2)/1e9
temp_p = tp(*,ny/2,nz/2)

;print,'bx...',b1(*,ny/2,nz/2,0)

x = x*1e3
;save,filename='KH_profiles.sav',x,bx,bz,ux,uz,dens,temp_p

stop

barr = fltarr(nx)
;for i = 0,nx-1 do begin
for k = ny/2-5,ny/2+5 do begin
      bx = reform(b1(*,k,nz/2-1,0))
      bz = reform(b1(*,k,nz/2-1,2))
      barr = barr + sqrt(bx^2 + bz^2);reform(b1(*,k,nz/2-1,2))
   endfor
;endfor

barr = barr/(11)


f = FFT_PowerSpectrum(barr, dx, FREQ=freq)

freq = 2*!pi*freq(1:*)
f = f(1:*)

; Plot the results
p = plot(freq, f, /YLOG, /xlog,XTITLE='$k_\perp$',/xsty,/ysty)
xrange=[min(freq),max(freq)]
yrange=[min(f),max(f)]
p.ytitle='Power'
p = plot([2*!pi/cwpi,2*!pi/cwpi],[min(f),max(f)],':',/overplot,/current)

fkx = 5e-9*freq^(-5./3.)

wh = where(freq lt 2*!pi/cwpi)
p = plot(freq(wh),fkx(wh),/overplot,/current,'2r')

fkx = 4e-11*freq^(-8./3.)

wh = where(freq gt 2*!pi/cwpi)
p = plot(freq(wh),fkx(wh),/overplot,/current,'2b')

fkx = 3e-10*freq^(-7./3.)

wh = where(freq gt 2*!pi/cwpi)
p = plot(freq(wh),fkx(wh),/overplot,/current,'2g')

p.save,'ps_kaw.png'

end
