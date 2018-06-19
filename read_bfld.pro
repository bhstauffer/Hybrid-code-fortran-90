
@get_const

dir = './run_36/'
read_para,dir
restore,filename=dir+'para.sav'
read_coords,dir,x,y,z

wpi = sqrt(q*q*np_top/1e9/(epo*mproton))
cwpi = 3e8/wpi
cwpi = cwpi/1e3

dx = x(1)-x(0)

nfrm = 17
c_read_3d_vec_m_32,dir,'c.b1',nfrm,b1
b1 = b1*mproton/1.6e-19/1e-9 ;nT

barr = fltarr(nx)
;for i = 0,nx-1 do begin
for k = ny/2-5,ny/2+5 do begin
      bx = reform(b1(*,k,nz/2-1,0))
      bz = reform(b1(*,k,nz/2-1,2))
      barr = barr + sqrt(bx^2 + bz^2);reform(b1(*,k,nz/2-1,2))
   endfor
;endfor

barr = barr/(11)

plot,barr
;stop

f = FFT_PowerSpectrum(barr, dx, FREQ=freq)

freq = 2*!pi*freq(1:*)
;freq = freq(1:*)
f = f(1:*)

; Plot the results
p = plot(freq, f*dx, /YLOG, /xlog,XTITLE='$k_\perp$ (km$^{-1}$)',/xsty,/ysty,font_size=14)
xrange=[min(freq),max(freq)]
yrange=[min(f),max(f)]
p.ytitle='Power (nT$^2$ km)'
p = plot([!pi/cwpi,!pi/cwpi],[min(f*dx),max(f*dx)],':',/overplot,/current)

fkx = 2e-5*freq^(-5./3.)

wh = where(freq lt !pi/cwpi)
p = plot(freq(wh),fkx(wh),/overplot,/current,'2r',name=red           )

fkx = 4e-8*freq^(-8./3.)

wh = where(freq gt !pi/cwpi)
p = plot(freq(wh),fkx(wh),/overplot,/current,'2b',name=blue  )

fkx = 1e-7*freq^(-7./3.)

wh = where(freq gt !pi/cwpi)
;p = plot(freq(wh),fkx(wh),/overplot,/current,'2b',name=blue  )

; = legend(target=[red,blue,green])

p.save,'ps_kaw.png'

end
