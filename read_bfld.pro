
@get_const

dir = './run_34/'
read_para,dir
restore,filename=dir+'para.sav'
read_coords,dir,x,y,z

wpi = sqrt(q*q*np_top/1e9/(epo*mproton))
cwpi = 3e8/wpi
cwpi = cwpi/1e3

dx = x(1)-x(0)

nfrm = 9
c_read_3d_vec_m_32,dir,'c.b1',nfrm,b1


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
