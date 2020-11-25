pro get_initial,z,v,d,Tel0,Npart,vth

  vth = sqrt(3*Tel0*1.6e-19/9.1e-31)
  z = randomu(seed,Npart)*max(d)
  Rt = randomu(seed,Npart)
  Rs = randomu(seed,Npart)
  vs = vth*sqrt(-2*alog(Rs))*cos(2*!pi*Rt)
  Rt = randomu(seed,Npart)
  Rs = randomu(seed,Npart)
  a = vs*sqrt(1-Rs^2)
  vx = a*cos(2*!pi*Rt) 
  vy = a*sin(2*!pi*Rt)
  v(*,0) = sqrt(vx^2 + vy^2) ;vperp
  v(*,1) = vs*Rs             ;vparallel

end

function vdot,z,dd,Epar,Btot
  q_over_m = -1.6e-19/9.1e-31
  dz = dd(1)-dd(0)
  i = floor(z/dz)
  ip = ceil(z/dz)
  wh = where(ip gt n_elements(dd)-1)
  ip(wh) = 0
  wh = where(i lt 0)
  i(wh) = n_elements(dd)-1

  ;gradB = (Btot(ip)-Btot(i))/dz
  return, vdot = q_over_m*(Epar(i) + (Epar(ip)-Epar(i))*(z-dd(i))/(dd(ip)-dd(i)))
end

pro move_part,z,v,dd,Epar,Btot

  dz = dd(1)-dd(0)
  dt = 0.5*dz/abs(max(v))
  vh = fltarr(n_elements(z),3)
  vh(*,1) = v(*,1) + 0.5*vdot(z,dd,Epar,Btot)*dt
  zh = z + 0.5*v(*,1)*dt

  v(*,1) = v(*,1) + vdot(zh,dd,Epar)*dt
  z = z + vh(*,1)*dt
  
  wh = where(z ge max(dd))
  if (wh(0) gt -1) then begin
     z(wh) = z(wh) - max(dd)
  endif
  wh = where(z lt min(dd))
  if (wh(0) gt -1) then begin
     z(wh) = max(dd) - abs(min(dd) - z(wh))
  endif
     
end
  
dir = './prl_nu_001_noscale/'
B0 = 5e-9

zarr = 0
varr = 0

;for n = 0,n_elements(files)-1 do begin
close,1
openr,1,dir+'files.txt'
file=' '
cnt=0
while not(eof(1)) do begin
   readf,1,file
   print,file
   
   restore,dir+file

;epar in V/m
   btot = btot*b0               ;T
   d = d*coverwpi*1e3           ;m
   
   dd = findgen(n_elements(d))*max(d)/n_elements(d)
   dz = dd(1)-dd(0)
   
   dd = [dd,max(dd)+dz+dd]
   epar = [epar,reverse(epar)]
   btot = [btot,reverse(btot)]   
   
   
;plot,dd,epar
;stop
;oplot,dd,epar,linestyle=1
   
   Tel0 = 30.0                  ;eV
   Npart = 50000
   z = fltarr(Npart)
   v = fltarr(Npart,2)
   
   get_initial,z,v,dd,Tel0,Npart,vth
   
   h0 =  histogram(v(*,1),nbins=20,locations=binvals0)
   
   for i = 0,500 do begin
      move_part,z,v,dd,Epar,Btot
;      plot,v(*,1),v(*,0),psym=1,/isotropic,xtitle='v_parallel',ytitle='v_perp'
;      plot,z,v(*,1),psym=3
     ; wait,0.01
   endfor
   p = plot(z/1e3/1e5,v(*,1)/1e3/1e4,'.',sym_size=0.1,font_size=14)
   p.xtitle='$y (10^5 km)$'
   p.ytitle='$v_\parallel (10^4 km/s)$'
   p.xtickdir=1
   p.ytickdir=1
   p.save,"test_part.png",resolution=300
   zarr = [zarr,z]
   varr = [varr,v(*,1)]
   cnt= cnt+1
endwhile

varr = varr[1:*]
zarr = zarr[1:*]

h = histogram(varr,nbins=20,locations=binvals)
b = barplot(binvals0/1e3/1e4,h0,ylog=1,index=0,nbars=2,fill_color='red',font_size=14)
b = plot(binvals0/1e3/1e4,max(h0)*exp(-binvals0^2/(vth)^2),/overplot)
b.xrange=[-2.5,2.5]
b.xtitle='$v_\parallel (10^4 km/s)$'
b.ytitle='counts'

b = barplot(binvals/1e3/1e4,h/cnt,ylog=1,index=1,nbars=2,color='blue',/overplot)
b = plot(binvals/1e3/1e4,0.5*max(h/cnt)*exp(-binvals^2/(2*vth)^2),/overplot)
;b.xtitle='$v_\parallel (10^4 km/s)$'
;b.ytitle='counts'

end
