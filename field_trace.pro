;----------------------------------------------------------------
pro plot_image,img,x,y,sclf,loc,tit
;----------------------------------------------------------------
  
   im = image(img(*,*),x,y,rgb_table=33,layout=[1,1,loc],font_size=12,aspect_ratio=1.0)

   xax = axis('x',axis_range=[min(x),max(x)],location=[0,min(y(*))],thick=2,tickdir=1,target=im,tickfont_size=12)
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
;----------------------------------------------------------------


;----------------------------------------------------------------
PRO read_part,file,nfrm,Ni_max,xp
;----------------------------------------------------------------

;Ni_max=long(0)
;nt=0l
;ntout=0l
frm=0l

file = file+'.dat'
print,' reading...',file
openr,1,file,/f77_unformatted
;readu,1,nt
;readu,1,ntout
;readu,1,Ni_max
;print,nt,ntout,Ni_max


xp=fltarr(Ni_max,3,/nozero)

readu,1,frm
print,'  image #.....',frm
readu,1,xp
frmcnt = 1

while (frmcnt lt nfrm) do begin

   readu,1,frm
   print,'  image #.....',frm
   readu,1,xp
   frmcnt = frmcnt + 1

endwhile

close,1

return
end

;----------------------------------------------------------------



;----------------------------------------------------------------
PRO read_part_scalar,file,nfrm,Ni_max,xp
;----------------------------------------------------------------

;Ni_max=long(0)
;nt=0
;ntout=0
frm=0

file = file+'.dat'
print,' reading...',file
openr,1,file,/f77_unformatted
;readu,1,nt
;readu,1,ntout
;readu,1,Ni_max
;print,nt,ntout,Ni_max

xp=fltarr(Ni_max,/nozero)

readu,1,frm
print,'  image #.....',frm
readu,1,xp
frmcnt = 1

while (frmcnt lt nfrm) do begin

   readu,1,frm
   print,'  image #.....',frm
   readu,1,xp
   frmcnt = frmcnt + 1

endwhile

close,1

return
end
;----------------------------------------------------------------

;-----------------------------------------------------------------------
pro get_edotb,a,b,c
;-----------------------------------------------------------------------
sz = size(a)
nz = sz(3)
a2 = sqrt(a(*,*,*,0)^2 +a(*,*,*,1)^2+a(*,*,*,2)^2)
ahat = a
;ahat(*,*,0)= ahat(*,*,1)
ahat(*,*,nz-1)= ahat(*,*,nz-2)
;for m = 0,2 do begin
;   ahat(*,*,*,m) = a(*,*,*,m);/a2(*,*,*)
;endfor

b2 = sqrt(b(*,*,*,0)^2 +b(*,*,*,1)^2+b(*,*,*,2)^2)
;stop
bhat = b
;bhat(*,*,0)= bhat(*,*,1)
bhat(*,*,nz-1)= bhat(*,*,nz-2)
for m = 0,2 do begin
   bhat(*,*,*,m) = b(*,*,*,m)/b2(*,*,*)
endfor

c(*,*,*) = ahat(*,*,*,0)*bhat(*,*,*,0) + ahat(*,*,*,1)*bhat(*,*,*,1) + ahat(*,*,*,2)*bhat(*,*,*,2)

return
end
;-----------------------------------------------------------------------


;----------------------------------------------------------------------
;main program
;----------------------------------------------------------------------
pro field_trace,nfrm,yslc

@clr_win
;yslc = 150
  
;dir = './run_29/'
;dir = './run_b3_v3_4/'

;dir = './run_va_0.8_beta_3/'
dir = '/data/KH3D/run_test/

nframe=nfrm
read_para,dir
restore,filename=dir+'para.sav'
read_coords,dir,x,y,z
@get_const

valfven = b0_top/sqrt(np_top*1e-9*mproton*4*3.14e-7)/1e3
wpi = sqrt(q*q*(np_top/1e9)/(epo*mproton))
coverwpi = 3e8/wpi/1e3
print,coverwpi
x = x/coverwpi
y = y/(coverwpi)
z = z/coverwpi

xx = [x-max(x),x,x+max(x)]
xx = xx(1:n_elements(xx)-2)

;np = np_top/1e9
;va = b0_bottom/sqrt(muo*np*mproton)/1e3

;c_read_3d_m_32,dir,'c.np',nfrm,np
c_read_3d_m_32,dir,'c.mixed',nfrm,mix
;c_read_3d_m_32,dir,'c.temp_p',nfrm,temp
c_read_3d_vec_m_32,dir,'c.b1',nfrm,b1
c_read_3d_vec_m_32,dir,'c.E',nfrm,Efld
c_read_3d_vec_m_32,dir,'c.up',nfrm,up

mixarr = fltarr(3*nx-2,ny,nz)
mixarr(0:nx-2,*,*) = mix(0:nx-2,*,*)
mixarr(nx-1:2*nx-2,*,*) = mix(*,*,*)
mixarr(2*nx-1:3*nx-3,*,*) = mix(1:nx-1,*,*)

edotb = fltarr(nx,ny,nz)
Efld(*,*,0,*) = Efld(*,*,1,*)
Efld(*,*,nz-1,*) = Efld(*,*,nz-2,*)
b1(*,*,0,*) = b1(*,*,1,*)
b1(*,*,nz-1,*) = b1(*,*,nz-2,*)
Efld0 = valfven*b0_top*1.6e-19/mproton
get_edotb,Efld,b1,edotb
edotb = 1e3*edotb*mproton/1.6e-19 ;E_parallel with units of V/m  ;/Efld0

edotbarr = fltarr(3*nx-2,ny,nz)
edotbarr(0:nx-2,*,*) = edotb(0:nx-2,*,*)
edotbarr(nx-1:2*nx-2,*,*) = edotb(*,*,*)
edotbarr(2*nx-1:3*nx-3,*,*) = edotb(1:nx-1,*,*)

uparr = fltarr(3*nx-2,ny,nz,3)
uparr(0:nx-2,*,*,*) = up(0:nx-2,*,*,*)
uparr(nx-1:2*nx-2,*,*,*) = up(*,*,*,*)
uparr(2*nx-1:3*nx-3,*,*,*) = up(1:nx-1,*,*,*)

plot_image,reform(mix(*,0,*))<1.0,x,z,0.8,1,'mix'
plot_image,reform(mix(*,yslc,*))<1.0,x,z,0.8,1,'mix'

data = FLTARR(3,3*n_elements(x)-2, n_elements(y), n_elements(z))
data(0,0:nx-2,*,*) = b1(0:nx-2,*,*,0)
data(0,nx-1:2*nx-2,*,*) = b1(*,*,*,0)
data(0,2*nx-1:3*nx-3,*,*) = b1(1:nx-1,*,*,0)
data(1,0:nx-2,*,*) = b1(0:nx-2,*,*,1)
data(1,nx-1:2*nx-2,*,*) = b1(*,*,*,1)
data(1,2*nx-1:3*nx-3,*,*) = b1(1:nx-1,*,*,1)
data(2,0:nx-2,*,*) = b1(0:nx-2,*,*,2)
data(2,nx-1:2*nx-2,*,*) = b1(*,*,*,2)
data(2,2*nx-1:3*nx-3,*,*) = b1(1:nx-1,*,*,2)
;data(1,*,*,*) = b1(*,*,*,1)
;data(2,*,*,*) = b1(*,*,*,2)

xstep = 1
zstep = 1
nseeds = 3*LONG((nx*ny)/(xstep*zstep))
seeds = FLTARR(nseeds)
iseed=0L

;k = round(nz/2.)-2
i = round(nx/2.)
for k = 1,nz-1 do begin
   for i = nx-1,2*nx-2 do begin
      print,i,k
      if ((xx(i) gt xx(1)) and (xx(i) lt max(xx)) $
         and (z(k) gt z(1)) and (z(k) lt max(z))) then begin
;         and (y(j) gt y(1)) and (y(j) lt max(y))) then begin
         seeds(iseed) = FLOAT(i)
         seeds(iseed+1) = FLOAT(0)
         seeds(iseed+2) = FLOAT(k)
;         print,seeds(iseed:iseed+3),iseed
         iseed = iseed+3
;         print,iseed,i,k
      endif
   endfor
endfor

maxIterations=3000
stepSize=0.5

PARTICLE_TRACE,data,seeds(0:iseed-1),verts,conn,outnormals, $
               MAX_ITERATIONS=maxIterations, MAX_STEPSIZE=stepSize,  $
               INTEGRATION=1,ANISOTROPY=[dx,dy,delz]/dx, SEED_NORMAL=[1, 0, 0]

epotarr = reform(edotbarr(*,0,*))
mixarr1=reform(mixarr(*,0,*))
uparr1=reform(uparr(*,0,*,*))

i = 0
sz = SIZE(verts, /STRUCTURE)
WHILE (i LT sz.dimensions[1]) DO BEGIN
   nverts = conn[i]
      
   xIndices = reform(verts[0, conn[i+1:i+nverts]])
   yIndices = reform(verts[1, conn[i+1:i+nverts]])
   zIndices = reform(verts[2, conn[i+1:i+nverts]])

   wh = where(abs(yIndices - yslc) eq min(abs((yIndices-yslc))))
   
   i0 = round(xIndices(0))
   j0 = round(yIndices(0))
   k0 = round(zIndices(0))
   
;   i1 = round(xIndices(nverts-1))-(nx-1)
;   j1 = round(yIndices(nverts-1))
;   k1 = round(zIndices(nverts-1))
 
   i1 = round(xIndices(wh(0)))
   j1 = round(yIndices(wh(0)))
   k1 = round(zIndices(wh(0)))
      
;   if (i1 gt nx-1) then i1 = nx-1
   
;   print,'indices...',wh(0),i0,i1,j0,j1,k0,k1,mixarr(i0,j0,k0),mixarr(i1,j1,k1)
   i += nverts + 1
   mixarr1(i0,k0) = mixarr(i1,j1,k1)
   uparr1(i0,k0,0) = uparr(i1,j1,k1,0)
   uparr1(i0,k0,1) = uparr(i1,j1,k1,1)
   uparr1(i0,k0,2) = uparr(i1,j1,k1,2)

   for ii = 1,n_elements(xIndices(0:wh(0)))-1 do begin
      x1 = xx(xIndices(ii))*coverwpi
      x0 = xx(xIndices(ii-1))*coverwpi
      y1 = y(yIndices(ii))*coverwpi
      y0 = y(yIndices(ii-1))*coverwpi
      z1 = z(zIndices(ii))*coverwpi
      z0 = z(zIndices(ii-1))*coverwpi
      dl = sqrt((x1-x0)^2 + (y1-y0)^2 + (z1-z0)^2)
      
;      dl = sqrt(xx(xIndices(ii))^2 + y(yIndices(ii))^2 + z(zIndices(ii))^2)
;      print,'dl..',dl
      epotarr(i0,k0) = epotarr(i0,k0) $
         + edotbarr(round(xIndices(ii)),round(yIndices(ii)),round(zIndices(ii)))*dl*1e3
   endfor
   
ENDWHILE

up1_2 = sqrt(uparr1(*,*,0)^2 + uparr1(*,*,1)^2 + uparr1(*,*,2)^2)
up0_2 = reform(sqrt(uparr(*,0,*,0)^2 + uparr(*,0,*,1)^2 + uparr(*,0,*,2)^2))

uparr0 = reform(uparr(*,0,*,*))

sz1 = size(uparr1)
for i = 0,sz1(1)-1 do begin
   for j = 0,sz1(2)-1 do begin
      uparr1(i,j,0) = uparr1(i,j,0)/up1_2(i,j)
      uparr1(i,j,1) = uparr1(i,j,1)/up1_2(i,j)
      uparr1(i,j,2) = uparr1(i,j,2)/up1_2(i,j)
      uparr0(i,j,0) = uparr0(i,j,0)/up0_2(i,j)
      uparr0(i,j,1) = uparr0(i,j,1)/up0_2(i,j)
      uparr0(i,j,2) = uparr0(i,j,2)/up0_2(i,j)
   endfor
endfor

;uparr0 = reform(uparr(*,0,*,*))/up0_2

up1dotup0 = uparr1(*,*,0)*uparr0(*,*,0) + uparr1(*,*,1)*uparr0(*,*,1) +uparr1(*,*,2)*uparr0(*,*,2)

plot_image,reform(mixarr1(nx-1:2*nx-2,*)<1.0),x,z,0.8,1,'mix'
plot_image,reform(abs(mixarr1(nx-1:2*nx-2,*)-mixarr(nx-1:2*nx-2,0,*)))<1.0,x,z,0.8,1,'mix'
plot_image,epotarr(nx-1:2*nx-2,*),x,z,0.8,1,'Field-aligned potential (V)'
plot_image,up1dotup0(nx-1:2*nx-2,*),x,z,0.8,1,'upxdiff'

return

end
