;----------------------------------------------------------------
pro plot_image,img,x,y,sclf,loc,tit,s,nplts,plt_dy
;----------------------------------------------------------------

  pltnum = 1
  im = image(img,x,y,rgb_table=33,layout=[1,nplts,pltnum],$
             font_size=10,$
;             aspect_ratio=1.0,$
             /current,margin=0.05,/buffer)
  im.position=[0.1,1-pltnum*plt_dy, 0.9,1-(pltnum-1)*plt_dy]
 
   xax = axis('x',axis_range=[min(s),max(s)],location=[0,min(s(*))],thick=2,tickdir=1,target=im,$
              tickfont_size=10)
   xax.showtext=0
   yax = axis('y',axis_range=[min(y),max(y)],location=[0,0],thick=2,tickdir=1,target=im,$
              tickfont_size=10)

   im.ytitle='$\alpha$'
;   im.xtitle='$s (c/\omega_{pi})$'

;   ct = colorbar(target = im,title=tit,orientation=1,textpos=1,font_size=12);,$
;                 position=[max(x), min(y),
;                 0.05*(max(x)-min(x)),max(y)],/data)
   im.scale,sclf,sclf
   ct = colorbar(target = im,title=tit,textpos=1,$
                 font_size=10,$
                 orientation=1,$
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
for m = 0,2 do begin
   ahat(*,*,*,m) = a(*,*,*,m);/a2(*,*,*)
endfor

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


;----------------------------------------------------------------
pro get_dist,xx,yy,zz,xp,vp,x,y,z,mrat,B,wfv,f_of_alpha,s,nplts,plt_dy
;----------------------------------------------------------------

wfv.SetCurrent


B = reform(B)
bx = B(0)
by = B(1)
bz = B(2)
B2 = bx^2 + by^2 + bz^2

bxhat = bx/sqrt(B2)
byhat = by/sqrt(B2)
bzhat = bz/sqrt(B2)

vscl=6
dx = x(1)-x(0)

xd = x(floor(xx))
yd = y(floor(yy))
zd = z(floor(zz))

r = sqrt((xp(*,0)-xd)^2 + (xp(*,1)-yd)^2 + (xp(*,2)-zd)^2)

wh = where(r lt 3*dx and mrat eq 1.0)
wh_part = wh
whh = where(r lt 3*dx and mrat lt 1.0)

if (whh(0) lt 0) then return

m = bzhat/bxhat

m = byhat/bxhat

nvx = 140
dx=0.2
vx_arr = (findgen(nvx)-nvx/2)*dx
vy_arr = (findgen(nvx)-nvx/2)*dx
fv_p =fltarr(nvx,nvx)

alpha_arr = fltarr(n_elements(whh))

for k = 0,n_elements(whh)-1 do begin

   vx = vp(whh(k),0)
   vy = vp(whh(k),1)
   vz = vp(whh(k),2)

   v2 = sqrt(vx*vx + vy*vy + vz*vz)
   
   vdotb = vx*bxhat + vy*byhat + vz*bzhat
   alpha = acos(vdotb/v2)
   alpha_arr(k) = alpha
   
   vperp = v2*sin(alpha)
   vpar = v2*cos(alpha)
   
endfor

alpha_arr = alpha_arr*!radeg
h = histogram(alpha_arr, binsize=5, LOCATIONS = xbin)

h1 = fltarr(36)
sz = size(h)
h1(0:sz(1)-1)= h

f_of_alpha = [[f_of_alpha],[h1]]

tvscl,f_of_alpha

wfv.erase

img = transpose(f_of_alpha)

sz = size(img)
plot_image,img,findgen(sz(1)),findgen(36)*5,0.9,1,'counts',s,nplts,plt_dy

return
end
;----------------------------------------------------------------------

;----------------------------------------------------------------------
pro field_trace_fv,nfrm,yslc
;----------------------------------------------------------------------
w = window()
@clr_win
  
dir = '/data/hybrid/run_1/'

nplts = 6
plt_dy = 1./(1.05*nplts)
plt_margin = 0.0

device,decomposed=0
loadct,33

wfv = window(dimensions=[700,900],/buffer)
wfv.erase

f_of_alpha = fltarr(36)

nframe=nfrm
read_para,dir
restore,filename=dir+'para.sav'
read_coords,dir,x,y,z
@get_const

valfven = b0_top/sqrt(np_top*1e-9*mproton*4*3.14e-7)/1e3   ;km/s
wpi = sqrt(q*q*(np_top/1e9)/(epo*mproton))
coverwpi = 3e8/wpi/1e3   ;km
print,coverwpi
x = x/coverwpi
y = y/(coverwpi)
z = z/coverwpi

xx = [x-max(x),x,x+max(x)]
xx = xx(1:n_elements(xx)-2)

c_read_3d_m_32,dir,'c.mixed',nfrm,mix
c_read_3d_vec_m_32,dir,'c.b1',nfrm,b1
c_read_3d_vec_m_32,dir,'c.E',nfrm,Efld   ;km/s^2
c_read_3d_vec_m_32,dir,'c.up',nfrm,up
c_read_3d_vec_m_32,dir,'c.aj',nfrm,aj

nfrm = nfrm/4
i=0
read_part,dir+'c.xp_'+strcompress(i,/remove_all),nfrm,Ni_max,xp
read_part,dir+'c.vp_'+strcompress(i,/remove_all),nfrm,Ni_max,vp
read_part_scalar,dir+'c.mrat_'+strcompress(i,/remove_all),nfrm,Ni_max,mrat
read_part_scalar,dir+'c.beta_p_'+strcompress(i,/remove_all),nfrm,Ni_max,pbeta

for i = 1,9 do begin
   read_part,dir+'c.xp_ '+strcompress(i,/remove_all),nfrm,Ni_max,xpp
   read_part,dir+'c.vp_ '+strcompress(i,/remove_all),nfrm,Ni_max,vpp
   read_part_scalar,dir+'c.mrat_ '+strcompress(i,/remove_all),nfrm,Ni_max,mratt
   read_part_scalar,dir+'c.beta_p_ '+strcompress(i,/remove_all),nfrm,Ni_max,pbta

   xp = [xp,xpp]
   vp = [vp,vpp]
   mrat = [mrat,mratt]
   pbeta = [pbeta,pbta]

endfor
for i = 10,11 do begin
   read_part,dir+'c.xp_'+strcompress(i,/remove_all),nfrm,Ni_max,xpp
   read_part,dir+'c.vp_'+strcompress(i,/remove_all),nfrm,Ni_max,vpp
   read_part_scalar,dir+'c.mrat_'+strcompress(i,/remove_all),nfrm,Ni_max,mratt
   xp = [xp,xpp]
   vp = [vp,vpp]
   mrat = [mrat,mratt]
endfor

xp = xp/coverwpi

mixarr = fltarr(3*nx-2,ny,nz)
mixarr(0:nx-2,*,*) = mix(0:nx-2,*,*)
mixarr(nx-1:2*nx-2,*,*) = mix(*,*,*)
mixarr(2*nx-1:3*nx-3,*,*) = mix(1:nx-1,*,*)

edotb = fltarr(nx,ny,nz)
Efld = 1e3*Efld*mproton/1.6e-19   ;V/m
Efld(*,*,0,*) = Efld(*,*,1,*)
Efld(*,*,nz-1,*) = Efld(*,*,nz-2,*)
b1(*,*,0,*) = b1(*,*,1,*)
b1(*,*,nz-1,*) = b1(*,*,nz-2,*)
Efld0 = 1e3*valfven*b0_top   ;V/m
get_edotb,Efld,b1,edotb
edotb = edotb

edotbarr = fltarr(3*nx-2,ny,nz)
edotbarr(0:nx-2,*,*) = edotb(0:nx-2,*,*)
edotbarr(nx-1:2*nx-2,*,*) = edotb(*,*,*)
edotbarr(2*nx-1:3*nx-3,*,*) = edotb(1:nx-1,*,*)

uparr = fltarr(3*nx-2,ny,nz,3)
uparr(0:nx-2,*,*,*) = up(0:nx-2,*,*,*)
uparr(nx-1:2*nx-2,*,*,*) = up(*,*,*,*)
uparr(2*nx-1:3*nx-3,*,*,*) = up(1:nx-1,*,*,*)

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

xstep = 1
zstep = 1
nseeds = 3*LONG((nx*ny)/(xstep*zstep))
seeds = FLTARR(nseeds)
iseed=0L

;k = round(nz/2.)-2
i = round(nx/2.)
for k = nz/2-5,nz/2,2 do begin
   for i = nx+10,2*nx-10,10 do begin
      print,i,k
      if ((xx(i) gt xx(1)) and (xx(i) lt max(xx)) $
         and (z(k) gt z(1)) and (z(k) lt max(z))) then begin
;         and (y(j) gt y(1)) and (y(j) lt max(y))) then begin
         seeds(iseed) = FLOAT(i)
         seeds(iseed+1) = FLOAT(1)
         seeds(iseed+2) = FLOAT(k)
;         print,seeds(iseed:iseed+3),iseed
         iseed = iseed+3
;         print,iseed,i,k
      endif
   endfor
endfor


maxIterations=3000
stepSize=1.0

PARTICLE_TRACE,data,seeds(0:iseed-1),verts,conn,outnormals, $
               MAX_ITERATIONS=maxIterations, MAX_STEPSIZE=stepSize,  $
               INTEGRATION=1,ANISOTROPY=[dx,dy,delz]/dx, SEED_NORMAL=[1, 0, 0]

epotarr = reform(edotbarr(*,0,*))
mixarr1=reform(mixarr(*,0,*))
uparr1=reform(uparr(*,0,*,*))
mixpfl = 0
eparpfl = 0
epotpfl = 0
bperp_pfl = 0
upxpfl = 0
ajypfl = 0
s = 0
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
   
   i1 = round(xIndices(wh(0)))
   j1 = round(yIndices(wh(0)))
   k1 = round(zIndices(wh(0)))
      
   print,'indices...',wh(0),i0,i1,j0,j1,k0,k1,mixarr(i0,j0,k0),mixarr(i1,j1,k1)
   i += nverts + 1
   mixarr1(i0,k0) = mixarr(i1,j1,k1)
   uparr1(i0,k0,0) = uparr(i1,j1,k1,0)
   uparr1(i0,k0,1) = uparr(i1,j1,k1,1)
   uparr1(i0,k0,2) = uparr(i1,j1,k1,2)

   xind0 = round(xIndices(0)) mod nx
   yind0 = round(yIndices(0))
   zind0 = round(zIndices(0))

   s = 0.0
   
   for ii = 0,n_elements(xIndices(0:wh(0)))-1 do begin
      xind = round(xIndices(ii)) mod nx
      yind = round(yIndices(ii))
      zind = round(zIndices(ii))
;      print,'indices...',xind,(xind mod nx),nx
;      epotarr(i0,k0) = epotarr(i0,k0) $
;                       + edotbarr(round(xIndices(ii)),round(yIndices(ii)),round(zIndices(ii)))
      s = [s,sqrt((x(xind)-x(xind0))^2 + (y(yind)-y(yind0))^2 + (z(zind)-z(zind0))^2)]
      ds=(s(-1)-s(-2))*coverwpi
      print,'ds...',ds
      if (ds gt 0) then begin
         get_dist,xind,yind,zind,xp,vp,x,y,z,mrat,b1(xind,yind,zind,*),wfv,f_of_alpha,s,nplts,plt_dy
      endif
         mixpfl = [mixpfl,mix(xind,yind,zind)]
         eparpfl = [eparpfl,edotbarr(xind,yind,zind)]
         dphi = ds*1e3*edotbarr(xind,yind,zind)
         epotpfl = [epotpfl,epotpfl(-1)+dphi]
         bperp_pfl = [bperp_pfl,sqrt(b1(xind,yind,zind,0)^2 + b1(xind,yind,zind,2)^2)/(q*b0_top/mproton)]
         upxpfl = [upxpfl,up(xind,yind,zind,0)/valfven]
         ajypfl = [ajypfl,aj(xind,yind,zind,1)]

         pltnum=2
         p = plot(s(1:*),mixpfl(1:*),/current,layout=[1,nplts,pltnum],/xsty,/ysty,margin=plt_margin,/buffer)
         p.position=[0.1,1-pltnum*plt_dy, 0.9,1-(pltnum-1)*plt_dy]
         ax = p.axes
         ax[0].showtext=0
         p.ytitle='mixing'

         pltnum=3
         p0 = plot(s(1:*),upxpfl(1:*),'r',/current,layout=[1,nplts,pltnum],/xsty,/ysty,margin=plt_margin,/buffer)
         p0.position=[0.1,1-pltnum*plt_dy, 0.9,1-(pltnum-1)*plt_dy]
         ax0 = p0.axes
         ax0[0].showtext=0
         p0.ytitle='$up_x/v_A$'
         p.Scale,0.9,0.9
         p0.Scale,0.9,0.9

         pltnum=4
         p2 = plot(s(1:*),bperp_pfl(1:*),/current,layout=[1,nplts,pltnum],/xsty,/ysty,margin=plt_margin,/buffer)
         p2.position=[0.1,1-pltnum*plt_dy, 0.9,1-(pltnum-1)*plt_dy]
         ax2 = p2.axes
         ax2[0].showtext=0
         p2.ytitle='$B_\perp/B_0$'
         p2.Scale,0.9,0.9

         pltnum=5
         p1 = plot(s(1:*), eparpfl(1:*),/current,layout=[1,nplts,pltnum],/xsty,/ysty,margin=plt_margin,/buffer)
         p1.position=[0.1,1-pltnum*plt_dy, 0.9,1-(pltnum-1)*plt_dy]
         ax1 = p1.axes
         ax1[0].showtext=0
         p1.ytitle='$E_\parallel$ (V/m)'
                                ;     p1.xtitle='s ($c/\omega_{pi})$'
         p1 = plot([min(s),max(s)],[0,0],':',/current,/overplot,layout=[1,6,5],/xsty,/ysty,/buffer)
         p1.Scale,0.9,0.9
         
         pltnum=6
         p3 = plot(s(1:*),epotpfl(1:*),/current,layout=[1,nplts,pltnum],/xsty,/ysty,margin=0.1,/buffer)
         p3.position=[0.1,1-pltnum*plt_dy, 0.9,1-(pltnum-1)*plt_dy]
         p3.ytitle='$\phi_\parallel$ (V)'
         p3.xtitle='s ($c/\omega_{pi})$'
         p3.Scale,0.9,0.9
         
         img = wfv.CopyWindow()
         help,img
         tv,img,true=1
         
      endfor

wfv.save,'./fv_figs/fv'+strcompress(string(i),/remove_all)+'.png'
f_of_alpha = fltarr(36)
eparpfl = 0
bperp_pfl = 0
mixpfl = 0
upxpfl = 0
epotpfl = 0

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
plot_image,smooth(epotarr(nx-1:2*nx-2,*),2),x,z,0.8,1,'Field-aligned potential'
plot_image,up1dotup0(nx-1:2*nx-2,*),x,z,0.8,1,'upxdiff'

stop


return

end
