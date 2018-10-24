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


;----------------------------------------------------------------
pro get_dist,xx,yy,zz,xp,vp,x,y,z,up,temp,va,mrat,B,wh_part
;----------------------------------------------------------------

B = reform(B)
bx = B(0)
by = B(1)
bz = B(2)
B2 = bx^2 + by^2 + bz^2

bxhat = bx/sqrt(B2)
byhat = by/sqrt(B2)
bzhat = bz/sqrt(B2)

vscl=12
dx = x(1)-x(0)

xd = x(floor(xx))
yd = y(floor(yy))
zd = z(floor(zz))

r = sqrt((xp(*,0)-xd)^2 + (xp(*,1)-yd)^2 + (xp(*,2)-zd)^2)

wh = where(r lt 3*dx and mrat lt 1.0)
wh_part = wh
whh = where(r lt 3*dx and mrat lt 1.0)

window,1
!p.multi=[0,2,1]
plot,[-vscl,vscl],[-vscl,vscl],/nodata,/isotropic,/xsty,/ysty,$
     xtitle='vx',ytitle='vz'
plots,vp(wh,0),vp(wh,2),psym=3
plots,vp(whh,0),vp(whh,2),psym=3,color=240

m = bzhat/bxhat
plots,[-vscl,vscl],[-vscl,vscl]*m,linestyle=1



u0 = up(fix(xx),fix(yy),fix(zz),*)/va
vth = sqrt(2*temp(fix(xx),fix(yy),fix(zz))*(1e-19/1.67e-27))/1e3/va

fv = findgen(100,100)

ux = -5+findgen(100)/10
uz = -5+findgen(100)/10
for i = 0,99 do begin
   for k = 0,99 do begin
      fv(i,k) = exp(-((ux(i)-u0(0))^2 + (uz(k)-u0(2))^2)/(vth^2))
   endfor
endfor


contour,fv,ux,uz,/overplot,levels=[1/exp(3),1/exp(2),1/(exp(1)),0.9]


plot,[-vscl,vscl],[-vscl,vscl],/nodata,/isotropic,/xsty,/ysty,$
          xtitle='vx',ytitle='vy'
plots,vp(wh,0),vp(wh,1),psym=3
plots,vp(whh,0),vp(whh,1),psym=3,color=240

m = byhat/bxhat
plots,[-vscl,vscl],[-vscl,vscl]*m,linestyle=1



nvx = 120
dx=0.2
vx_arr = (findgen(nvx)-nvx/2)*dx
vy_arr = (findgen(nvx)-nvx/2)*dx
fv_p =fltarr(nvx,nvx)


for k = 0,n_elements(wh)-1 do begin

   vx = vp(wh(k),0)
   vy = vp(wh(k),1)

   i = round(vx/dx)+nvx/2
   j = round(vy/dx)+nvx/2

   fv_p(i,j) = fv_p(i,j) + 1.0

endfor

wh = where(fv_p ge 1)
fv_p = alog(fv_p)

sz = size(fv_p)
scl = round(1000/sz(1))

img = image(fv_p,vx_arr,vy_arr,rgb_table=39,aspect_ratio=1.0,/current,font_size=20)
xax = axis('x',axis_range=[-max(vx_arr),max(vx_arr)],location=[-max(vx_arr),-max(vx_arr)],thick=2,tickdir=1,target=img,tickfont_size=12)
yax = axis('y',axis_range=[-max(vy_arr),max(vy_arr)],location=[-max(vy_arr),-max(vy_arr)],thick=2,tickdir=1,target=img,tickfont_size=12)
img.scale,scl,scl
img.xtitle='$v_\perp$'
img.ytitle='$v_\parallel$'
img = plot([-vscl,vscl],[-vscl,vscl]*m,'2--',/current,/overplot)
img.xrange=[-4,4]
img.yrange=[-4,4]

img.save,'test.png'


;nvr = 150
;nphi = 36
;fv_p = fltarr(nvr,nphi)
;dvr = 0.05
;dphi=10*!dtor

;vr_arr = findgen(nvr)*dvr
;phi_arr = findgen(nphi)*dphi
;for k = 0,n_elements(whh)-1 do begin

;   if (abs(vr) lt 5) then begin
;      vr = sqrt(vp(whh(k),0)^2 + vp(whh(k),1)^2)
;      phi = atan(vp(whh(k),1),vp(whh(k),0)) + !pi
      
;      i = floor(vr/dvr)
;      j = floor(phi/dphi)
      
;      fv_p(i,j) = fv_p(i,j) + 1.0
;   endif
;endfor

;img = polar_surface(fv_p,vr_arr,phi_arr,/grid)

;;wh = where(img gt 1)
;;img = alog(img(wh))

;im = bytscl(img)
;wh = where(im eq 0)
;im(wh) = 255
;img = image(im,rgb_table=39)
;; xax = axis('x',axis_range=[min(vx),max(vx)],location=[0,0],thick=2,tickdir=1,target=img,tickfont_size=12)
;;   yax = axis('y',axis_range=[min(vy),max(vy)],location=[0,0],thick=2,tickdir=1,target=img,tickfont_size=12)
;img.scale,5,5

contour,fv,ux,uz,/overplot,levels=[1/exp(3),1/exp(2),1/(exp(1)),0.9]



return
end
;----------------------------------------------------------------------

;----------------------------------------------------------------------
;main program
;----------------------------------------------------------------------
pro get_fv,nfrm

;nfrm = v
;dir = '/Volumes/Scratch/hybrid/KH_new/run_3d_29/'
dir = './run_test/'

nframe=nfrm
read_para,dir
restore,filename=dir+'para.sav'
read_coords,dir,x,y,z
@get_const

np = np_top/1e9
va = b0_bottom/sqrt(muo*np*mproton)/1e3

c_read_3d_m_32,dir,'c.np',nfrm,np
c_read_3d_m_32,dir,'c.mixed',nfrm,mix
c_read_3d_m_32,dir,'c.temp_p',nfrm,temp
c_read_3d_vec_m_32,dir,'c.b1',nfrm,b1
c_read_3d_vec_m_32,dir,'c.up',nfrm,up

nfrm = nfrm/4
i=0
read_part,dir+'c.xp_'+strcompress(i,/remove_all),nfrm,Ni_max,xp
read_part,dir+'c.vp_'+strcompress(i,/remove_all),nfrm,Ni_max,vp
read_part_scalar,dir+'c.mrat_'+strcompress(i,/remove_all),nfrm,Ni_max,mrat


for i = 1,9 do begin
   read_part,dir+'c.xp_ '+strcompress(i,/remove_all),nfrm,Ni_max,xpp
   read_part,dir+'c.vp_ '+strcompress(i,/remove_all),nfrm,Ni_max,vpp
   read_part_scalar,dir+'c.mrat_ '+strcompress(i,/remove_all),nfrm,Ni_max,mratt

   xp = [xp,xpp]
   vp = [vp,vpp]
   mrat = [mrat,mratt]
;read_part_scalar,dir+'c.mrat_0',nfrm,Ni_max,mrat

endfor
for i = 10,11 do begin
   read_part,dir+'c.xp_'+strcompress(i,/remove_all),nfrm,Ni_max,xpp
   read_part,dir+'c.vp_'+strcompress(i,/remove_all),nfrm,Ni_max,vpp
   read_part_scalar,dir+'c.mrat_'+strcompress(i,/remove_all),nfrm,Ni_max,mratt
   xp = [xp,xpp]
   vp = [vp,vpp]
   mrat = [mrat,mratt]
;read_part_scalar,dir+'c.mrat_0',nfrm,Ni_max,mrat

endfor


i=0
nfrm = 1
read_part,dir+'c.xp_'+strcompress(i,/remove_all),nfrm,Ni_max,xp1
read_part,dir+'c.vp_'+strcompress(i,/remove_all),nfrm,Ni_max,vp1
read_part_scalar,dir+'c.mrat_'+strcompress(i,/remove_all),nfrm,Ni_max,mrat1
for i = 1,9 do begin
   read_part,dir+'c.xp_ '+strcompress(i,/remove_all),nfrm,Ni_max,xpp
   read_part,dir+'c.vp_ '+strcompress(i,/remove_all),nfrm,Ni_max,vpp
   read_part_scalar,dir+'c.mrat_ '+strcompress(i,/remove_all),nfrm,Ni_max,mratt

   xp1 = [xp1,xpp]
   vp1 = [vp1,vpp]
   mrat1 = [mrat1,mratt]
;read_part_scalar,dir+'c.mrat_0',nfrm,Ni_max,mrat
endfor
for i = 10,11 do begin
   read_part,dir+'c.xp_'+strcompress(i,/remove_all),nfrm,Ni_max,xpp
   read_part,dir+'c.vp_'+strcompress(i,/remove_all),nfrm,Ni_max,vpp
   read_part_scalar,dir+'c.mrat_'+strcompress(i,/remove_all),nfrm,Ni_max,mratt
   xp1 = [xp1,xpp]
   vp1 = [vp1,vpp]
   mrat1 = [mrat1,mratt]
;read_part_scalar,dir+'c.mrat_0',nfrm,Ni_max,mrat
endfor


device,decomposed=0
loadct,39


!p.multi=[0,1,1]
yy = ny/2
zz = nz/2+10
w=window(DIMENSIONS=[700,700])
      
while (yy ge 1) do begin
   window,0
   !p.multi=[0,1,1]
   zz = nz/2
   contour,reform(mix(*,*,zz)),/fill,/isotropic,nlev=50
   cursor,xx,yy,/data
   
   while (zz gt 1) do begin
      
      window,0
      !p.multi=[0,1,1]
      contour,reform(mix(*,yy,*)<1.1),/fill,/isotropic,nlev=50
      cursor,xx,zz,/data
      
      w.erase
      get_dist,xx,yy,zz,xp,vp/va,x,y,z,up,temp,va,mrat,b1(xx,yy,zz,*),wh_part

      window,2
      !p.multi=[0,2,1]
      contour,reform(mix(*,yy,*)<1.1),x,z,/fill,/isotropic,nlev=50
      oplot,xp1(wh_part,0),xp1(wh_part,2),psym=5,color=200

      vscl=9
      plot,[-vscl,vscl],[-vscl,vscl],/nodata,/isotropic,/xsty,/ysty,$
           xtitle='vx',ytitle='vz'
      plots,vp1(wh_part,0)/va,vp1(wh_part,1)/va,psym=3,color=240
      

   endwhile
endwhile

return
end
;----------------------------------------------------------------------
