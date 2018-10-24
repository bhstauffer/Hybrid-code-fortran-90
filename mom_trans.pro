pro plot_image,img,x,y,sclf,loc,tit

   im = image(img(*,*),x,y,/current,rgb_table=33,layout=[2,1,loc],font_size=14,aspect_ratio=1.0)

   xax = axis('x',axis_range=[min(x),max(x)],location=[0,min(y(*))],thick=2,tickdir=1,target=im,tickfont_size=14)
   yax = axis('y',axis_range=[min(y),max(y)],location=[0,0],thick=2,tickdir=1,target=im,tickfont_size=14)

   im.xtitle='$x (c/\omega_{pi})$'
   im.ytitle='$y (c/\omega_{pi})$'

;   ct = colorbar(target = im,title=tit,orientation=1,textpos=1,font_size=12);,$
;                 position=[max(x), min(y),
;                 0.05*(max(x)-min(x)),max(y)],/data)
   im.scale,sclf,sclf
   ct = colorbar(target = im,title=tit,textpos=1,font_size=14,orientation=1,$
                 position=[1.01,0,1.06,1.0],/relative)

   

end


pro mom_trans,nf,slc,rhouu,bb,TM,TR,x,y,w,dt,omega_p

;1 = xy
;2 = yz
;3 = xz

pln = 1


  sclf = 1.0
  xsz = 2000.
  ysz = 1200.
  file = 'ion_cyclotron.mp4'
  width = xsz
  height = ysz
  frames = 180
  fps = 30
  speed = 2
  samplerate = 22050L
 
device,decomposed=0
loadct,27

;dir = '/Volumes/Scratch/hybrid/KH_new/run_3d_30/'
;dir = '/Volumes/Scratch/hybrid/KH3d/run_3_periodic/'
dir='./run1_3d/'

nframe=nf
read_para,dir
restore,filename=dir+'para.sav'
read_coords,dir,x,y,z
@get_const

omega_p = q*B0_top/mp
print,omega_p,dt,0.08/omega_p
;stop

;XINTERANIMATE, SET=[xsz,ysz, nframe], /SHOWLOAD 
;w = window(dimensions=[xsz,ysz],/buffer)   

;video_file = 'KAW.mp4'
;video = idlffvideowrite(video_file)
;framerate = 7.5
;;wdims = w.window.dimensions
;stream = video.addvideostream(xsz, ysz, framerate)

;openr,1,'./tmp1/c.test_part.dat',/f77_unformatted
;xpart = fltarr(3,3)
;readu,1,xpart
;xp = xpart
;while not(eof(1)) do begin
;   readu,1,xpart
;   xp = [xp,xpart]
;endwhile
;close,1

case pln of
   1: begin
      x = x 
      y = y
      dy = dy
   end
   2: begin
      x = y
      y = z
      dy = delz
   end
   3: begin 
      x = x 
      y = z
      dy = delz
   end
endcase

wpi = sqrt(q*q*(np_top/1e9)/(epo*mproton))
coverwpi = 3e8/wpi/1e3
print,coverwpi
x = x/coverwpi
y = y/(coverwpi)
z = z/coverwpi


;mrestart = '_0'
mrestart = ''

;for nfrm = 1,nframe,1 do begin
nfrm = nf

   c_read_3d_vec_m_32,dir,'c.b1'+mrestart,nfrm,b1
   c_read_3d_vec_m_32,dir,'c.up'+mrestart,nfrm,up
   c_read_3d_m_32,dir,'c.np'+mrestart,nfrm,np
   c_read_3d_m_32,dir,'c.mixed'+mrestart,nfrm,mixed
   c_read_3d_m_32,dir,'c.temp_p'+mrestart,nfrm,tp

   comegapi = (z(1)-z(0))*2

;   bx = (reform(b1(*,*,1,1))*16*mp/q)/1e-9
;   b1(0,0,0,0) = 120.
;   print,max(bx)

   w.erase

   case pln of
      1: begin ;xy
         
         bx = (reform(b1(*,*,slc,0))*mp/q)/1e-9
         by = (reform(b1(*,*,slc,1))*mp/q)/1e-9
         bz = (reform(b1(*,*,slc,2))*mp/q)/1e-9
         
         nparr = reform(np(*,*,slc))/1e15
         upx = reform(up(*,*,slc,0))
         upy = reform(up(*,*,slc,1))
         upz = reform(up(*,*,slc,2))
         mixarr = reform(mixed(*,*,slc))
         tparr = reform(tp(*,*,slc))
         
      end
      2: begin ;yz
         
         bx = (reform(b1(slc,*,*,0))*mp/q)/1e-9
         by = (reform(b1(slc,*,*,1))*mp/q)/1e-9
         bz = (reform(b1(slc,*,*,2))*mp/q)/1e-9
         
         nparr = reform(np(slc,*,*))/1e15
         upx = reform(up(slc,*,*,0))
         upy = reform(up(slc,*,*,1))
         upz = reform(up(slc,*,*,2))
         mixarr = reform(mixed(slc,*,*))
         tparr = reform(tp(slc,*,*))
         
      end
      3: begin ;xz
         
         bx = (reform(b1(*,slc,*,0))*mp/q)/1e-9
         by = (reform(b1(*,slc,*,1))*mp/q)/1e-9
         bz = (reform(b1(*,slc,*,2))*mp/q)/1e-9
         
         nparr = reform(np(*,slc,*))/1e15
         upx = reform(up(*,slc,*,0))
         upy = reform(up(*,slc,*,1))
         upz = reform(up(*,slc,*,2))
         mixarr = reform(mixed(*,slc,*))
         tparr = reform(tp(*,slc,*))
         
      end
   endcase

; calculate stresses on boundary

   rhouu = 0.0
   npave = total(np(*,1,2))/nx
   upave = total(up(*,1,2,0))/nx
   rhouu0 = mp*npave/1e9*(upave*1e3)^2
   bb = 0.0

   TM = bx(*,*)*1e-9*bz(*,*)*1e-9/muo/rhouu0
   TR = -mp*nparr(*,*)*1e6*upx(*,*)*upz(*,*)*1e6/rhouu0
   TM(0,0) = 0.3
   TR(0,0) = 0.3

   for i = 0,nx-1 do begin
      for j = 0,ny-1 do begin
;         rhouu = rhouu-mp*nparr(i,j)*1e6*upx(i,j)*upz(i,j)*1e6/rhouu0
;         bb = bb+bx(i,j)*1e-9*bz(i,j)*1e-9/muo/rhouu0
         rhouu = rhouu-mp*nparr(i,j)*1e6*upx(i,j)*upz(i,j)*1e6
         bb = bb+bx(i,j)*1e-9*bz(i,j)*1e-9/muo

      endfor
   endfor
   print,'rhouu...',rhouu/(nx*ny),bb/(nx*ny),(rhouu+bb)/(nx*ny)



;   plot_image,bx,x,y,sclf,1,'Bx (nT)'

;   plot_image,by,x,y,sclf,2,'By (nT)'

;   plot_image,bz,x,y,sclf,3,'Bz (nT)'

;   plot_image,upx,x,y,sclf,4,'ux (km/s)'

;   plot_image,upy,x,y,sclf,5,'uy (km/s)'      

;   plot_image,upz,x,y,sclf,6,'uz (km/s)'

;   plot_image,tparr,x,y,sclf,7,'Tp (eV)'

;   plot_image,nparr,x,y,sclf,8,'np (cm$^{-3}$)'

;   wh = where(mixarr gt 1.0)
;   if(wh(0) gt -1.0) then begin
;      mixarr(wh) = 1.0;dir = './tmp1/'

;   endif

;   plot_image,mixarr,x,y,sclf,9,'mixing'

;   img = w.CopyWindow()

;   print, 'Time:', video.put(stream, w.copywindow())
;   xinteranimate, frame = nfrm-1, image = img


;endfor

;video.cleanup
;xinteranimate,/keep_pixmaps

return
end


;main program

xsz = 700.
ysz = 800.

;XINTERANIMATE, SET=[xsz,ysz, nframe], /SHOWLOAD 
w = window(dimensions=[xsz,ysz],/buffer)   

rhouu_arr = 0.0
bb_arr = 0.0

;ntot = 109.*99.

ntot = 209*129

for i = 1,15 do begin
   mom_trans,i,55,rhouu,bb,TM,TR,x,y,w,dt,omega_p
   plot_image,TM,x,y,0.8,1,'$T^M_{xz}$'
   plot_image,TR,x,y,0.8,2,'$T^R_{xz}$'

   w.save,'gifs_mom_trans/gif'+strcompress(10+i)+'.gif'

   rhouu_arr = [rhouu_arr,rhouu]
   bb_arr = [bb_arr,bb]
endfor

tm = dt*100*findgen(n_elements(bb_arr))*omega_p

p=plot(tm,bb_arr/ntot/1e-13,'b',name='$B_xB_z \mu_o^{-1}$',/xsty,font_size=18)
p1=plot(tm,rhouu_arr/ntot/1e-13,'r',name='$\rho u_x u_z$',/overplot)
p2 = plot(tm,(bb_arr+rhouu_arr)/ntot/1e-13,/overplot,name='Total')
print,'max mom trans...',(bb_arr+rhouu_arr)/ntot
p2.ytitle='Momentum flux (10$^{-13}$ N/m$^2$)'
p2.xtitle='time ($\Omega_i^{-1}$)'

print,'dt...',dt

l = legend(target=[p,p1,p2],position=[0.33,0.83],/normal,font_size=16)

l.save,'mom_trans_stress.png'

end
