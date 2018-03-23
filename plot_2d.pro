pro plot_2d,nf


  xsz = 1800.
  ysz = 400.
  file = 'ion_cyclotron.mp4'
  width = xsz
  height = ysz
  frames = 180
  fps = 30
  speed = 2
  samplerate = 22050L
 
  ; Create object and initialize video/audio streams
  oVid = IDLffVideoWrite(file)
  vidStream = oVid.AddVideoStream(width, height, fps)
 ; audStream = oVid.AddAudioStream(samplerate)

device,decomposed=0
loadct,27

dir = './tmp1/'

nframe=nf
read_para,dir
restore,filename=dir+'para.sav'
read_coords,dir,x,y,z
@get_const

;xarr = fltarr(nx,nz)
;zarr = fltarr(nx,nz)

;for i = 0,nx-1 do begin
;   zarr(i,*) = (z-z(nz/2))/1800.
;endfor

;for k = 0,nz-1 do begin
;   xarr(*,k) = (x-x(nx/2))/1800.
;endfor


;window,xsize=xsz,ysize=ysz
XINTERANIMATE, SET=[xsz,ysz, nframe], /SHOWLOAD 
w = window(dimensions=[xsz,ysz],/buffer)   

!p.multi=[0,3,1]
for nfrm = 1,nframe,1 do begin

   c_read_3d_vec_m_32,dir,'c.b1',nfrm,b1
   c_read_3d_vec_m_32,dir,'c.up',nfrm,up
   c_read_3d_m_32,dir,'c.np',nfrm,np
   c_read_3d_m_32,dir,'c.temp_p',nfrm,tp


   comegapi = (z(1)-z(0))*2

;   plot,z,b1(nx/2+1,1,*,0),yrange=[-1,1]
;   bx = (reform(b1(*,1,*,1))*16*mp/q)/1e-9
   bx = (reform(b1(*,*,1,2))*16*mp/q)/1e-9
   b1(0,0,0,0) = 120.
   print,max(bx)
;   contour,bytscl(bx),x,z,/ysty,/fill,$
;           /isotropic,/xsty,nlev=10
;   for i = 0,nx-1,6 do begin
;      plots,[x(i),x(i)],[!y.crange(0),!y.crange(1)],linestyle=0,/data
;   endfor
;   for k = 0,nz-1,6 do begin
;      plots,[!x.crange(0),!x.crange(1)],[z(k),z(k)],linestyle=0,/data
;   endfor
;   nparr = median(reform(np(*,1,*))/1e15,3)
;   contour,bytscl(nparr),x,z,/ysty,/fill,$
;           /isotropic,/xsty,nlev=10

;   uparr = reform(up(*,1,*,0))/1e15
;   contour,bytscl(uparr),x,z,/ysty,/fill,$
;           /isotropic,/xsty,nlev=10


   ;plot,z,(b1(nx/2,1,*,1)*16*mp/q)/1e-9

   w.erase
;   arr = reform(np(*,1,*))/1e18

   arr = (reform(b1(*,1,*,1))*mp/q)/1e-9
;   arr = (reform(b1(*,1,*,1))*mp/q)/1e-9
;   bmx = 40
;   arr(0,0) = bmx
;   arr(0,1) = -bmx
;   wh = where(arr gt bmx)
;   arr(wh) = bmx
;   wh = where(arr lt -bmx)
;   arr(wh) = -bmx

;   nparr = reform(np(*,1,*))/1e15
;   uparr = reform(up(*,1,*,0))
   nparr = reform(np(*,1,*))/1e15
   uparr = reform(up(*,1,*,0))
   tparr = reform(tp(*,1,*))


   im = image(arr,/current,rgb_table=33,layout=[3,1,1],/buffer)
;   im = image(temparr,xarr,zarr,/current,rgb_table=33,layout=[3,1,1],/buffer)
   xax = axis('x',axis_range=[0,max(x)],location=[0,0],thick=2)
   yax = axis('y',axis_range=[0,max(z)],location=[0,0],thick=2)

   ct = colorbar(target = im,title='$B_x$ (nT)',orientation=1,textpos=1)

   im1 = image(tparr,/current,rgb_table=33,layout=[3,1,2],/buffer)
   xax = axis('x',axis_range=[0,max(x)],location=[0,0],thick=2,target=im1)
   yax = axis('y',axis_range=[0,max(z)],location=[0,0],thick=2,target=im1)

;   ct = colorbar(target = im1,title='$n_p$ (cm$^{-3}$)',orientation=1,textpos=1)
   ct = colorbar(target = im1,title='$T_p$ (eV)',orientation=1,textpos=1)

   im2 = image(uparr,/current,rgb_table=33,layout=[3,1,3],/buffer)
   xax = axis('x',axis_range=[0,max(x)],location=[0,0],thick=2,target=im2)
   yax = axis('y',axis_range=[0,max(z)],location=[0,0],thick=2,target=im2)

   ct = colorbar(target = im2,title='$u_p$ (km/s)',orientation=1,textpos=1)

;   im.scale,1.2,1.2

;   xdot = reform(up(*,1,*,0))
;   zdot = reform(up(*,1,*,2))

;   p = contour(arr,$
;               xarr,zarr,layout=[3,1,1],$
;               /current,/fill,rgb_table=33,n_levels=10,$
;               xtitle='x ($R_{Io}$)',ytitle='z ($R_{Io}$)',font_size=16,$
;               xstyle=1,xthick=2,ystyle=1,ythick=2,aspect_ratio=1.0,/buffer)

;   ct = colorbar(target = p,title='$B_x$ (nT)',orientation=1)

;   p1 = contour(arr,$
;               xarr,zarr,layout=[3,1,2],$
;               /current,/fill,rgb_table=33,n_levels=10,$
;               xtitle='x ($R_{Io}$)',ytitle='z ($R_{Io}$)',font_size=16,$
;               xstyle=1,xthick=2,ystyle=1,ythick=2,aspect_ratio=1.0,/buffer)

;   p2 = contour(arr,$
;               xarr,zarr,layout=[3,1,3],$
;               /current,/fill,rgb_table=33,n_levels=10,$
;               xtitle='x ($R_{Io}$)',ytitle='z ($R_{Io}$)',font_size=16,$
;               xstyle=1,xthick=2,ystyle=1,ythick=2,aspect_ratio=1.0,/buffer)

;   v = streamline(xdot,zdot,xarr,zarr,/current,overplot=1,$
;                 x_streamparticles=11,y_streamparticles=11)

;   time = oVid.Put(vidStream, w.CopyWindow())
;   w.close

;   plot,tm,(by*16*mp/q)/1e-9,xtitle='time (s)',xrange=[0,dt*nframe]

;   img = tvrd(0,0,xsz,ysz,true=1)
   img = w.CopyWindow()

;   tvscl,img
   xinteranimate, frame = nfrm-1, image = img

;   plot,tm,by,xtitle='time (s)'

endfor

oVid.Cleanup

xinteranimate,/keep_pixmaps

return
end
