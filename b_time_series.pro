pro b_time_series,nf


  xsz = 1000.
  ysz = 500.
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



dir = './tmp/'

nframe=nf
read_para,dir
restore,filename=dir+'para.sav'
read_coords,dir,x,y,z
@get_const

dt = 2.0*dt

;xsz = 1000.
;ysz = 1000.
;window,xsize=xsz,ysize=ysz
XINTERANIMATE, SET=[xsz,ysz, nframe], /SHOWLOAD 

by = 0.0
tm = 0.0
!p.multi=[0,1,2]
for nfrm = 1,nframe do begin

   c_read_3d_vec_m_32,dir,'c.b1',nfrm,b1

   by = [by,b1(1,1,nz/2+10,1)]
   tm = [tm,dt*nfrm]
   
   comegapi = (z(1)-z(0))*2

;   plot,z,b1(1,1,*,0),yrange=[-3,3]
   plot,z,(b1(1,1,*,0)*23*mp/q)/1e-9,yrange=[-100,100],/ysty
   w = window(dimensions=[xsz,ysz],/buffer)   
   p = plot((z(nz/2)-z)/1800.,(b1(1,1,*,1)*16*mp/q)/1e-9,$
       /current,/buffer,yrange=[-100,100],xrange=[-3,3],$
       xtitle='z ($R_{Io}$)',ytitle='$B_y$ (nT)',font_size=18,thick=2,$
           xstyle=1,xthick=2,ythick=2)
   time = oVid.Put(vidStream, w.CopyWindow())
   w.close

   plot,tm,(by*23*mp/q)/1e-9,xtitle='time (s)',xrange=[0,dt*nframe]

   img = tvrd(0,0,xsz,ysz)

;   tvscl,img
   xinteranimate, frame = nfrm-1, image = img

;   plot,tm,by,xtitle='time (s)'

endfor

oVid.Cleanup


xinteranimate,/keep_pixmaps


return
end
