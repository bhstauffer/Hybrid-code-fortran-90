pro plot_image,img,x,y,sclf,loc,tit

   im = image(img,x,y,/current,rgb_table=33,layout=[2,2,loc],/buffer,font_size=12)

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


pro plot_2x2,nf,yslc


  sclf = 1.3
  xsz = 2000.
  ysz = 1000.
  file = 'ion_cyclotron.mp4'
  width = xsz
  height = ysz
  frames = 180
  fps = 30
  speed = 2
  samplerate = 22050L

  
  device,decomposed=0
  loadct,27
  
  dir = '/Volumes/Scratch/hybrid/KH_new/run_3d_14/'
  
  nframe=nf
  read_para,dir
  restore,filename=dir+'para.sav'
  read_coords,dir,x,y,z
@get_const
  wpi = sqrt(q*q*(np_top/1e9)/(epo*mproton))
  coverwpi = 3e8/wpi/1e3
  print,coverwpi
  x = x/coverwpi
  y = y/coverwpi
  z = z/coverwpi

  minz = 1
  maxz = nz-1 
  z = z - z(nz/2)
  z = z(minz:maxz)
  
  XINTERANIMATE, SET=[xsz,ysz, nframe], /SHOWLOAD 
  w = window(dimensions=[xsz,ysz],/buffer)   
  
  video_file = 'KAW.mp4'
  video = idlffvideowrite(video_file)
  framerate = 7.5
;wdims = w.window.dimensions
  stream = video.addvideostream(xsz, ysz, framerate)
  
  
  openr,1,'./tmp1/c.test_part.dat',/f77_unformatted
  xpart = fltarr(3,3)
  readu,1,xpart
  xp = xpart
  while not(eof(1)) do begin
     readu,1,xpart
     xp = [xp,xpart]
  endwhile
  close,1
  

  
  !p.multi=[0,3,1]
  for nfrm = 1,nframe,1 do begin
     
     c_read_3d_vec_m_32,dir,'c.b1',nfrm,b1
     c_read_3d_vec_m_32,dir,'c.up',nfrm,up
     c_read_3d_m_32,dir,'c.np',nfrm,np
     c_read_3d_m_32,dir,'c.mixed',nfrm,mixed
     c_read_3d_m_32,dir,'c.temp_p',nfrm,tp
     
     print,'Average temp....',total(tp)/n_elements(tp)
     
     w.erase
     
     bx = (reform(b1(*,yslc,minz:maxz,0))*mp/q)/1e-9
     by = (reform(b1(*,yslc,minz:maxz,1))*mp/q)/1e-9
     bz = (reform(b1(*,yslc,minz:maxz,2))*mp/q)/1e-9

     binplane = sqrt(bx^2 + bz^2)
     
     nparr = reform(np(*,yslc,minz:maxz))/1e15
     upx = reform(up(*,yslc,minz:maxz,0))
     upy = reform(up(*,yslc,minz:maxz,1))
     upz = reform(up(*,yslc,minz:maxz,2))
     mixarr = reform(mixed(*,yslc,minz:maxz))
     tparr = reform(tp(*,yslc,minz:maxz))
     tparr = tparr/(total(tparr(*,50))/nx)
     
     
;     plot_image,bx,x,z,sclf,1,'$B_x$ (nT)'
     
;     plot_image,by,x,z,sclf,2,'$B_y$ (nT)'
     
     plot_image,binplane,x,z,sclf,1,'$B_{inplane}$ (nT)'
     
;     plot_image,upx,x,z,sclf,4,'$u_x$ (km/s)'
     
;     plot_image,upy,x,z,sclf,5,'$u_y$ (km/s)'
     
;     plot_image,upz,x,z,sclf,6,'$u_z$ (km/s)'
     
;;   p = plot(xp(*,0),xp(*,2),'w.',/overplot,/current,/buffer,axis_style=0)
     
     plot_image,tparr,x,z,sclf,2,'$T_p/<T_p>$'
     
     
     nparr = nparr*tparr/double(mp*nparr)^(5./3.)
     nparr = nparr/(total(nparr(*,50))/nx)
     
     plot_image,nparr,x,z,sclf,3,'$(p/\rho^\gamma)/<p/\rho^\gamma>$'
     
     wh = where(mixarr gt 1.0)
     if(wh(0) gt -1.0) then begin
        mixarr(wh) = 1.0
     endif
     
     plot_image,mixarr,x,z,sclf,4,'mixing'
     
     img = w.CopyWindow()
     
     print, 'Time:', video.put(stream, w.copywindow())
     xinteranimate, frame = nfrm-1, image = img
     
     
  endfor
  
  video.cleanup
  xinteranimate,/keep_pixmaps
  
  
  npim = image(nparr(*,*),rgb_table=33,/buffer)
  xax = axis('x',axis_range=[0,max(x)],location=[0,0],thick=2,target=npim,tickdir=1)
  yax = axis('y',axis_range=[0,max(z)],location=[0,0],thick=2,target=npim,tickdir=1)
  
  ct = colorbar(target = npim,title='$p/\rho^\gamma$',orientation=1,textpos=1)
  npim.scale,0.9,0.9
  
  npim.xtitle='$x (c/\omega_{pi})$'
  npim.ytitle='$z (c/\omega_{pi})$'
  
  npim.save,'entropy.tif'
  
  
  return
end
