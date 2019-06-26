;--------------------------------------------------------------------------------------------
pro plot_image,img,x,y,sclf,loc,tit
;--------------------------------------------------------------------------------------------
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
;--------------------------------------------------------------------------------------------

;--------------------------------------------------------------------------------------------
pro mom_trans,dir,nf,slc,rhouu,bb,TM,TR,x,y,w,dt,omega_p,dudz,rho0,dV0,coverwpi
;--------------------------------------------------------------------------------------------
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
;dir='./run_b1_v3_4/'

nframe=nf
read_para,dir
restore,filename=dir+'para.sav'
read_coords,dir,x,y,z
@get_const

omega_p = q*B0_top/mp
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

rho0 = mproton*(np_top/1e9)

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

   dV0 = total(up(*,1,nz-5,0)-up(*,1,5,0))/nx
   print,dV0
   
; calculate stresses on boundary

   rhouu = 0.0
   npave = total(np(*,1,2))/nx
   upave = total(up(*,1,2,0))/nx
   rhouu0 = mp*npave/1e9*(upave*1e3)^2
   bb = 0.0
   dudz = 0.0
   
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
         dudz = dudz + (up(i,j,slc+1,0) - up(i,j,slc-1,0))/(3*dx)
      endfor
   endfor

   rhouu = rhouu/(nx*ny)
   bb = bb/(nx*ny)
   dudz = dudz/(nx*ny)

return
end

;--------------------------------------------------------------------------------------------
;main program
;--------------------------------------------------------------------------------------------

files=['./run_b0.5_v3_4/','./run_b1_v3_4/','./run_b2_v3_4/','./run_b3_v3_4/','./run_b4_v3_4/','./run_b5_v3_4/','./run_b6_v3_4/','./run_b7_v3_4/']


maxnu = 0.0
maxnu_bb = 0.0
maxnu_rhouu = 0.0

for nf = 0,n_elements(files)-1 do begin

   xsz = 700.
   ysz = 800.
   
;XINTERANIMATE, SET=[xsz,ysz, nframe], /SHOWLOAD 
   w = window(dimensions=[xsz,ysz],/buffer)   
   
   rhouu_arr = 0.0
   bb_arr = 0.0
   dudz_arr = 0.0
   
;ntot = 109.*99.
   
   ntot = 129.*199.
   
   for i = 1,8 do begin
      mom_trans,files(nf),i,55,rhouu,bb,TM,TR,x,y,w,dt,omega_p,dudz,rho0,dV0,coverwpi
      plot_image,TM,x,y,0.8,1,'$T^M_{xz}$'
      plot_image,TR,x,y,0.8,2,'$T^R_{xz}$'
      
      w.save,'gifs_mom_trans/gif'+strcompress(10+i)+'.gif'
      
      rhouu_arr = [rhouu_arr,rhouu]
      bb_arr = [bb_arr,bb]
      dudz_arr = [dudz_arr, dudz]
   endfor
   
   tm = dt*100*findgen(n_elements(bb_arr))*omega_p
   
;dudz_avg = dudz_arr
   
;print,'dudz_avg...',dudz_avg
;dudz_avg(0) = dudz_avg(1)
;stop
   
   nu_anom = (1/rho0)*(bb_arr+rhouu_arr)/dudz_arr
   nu_anom_bb = (1/rho0)*(bb_arr)/dudz_arr
   nu_anom_rhouu = (1/rho0)*(rhouu_arr)/dudz_arr
   nu_anom = nu_anom/(2*1.5*coverwpi*1e3*dV0*1e3)
   nu_anom_bb = nu_anom_bb/(2*1.5*coverwpi*1e3*dV0*1e3)
   nu_anom_rhouu = nu_anom_rhouu/(2*1.5*coverwpi*1e3*dV0*1e3)

   
   p=plot(tm,bb_arr/1e-13,'b',name='$B_xB_z \mu_o^{-1}$',/xsty,font_size=18)
   p1=plot(tm,rhouu_arr/1e-13,'r',name='$\rho u_x u_z$',/overplot)
   p2 = plot(tm,(bb_arr+rhouu_arr)/1e-13,/overplot,name='Total')
   print,'max mom trans...',(bb_arr+rhouu_arr)
   p2.ytitle='Momentum flux (10$^{-13}$ N/m$^2$)'
   p2.xtitle='time ($\Omega_i^{-1}$)'

   maxnu = [maxnu,max(nu_anom)]
   maxnu_bb = [maxnu_bb,max(nu_anom_bb)]
   maxnu_rhouu = [maxnu_rhouu,max(nu_anom_rhouu)]
   
;print,'dt...',dt
   
   l = legend(target=[p,p1,p2],position=[0.33,0.83],/normal,font_size=16)
   
   pnu = plot(tm,nu_anom,/xsty,font_size=18)
   pnu = plot(tm,nu_anom_bb,'b',/overplot)
   pnu = plot(tm,nu_anom_rhouu,'r',/overplot)
   pnu.xtitle='time ($\Omega_i^{-1}$)'
   pnu.ytitle='$\nu_{anom}$ ($2L_0v_0$)'
   
   l.save,'mom_trans_stress.png'

endfor

xarr = [0.5,1,2,3,4,5,6,7]
maxnu = maxnu(1:*)
maxnu_bb = maxnu_bb(1:*)
maxnu_rhouu = maxnu_rhouu(1:*)

b1 = barplot(xarr,maxnu, bottom_color="white",fill_color='black', nbars=3,index=0,xtickdir=1,ytickdir=1,xminor=0,$
           name='total',width=1.2)
b2 = barplot(xarr,maxnu_bb, bottom_color="white",fill_color='blue', nbars=3,index=1,xtickdir=1,ytickdir=1,xminor=0,$
            name='Maxwell',/overplot,width=1.2)
b3 = barplot(xarr,maxnu_rhouu, bottom_color="white",fill_color='red', nbars=3,index=2,xtickdir=1,ytickdir=1,xminor=0,$
            name='Reynolds',/overplot,width=1.2)
b1.xtitle='plasma $\beta$'
b1.ytitle='max($\nu_{anom}$) $(2 L_0 v_0)$'
l = legend(target=[b1,b2,b3],position=[0.33,0.83],/normal,font_size=16)

xarr = 0.8/sqrt(1 + ((5./3.)/2)*xarr)
b1 = barplot(xarr,maxnu, bottom_color="white",fill_color='black', nbars=3,index=0,xtickdir=1,ytickdir=1,xminor=0,$
           name='total',width=0.3)
b2 = barplot(xarr,maxnu_bb, bottom_color="white",fill_color='blue', nbars=3,index=1,xtickdir=1,ytickdir=1,xminor=0,$
            name='Maxwell',/overplot,width=0.3)
b3 = barplot(xarr,maxnu_rhouu, bottom_color="white",fill_color='red', nbars=3,index=2,xtickdir=1,ytickdir=1,xminor=0,$
            name='Reynolds',/overplot,width=0.3)
b1.xtitle='$M_f = 0.8/(1-\gamma \beta/2)^{1/2}$'
b1.ytitle='max($\nu_{anom}$) $(2 L_0 v_0)$'
b1.xrange=[0.3,0.7]
l = legend(target=[b1,b2,b3],position=[0.33,0.83],/normal,font_size=16)


end
