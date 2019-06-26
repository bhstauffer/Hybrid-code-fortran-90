 
;------------------------------------------------------------------------
pro read_data,dir,f,s
;------------------------------------------------------------------------

close,1

read_para,dir
restore,filename=dir+'para.sav'
read_coords,dir,x,y,z
@get_const

va = b0_bottom/sqrt(muo*mproton*np_bottom/1e9)

vo = 0.8*va/1e3
Lo = 1.5*dx

Am = 1.5*(nx-1)*(ny-1)

nfrm = 16

lnvz = fltarr(nfrm)
mix_arr = fltarr(nfrm)

nout=200
;dt = 0.5
tm = dt*nout+findgen(nfrm)*nout*dt

t1 = 280;min(tm)
t2 = max(tm)

t1 = 280;min(tm)
t2 = max(tm)

for i = 1,nfrm do begin 
   c_read_3d_vec_m_32,dir,'c.up',i,up
   c_read_3d_m_32,dir,'c.mixed',i,mixed

   wh = where((mixed le 0.75) and (mixed ge 0.25))
   mix_arr(i-1) = dx*1e3*n_elements(wh)/Am

   lnvz(i-1) = alog(max(abs(smooth(up(*,1,*,2),2))))
endfor

!x.title='time (s)'
!p.multi=[0,2,1]

p = plot(tm,mix_arr^2,xrange=[t1,t2],yrange=[1,max(mix_arr^2)],/xsty,/ysty,$
        layout=[2,1,1])
p.xtitle='time (s)'
p.ytitle='Mixing/(Lo*nx*ny)'

wh = where((tm ge t1) and (tm le t2))
fit = poly_fit(tm(wh),mix_arr(wh)^2,1,sigma=sigma)

f = fit(1)
s = sigma(1)

;f = f*Lo/vo

p = plot(tm,fit(0)+fit(1)*tm,'b',layout=[2,1,2],/current,/overplot)


p = plot(tm,lnvz,yrange=[2.0,5.0],xrange=[t1,t2],/ysty,/xsty,$
        layout=[2,1,2],/current)
p.xtitle='time (s)'
p.ytitle='ln(v!dmax!n)'


return
end
;-------------------------------------------------------------------------

;dir = '/Volumes/Scratch/hybrid/KH_new/run_3d_19/'
;beta runs
files=['./run1_3d/','./run2_3d/','./run3_3d/']

;heavy fraction
files=['./run1_3d/','./run4_3d/','./run5_3d/']

;mach number
files=['./run1_3d/','./run6_3d/','./run7_3d/']


;files=['./run1_3d/','./run2_3d/','./run3_3d/','./run4_3d/','./run5_3d/','./run6_3d/','./run7_3d/']

files = ['./run_va_0.8_121/','./run_va_0.8_159/']

dir = files(0)
read_data,dir,f,s

for i = 1,n_elements(files)-1 do begin

   print,'file...',files(i)
   dir = files(i)
   read_data,dir,ff,ss
   
   f = [f,ff]
   s = [s,ss]

endfor

xarr = indgen(n_elements(files))+1

xarr = [121,159]

;b2 = barplot(xarr,(f-s)/1e9, bottom_color="white", fill_color='yellow', index=0,name='beta 1.0',xtick;dir=1,ytickdir=1)
b = barplot(xarr,f/1e9, bottom_color="white",fill_color='red', index=0,name='beta 1.0',xtickdir=1,ytickdir=1,xminor=0)
;b3 = barplot(xarr,(f+s)/1e9, bottom_color="white", fill_color='yellow', index=0,name='beta 1.0',xtick;dir=1,ytickdir=1,/overplot,bottom_values=f/1e9)
b = errorplot(xarr,f/1e9,s/1e9,/overplot,linestyle=6,thick=4)

;b.title='$\beta$ = 1.0'
;b.xtitle='plasma $\beta$'
;b.xtitle='Heavy ion fraction'
b.xtitle='$\Delta v$ (M_A)'
b.ytitle='$D (10^9 m^2/s)$'
b.yrange=[0,max(f+2*s)/1e9]

;b1 = barplot(xarr,f1/1e9, nbars = 2, fill_color='blue',index=1,/overplot,name='beta 0.5')

;l = legend(target=[b,b1])

b.save,"diffusion_coeff_mach.png"

stop
end
