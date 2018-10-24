 
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

nfrm = 12

lnvz = fltarr(nfrm)
mix_arr = fltarr(nfrm)

nout=500
;dt = 0.5
tm = dt*nout+findgen(nfrm)*nout*dt

t1 = min(tm)
t2 = max(tm)

t1 = 250
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
dir = '../run2/'
read_data,dir,f,s


;dir = '/Volumes/Scratch/hybrid/KH_new/run_3d_20/'  
dir = '../run3/'
read_data,dir,ff,ss

f = [f,ff]
s = [s,ss]


dir = '../run4/'
read_data,dir,ff,ss

f = [f,ff]
s = [s,ss]

dir = '../run5/'
read_data,dir,ff,ss

f = [f,ff]
s = [s,ss]


dir = '../run6/'
read_data,dir,ff,ss

f = [f,ff]
s = [s,ss]

dir = '../run7/'
read_data,dir,ff,ss

f = [f,ff]
s = [s,ss]




;dir = '/Volumes/Scratch/hybrid/KH_new/run_3d_21/'  
;read_data,dir,ff,ss

;f = [f,ff]
;s = [s,ss]



;dir = '/Volumes/Scratch/hybrid/KH_new/run_3d_24/'  
;read_data,dir,ff,ss

;f = [f,ff]
;s = [s,ss]


;dir = '/Volumes/Scratch/hybrid/KH_new/run_3d_25/'  
;read_data,dir,ff,ss


;f = [f,ff]
;s = [s,ss]

;dir = '/Volumes/Scratch/hybrid/KH_new/run_3d_24/'  
;read_data,dir,f1,s1


;dir = '/Volumes/Scratch/hybrid/KH_new/run_3d_23/'  
;read_data,dir,ff,ss

;f1 = [f1,ff]
;s1 = [s1,ss]


;dir = '/Volumes/Scratch/hybrid/KH_new/run_3d_22/'  
;read_data,dir,ff,ss

;f1 = [f1,ff]
;s1 = [s1,ss]

xarr = [2,3,4,5,6,7]

b = barplot(xarr,f/1e9, nbars = 2, fill_color='red', index=0,name='beta 1.0',xtickdir=1,ytickdir=1)
;b.title='$\beta$ = 1.0'
b.xtitle='run'
b.ytitle='$D (10^9 m^2/s)$'

;b1 = barplot(xarr,f1/1e9, nbars = 2, fill_color='blue',index=1,/overplot,name='beta 0.5')

;l = legend(target=[b,b1])

;l.save,"diffusion_coeff.pdf"

stop
end
