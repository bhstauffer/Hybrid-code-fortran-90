close,1


dir = '/Volumes/Scratch/hybrid/KH_new/run_3d_15/'
;rundir = '/Volumes/Scratch/hybrid/KHI/'+'run_test'
;f_read_coord,rundir+'/coord.dat',x,y,z,dzc,dzg,nx,ny,nz
read_para,dir
restore,filename=dir+'para.sav'
read_coords,dir,x,y,z

nfrm=1

;f_read_3d_m_32,rundir+'/temp_p_1',nfrm,tp
;f_read_3d_m_32,rundir+'/npall_1',nfrm,np_t
;f_read_3d_m_32,rundir+'/np_b_1',nfrm,np_b
c_read_3d_vec_m_32,dir,'c.b1',nfrm,b1
c_read_3d_vec_m_32,dir,'c.up',nfrm,up
c_read_3d_m_32,dir,'c.np',nfrm,np
c_read_3d_m_32,dir,'c.mixed',nfrm,mixed
c_read_3d_m_32,dir,'c.temp_p',nfrm,tp



nt = np*tp/1e15 

!p.multi=[0,1,3]
;surface,reform(nt(*,1,*)),charsize=2.0
im = image(reform(tp(*,1,*)), rgb_table=33)
im.scale,2.0,2.0

plot,np(10,1,*)
plot,tp(10,1,*)
plot,(b1(10,1,*,0)^2 + b1(10,1,*,1)^2 + b1(10,1,*,2)^2)

end
