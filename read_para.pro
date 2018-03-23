pro read_para,dir

  close,1
  openr,1,dir+'para.dat',/f77_unformatted

  nx = 0l
  ny = 0l
  nz = 0l
  readu,1,nx,ny,nz,dx,dy,delz
  nt = 0l
  ntsub = 0l
  nout = 0l
  readu,1,nt,dtsub_init,ntsub,dt,nout

  out_dir = '                               '
  readu,1,out_dir
  print,out_dir

  readu,1,vtop,vbottom
  print,vtop,vbottom

  Ni_max = 0l
  readu,1,Ni_max

  readu,1,mproton,m_pu,m_heavy

  readu,1,np_top,np_bottom
  readu,1,b0_top,b0_bottom
  readu,1,vth_top,vth_bottom

  readu,1,alpha,beta

  close,1
  

  save,filename=dir+'para.sav',nx,ny,nz,dx,dy,delz,$
       nt,dtsub_init,ntsub,dt,nout,$
       out_dir,vtop,vbottom,Ni_max,mproton,m_pu,m_heavy,$
       np_top,np_bottom,b0_top,b0_bottom,$
       vth_top,vth_bottom,alpha,beta



return
end
