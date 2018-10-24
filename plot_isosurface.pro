;-----------------------------------------------------------------------
pro swap_arr_s,a,nx,ny,nz
;-----------------------------------------------------------------------

b = fltarr(nx,nz,ny)
for k = 0,nz-1 do begin
   for j = 0,ny-1 do begin
      b(*,k,j) = a(*,j,k)      
   endfor
endfor
a = b

return
end
;-----------------------------------------------------------------------

;-----------------------------------------------------------------------
pro swap_arr_v,a,nx,ny,nz
;-----------------------------------------------------------------------

b = fltarr(nx,nz,ny,3)
for k = 0,nz-1 do begin
   for j = 0,ny-1 do begin
      for m = 0,2 do begin
         b(*,k,j,m) = a(*,j,k,m)      
      endfor
   endfor
endfor
a = b

return
end
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
pro get_edotb,a,b,c
;-----------------------------------------------------------------------
;a = efld
;b = bfld

;a2 = sqrt(a(*,*,*,0)^2 + a(*,*,*,1)^2 + a(*,*,*,2)^2)
;b2 = sqrt(b(*,*,*,0)^2 + b(*,*,*,1)^2 + b(*,*,*,2)^2)

;for m = 0,2 do begin
;   a(*,*,*,m) = a(*,*,*,m)/a2(*,*,*)
;   b(*,*,*,m) = b(*,*,*,m)/b2(*,*,*)
;endfor

sz = size(a)
nz = sz(3)
a2 = sqrt(a(*,*,*,0)^2 +a(*,*,*,1)^2+a(*,*,*,2)^2)
ahat = a
;ahat(*,*,0)= ahat(*,*,1)
ahat(*,*,nz-1)= ahat(*,*,nz-2)
for m = 0,2 do begin
   ahat(*,*,*,m) = a(*,*,*,m);/a2(*,*,*)
endfor

b2 = sqrt(b(*,*,*,0)^2 +b(*,*,*,1)^2+b(*,*,*,2)^2)
;stop
bhat = b
;bhat(*,*,0)= bhat(*,*,1)
bhat(*,*,nz-1)= bhat(*,*,nz-2)
for m = 0,2 do begin
   bhat(*,*,*,m) = b(*,*,*,m)/b2(*,*,*)
endfor

c(*,*,*) = ahat(*,*,*,0)*bhat(*,*,*,0) + ahat(*,*,*,1)*bhat(*,*,*,1) + ahat(*,*,*,2)*bhat(*,*,*,2)

;c = abs(c)

;c = c/max(c)

return
end
;-----------------------------------------------------------------------

;-----------------------------------------------------------------------
pro plot_isosurface,nfrm,LINES=lines, TUBES=tubes
;-----------------------------------------------------------------------

temp_val = 300 ;eV
den_val = 0.5 ;c
b_val = 0.0e-9
edotb_val = 0.025

Rio = 1200.
ctbl = 13
ach = 0.7

;file_dir = '/Volumes/MacD97-2/hybrid/3d_buf/run_test/'

;dir = '/Volumes/Scratch/hybrid/KH_new/run_3d_9/'
dir = './run2_3d/'
;dir = '/Volumes/Scratch/hybrid/KH3D/run_3_periodic/'
;dir = './tmp1/'

mrestart=''

read_para,dir
restore,filename=dir+'para.sav'
;stop
read_coords,dir,x,y,z
dx = x(1)-x(0)
dy = y(1)-y(0)
dz = z(1)-z(0)
x = x/dx
y = y/dy
z = z/dz

@get_const
c_read_3d_vec_m_32,dir,'c.b1'+mrestart,nfrm,b1
c_read_3d_vec_m_32,dir,'c.up'+mrestart,nfrm,up
c_read_3d_vec_m_32,dir,'c.E'+mrestart,nfrm,Efld
c_read_3d_m_32,dir,'c.np'+mrestart,nfrm,np
c_read_3d_m_32,dir,'c.mixed'+mrestart,nfrm,mixed
c_read_3d_m_32,dir,'c.temp_p'+mrestart,nfrm,temp_p

wh = where(mixed gt 1.0)
mixed(wh) = 1.0

edotb = fltarr(nx,ny,nz)
Efld(*,*,0,*) = Efld(*,*,1,*)
Efld(*,*,nz-1,*) = Efld(*,*,nz-2,*)
b1(*,*,0,*) = b1(*,*,1,*)
b1(*,*,nz-1,*) = b1(*,*,nz-2,*)
Efld0 = abs(vtop)*b0_top*1.6e-19/mproton
get_edotb,Efld,b1,edotb
edotb = edotb/Efld0

;swap y and z directions

swap_arr_v,b1,nx,ny,nz
swap_arr_v,up,nx,ny,nz

swap_arr_s,mixed,nx,ny,nz
swap_arr_s,np,nx,ny,nz
swap_arr_s,temp_p,nx,ny,nz
swap_arr_s,edotb,nx,ny,nz

xx = x
yy = z
zz = y

x = xx
y = yy
z = zz

nnx = nx
nny = nz
nnz = ny
nz = nnz
ny = nny

np = mixed

minx = min(x)
maxx = max(x)
miny = min(y)
maxy = max(y)
minz = min(z)
maxz = max(z)

oModel = OBJ_NEW('IDLgrModel')  

xax = OBJ_NEW('IDLgrAxis',direction=0,range=[minx,maxx],location=[0,0,0],$
              /notext,/exact,thick=3,ticklen=1)
yax = OBJ_NEW('IDLgrAxis',direction=1,range=[miny,maxy],location=[0,0,0],$
              /notext,/exact,thick=3,ticklen=1)
zax = OBJ_NEW('IDLgrAxis',direction=2,range=[minz,maxz],location=[0,0,0],$
              /notext,/exact,thick=3,ticklen=1)

xtit = OBJ_NEW('IDLgrText','x',baseline = [1.0,0,0], $
               locations=[nx-5,-10,0],$
               alignment=0.5,/enable_formatting,color=!color.black,$
              char_dimensions = [10,10])

ytit = OBJ_NEW('IDLgrText','z',baseline = [0,1.0,0], $
               updir = [-1.0,0,0], locations=[-10,ny-5,0],$
               char_dimensions=[10,10],$
               alignment=0.5,/enable_formatting,color=!color.black)


ztit = OBJ_NEW('IDLgrText','y',baseline = [1.0,0,0.0], $
               updir = [0.0,0.0,1.0],locations=[-10,0,nz-5],$
               char_dimensions=[10,10],$
               alignment=0.5,/enable_formatting,color=!color.black)


oPalette = OBJ_NEW('IDLgrPalette')  
oPalette -> LOADCT, 33  
oPalette -> SetRGB, 255, 255, 255, 255  


oPalette2 = OBJ_NEW('IDLgrPalette')
oPalette2->LOADCT, 33
oPalette2 -> SetRGB, 255, 255, 255, 255  

;oImage = OBJ_NEW('IDLgrImage', image, PALETTE = oPalette)  

vector = FINDGEN(101)/100.  
;vector = shift(vector,20)
texure_coordinates = FLTARR(2, 101, 101)  
texure_coordinates[0, *, *] = vector # shift(REPLICATE(1., 101),0)  
texure_coordinates[1, *, *] = REPLICATE(1., 101) # vector  

oPolygons = OBJ_NEW('IDLgrPolygon', SHADING = 1, $  
   DATA = vertices, POLYGONS = polygons, $  
   COLOR = [255, 255, 255], $  
   TEXTURE_COORD = texure_coordinates, $  
   TEXTURE_MAP = oImage, /TEXTURE_INTERP)  

data = FLTARR(3,n_elements(x), n_elements(y), n_elements(z))

;data(0,*,*,*) = smooth(b1(*,*,*,0),2)
;data(1,*,*,*) = smooth(b1(*,*,*,1),2)
;data(2,*,*,*) = smooth(b1(*,*,*,2),2)

data(0,*,*,*) = b1(*,*,*,0)
data(2,*,*,*) = b1(*,*,*,1)
data(1,*,*,*) = b1(*,*,*,2)


xstep = 2
zstep = 2
nseeds = 3*LONG((nx*ny)/(xstep*zstep))
seeds = FLTARR(nseeds)
iseed=0L
;whx = where(x gt minx)

;for k = nz/2+1,nz/2+2,4 do begin
;for j = ny/2-2,ny/2+2,3 do begin
for j = 52,56,3 do begin
;j = 29
   for i = 10,nx-10,4 do begin
      print,i,j
      if ((x(i) gt x(1)) and (x(i) lt max(x)) $
;         and (z(k) gt z(1)) and (z(k) lt max(z))) then begin
         and (y(j) gt y(1)) and (y(j) lt max(y))) then begin
         seeds(iseed) = FLOAT(i)
         seeds(iseed+1) = FLOAT(j)
         seeds(iseed+2) = FLOAT(0)
;         print,seeds(iseed:iseed+3),iseed
         iseed = iseed+3
;         print,iseed,i,k
      endif
   endfor
;endfor

maxIterations=4000
stepSize=0.4

PARTICLE_TRACE,data,seeds(0:iseed-1),outverts,outconn,outnormals, $
               MAX_ITERATIONS=maxIterations, MAX_STEPSIZE=stepSize,  $
               INTEGRATION=1,ANISOTROPY=[1.0,1.0,1.0], SEED_NORMAL=[1, 0, 0]

;outnormals(0,*)= median(reform(outnormals(0,*)),5)
;outnormals(1,*)= median(reform(outnormals(1,*)),5)
;outnormals(2,*)= median(reform(outnormals(2,*)),5)

mon = sqrt(outnormals(0,*)^2 + outnormals(1,*)^2 + outnormals(2,*)^2)
outnormals(0,*) = outnormals(0,*)/mon
outnormals(1,*) = outnormals(1,*)/mon
outnormals(2,*) = outnormals(2,*)/mon


magdata = SQRT(data[0,*, *,*]^2 + 0.0*data[1,*, *,*]^2 + 0.0*data[2,*, *,*]^2)
vertX =  REFORM(outverts[0,*],N_ELEMENTS(outverts)/3)
vertY =  REFORM(outverts[1,*],N_ELEMENTS(outverts)/3)
vertZ =  REFORM(outverts[2,*],N_ELEMENTS(outverts)/3)
vertcolors = BYTSCL(INTERPOLATE(magdata,vertX, vertY, vertZ))

endfor

if(KEYWORD_SET(tubes)) then begin    ;square profile for stream-tubes.
   width=2.0
;   profile = [[-1,-1],[-1,1],[1,1],[1,-1],[-1,-1]]
   profile=fltarr(2,10)
   for i = 0,9 do begin
      profile(0,i) = cos(i*40*!dtor)
      profile(1,i) = sin(i*40*!dtor)
   endfor
   profile(1,*) = 0.5*stepsize*profile(1,*)

   nverts = N_ELEMENTS(outverts)/3
   STREAMLINE, TEMPORARY(outverts),TEMPORARY(outconn), $
               outnormals*width,outverts,outconn, PROFILE=profile

   magdata = SQRT(data[0,*, *,*]^2 + data[1,*, *,*]^2 + data[2,*, *,*]^2)
   vertX =  REFORM(outverts[0,*],N_ELEMENTS(outverts)/3)
   vertY =  REFORM(outverts[1,*],N_ELEMENTS(outverts)/3)
   vertZ =  REFORM(outverts[2,*],N_ELEMENTS(outverts)/3)
;   vertcolors = BYTSCL(INTERPOLATE(magdata,vertX, vertY, vertZ))
   vertcolors = BYTSCL(INTERPOLATE(smooth(mixed<1.0,2),vertX, vertY, vertZ))
;   vertcolors = BYTSCL(INTERPOLATE(reform(up(*,*,*,0)),vertX, vertY, vertZ))

;   ovx = reform(outverts(0,*))
;   ovy = reform(outverts(1,*))
;   ovz = reform(outverts(2,*))
   
;;truncated float to interger (e.g. 2.56 = 2)
;   fovx = fix(ovx)
;   fovy = fix(ovy)
;   fovz = fix(ovz)
   
;;interpolate floating point indices to intermediate grid locations
;;switch x and z directions
;   px = (ovx - fovx)*(x(fovx+1)-x(fovx)) + x(fovx) 
;   py = (ovy - fovy)*(y(fovy+1)-y(fovy)) + y(fovy) 
;   pz = (ovz - fovz)*(z(fovz+1)-z(fovz)) + z(fovz) 
   
;   outverts(0,*) = px
;   outverts(1,*) = py
;   outverts(2,*) = pz
   
   
   oStreamlines = OBJ_NEW('IDLgrPolygon',outverts, POLYGONS=outconn, $
                          SHADING = 1,alpha_channel=1.0)



endif


if(KEYWORD_SET(lines)) then begin $   ;square profile for stream-tubes.

;   ovx = reform(outverts(0,*))
;   ovy = reform(outverts(1,*))
;   ovz = reform(outverts(2,*))
   
;;truncated float to interger (e.g. 2.56 = 2)
;   fovx = fix(ovx)
;   fovy = fix(ovy)
;   fovz = fix(ovz)
   
;;interpolate floating point indices to intermediate grid locations
;;switch x and z directions
;   px = (ovx - fovx)*(x(fovx+1)-x(fovx)) + x(fovx) 
;   py = (ovy - fovy)*(y(fovy+1)-y(fovy)) + y(fovy) 
;   pz = (ovz - fovz)*(z(fovz+1)-z(fovz)) + z(fovz) 
   
;   outverts(0,*) = px
;   outverts(1,*) = py
;   outverts(2,*) = pz
   
   outverts=MESH_SMOOTH(outverts,outconn)

   oStreamlines = OBJ_NEW('IDLgrPolyline',outverts, $
                          POLYLINES=outconn) ;, $
;                       LABEL_OBJECTS=oLblSymbols, $
;                       LABEL_POLYLINES=lblPolys, $
;                       LABEL_OFFSETS=lblOffsets, $
;                       /LABEL_USE_VERTEX_COLOR, $
;                       /LABEL_NOGAPS)

endif


oStreamlines->SetProperty, PALETTE = oPalette, VERT_COLORS = vertcolors

;oModel->Add, oStreamlines 

ntot = mixed
ntot = smooth(ntot,2)
;isosurface,ntot(*,*,*)/1e15,den_val,$
;           outverts,outconn
isosurface,ntot(*,*,*),den_val,$
           outverts,outconn

outverts=MESH_SMOOTH(outverts,outconn)

;ovx = reform(outverts(0,*))
;ovy = reform(outverts(1,*))
;ovz = reform(outverts(2,*))

;;truncated float to interger (e.g. 2.56 = 2)
;fovx = fix(ovx)
;fovy = fix(ovy)
;fovz = fix(ovz)

;px = (ovx - fovx)*(x(fovx+1)-x(fovx)) + x(fovx) 
;py = (ovy - fovy)*(y(fovy+1)-y(fovy)) + y(fovy) 
;pz = (ovz - fovz)*(z(fovz+1)-z(fovz)) + z(fovz) 

;outverts(0,*) = px
;outverts(1,*) = py
;outverts(2,*) = pz

oden = OBJ_NEW('IDLgrPolygon',outverts, POLYGONS=outconn, $
                          SHADING = 1,alpha_channel=0.7,palette=opalette,$
			color=!color.grey)

;oModel->Add, oden

ti = temp_p>0
ti = smooth(ti,2)
isosurface,ti(*,*,*),temp_val,$
           outverts,outconn


outverts=MESH_SMOOTH(outverts,outconn)

;ovx = reform(outverts(0,*))
;ovy = reform(outverts(1,*))
;ovz = reform(outverts(2,*))

;;truncated float to interger (e.g. 2.56 = 2)
;fovx = fix(ovx)
;fovy = fix(ovy)
;fovz = fix(ovz)

;px = (ovx - fovx)*(x(fovx+1)-x(fovx)) + x(fovx) 
;py = (ovy - fovy)*(y(fovy+1)-y(fovy)) + y(fovy) 
;pz = (ovz - fovz)*(z(fovz+1)-z(fovz)) + z(fovz) 

;outverts(0,*) = px
;outverts(1,*) = py
;outverts(2,*) = pz

otemp = OBJ_NEW('IDLgrPolygon',outverts, POLYGONS=outconn, $
                          SHADING = 1,alpha_channel=ach,palette=opalette,$
	  color=!color.red)

edotb = smooth(edotb,2)
;edotb = abs(edotb)

isosurface,edotb(2:nx-2,*,*),edotb_val,$
           outverts,outconn

sz = size(outverts)
if(sz(0) gt 1) then begin
outverts=MESH_SMOOTH(outverts,outconn)

;ovx = reform(outverts(0,*))
;ovy = reform(outverts(1,*))
;ovz = reform(outverts(2,*))

;;truncated float to interger (e.g. 2.56 = 2)
;fovx = fix(ovx)
;fovy = fix(ovy)
;fovz = fix(ovz)

;px = (ovx - fovx)*(x(fovx+1)-x(fovx)) + x(fovx) 
;py = (ovy - fovy)*(y(fovy+1)-y(fovy)) + y(fovy) 
;pz = (ovz - fovz)*(z(fovz+1)-z(fovz)) + z(fovz) 

;outverts(0,*) = px
;outverts(1,*) = py
;outverts(2,*) = pz

oedotb_p = OBJ_NEW('IDLgrPolygon',outverts, POLYGONS=outconn, $
                          SHADING = 1,alpha_channel=ach,palette=opalette,$
	  color=!color.red)

oModel->Add, oedotb_p
endif

isosurface,edotb(2:nx-2,*,*),-edotb_val,$
           outverts,outconn

sz = size(outverts)
if(sz(0) gt 1) then begin
outverts=MESH_SMOOTH(outverts,outconn)

;ovx = reform(outverts(0,*))
;ovy = reform(outverts(1,*))
;ovz = reform(outverts(2,*))

;;truncated float to interger (e.g. 2.56 = 2)
;fovx = fix(ovx)
;fovy = fix(ovy)
;fovz = fix(ovz)

;px = (ovx - fovx)*(x(fovx+1)-x(fovx)) + x(fovx) 
;py = (ovy - fovy)*(y(fovy+1)-y(fovy)) + y(fovy) 
;pz = (ovz - fovz)*(z(fovz+1)-z(fovz)) + z(fovz) 

;outverts(0,*) = px
;outverts(1,*) = py
;outverts(2,*) = pz

oedotb_m = OBJ_NEW('IDLgrPolygon',outverts, POLYGONS=outconn, $
                          SHADING = 1,alpha_channel=ach,palette=opalette,$
	  color=!color.red)

oModel->Add, oedotb_m
endif


;oModel->Add, otemp

bx = reform(smooth(b1(*,*,*,2),2))*mp/q

isosurface,(bx(*,*,*)),b_val,$
           outverts,outconn


outverts=MESH_SMOOTH(outverts,outconn)

;ovx = reform(outverts(0,*))
;ovy = reform(outverts(1,*))
;ovz = reform(outverts(2,*))

;;truncated float to interger (e.g. 2.56 = 2)
;fovx = fix(ovx)
;fovy = fix(ovy)
;fovz = fix(ovz)

;px = (ovx - fovx)*(x(fovx+1)-x(fovx)) + x(fovx) 
;py = (ovy - fovy)*(y(fovy+1)-y(fovy)) + y(fovy) 
;pz = (ovz - fovz)*(z(fovz+1)-z(fovz)) + z(fovz) 

;outverts(0,*) = px
;outverts(1,*) = py
;outverts(2,*) = pz

obx = OBJ_NEW('IDLgrPolygon',outverts, POLYGONS=outconn, $
                          SHADING = 1,alpha_channel=0.7,palette=opalette,$
			color=!color.green)


by = reform(smooth(b1(*,*,*,2),2))*mp/q

isosurface,by(*,*,*),-b_val,$
           outverts,outconn


outverts=MESH_SMOOTH(outverts,outconn)

;ovx = reform(outverts(0,*))
;ovy = reform(outverts(1,*))
;ovz = reform(outverts(2,*))

;;truncated float to interger (e.g. 2.56 = 2)
;fovx = fix(ovx)
;fovy = fix(ovy)
;fovz = fix(ovz)

;px = (ovx - fovx)*(x(fovx+1)-x(fovx)) + x(fovx) 
;py = (ovy - fovy)*(y(fovy+1)-y(fovy)) + y(fovy) 
;pz = (ovz - fovz)*(z(fovz+1)-z(fovz)) + z(fovz) 

;outverts(0,*) = px
;outverts(1,*) = py
;outverts(2,*) = pz

oby = OBJ_NEW('IDLgrPolygon',outverts, POLYGONS=outconn, $
                          SHADING = 1,alpha_channel=0.7,palette=opalette,$
			color=!color.red)

;oModel->Add, oby

wh = where(mixed gt 1.0)
mixed(wh) = 1.0
;mixed = mixed+1.0
mixed(0,0,*) = -0.001


oCont = OBJ_NEW('IDLgrContour',reform(mixed(*,*,nz/2)), geomx = x, geomy = y,geomz = nz/2,$
                planar=1,fill=1,n_levels=255,palette=oPalette2,shading=1,$
               alpha_channel=ach,c_color=bytscl(indgen(255)))

oCont1 = OBJ_NEW('IDLgrContour',reform(mixed(*,*,0)), geomx = x, geomy = y,geomz = 0,$
                planar=1,fill=1,n_levels=255,palette=oPalette2,shading=1,$
               alpha_channel=ach,c_color=bytscl(indgen(255)))

oCont2 = OBJ_NEW('IDLgrContour',reform(mixed(*,*,nz-1)), geomx = x, geomy = y,geomz = nz-1,$
                planar=1,fill=1,n_levels=255,palette=oPalette2,shading=1,$
               alpha_channel=ach,c_color=bytscl(indgen(255)))


oContup = OBJ_NEW('IDLgrContour',reform(up(*,*,nz/2,0)), geomx = x, geomy = y,geomz = nz/2,$
                planar=1,fill=1,n_levels=255,palette=oPalette2,shading=1,$
               alpha_channel=ach,c_color=bytscl(indgen(255)))

oContup1 = OBJ_NEW('IDLgrContour',reform(up(*,*,0,0)), geomx = x, geomy = y,geomz = 0,$
                planar=1,fill=1,n_levels=255,palette=oPalette2,shading=1,$
               alpha_channel=ach,c_color=bytscl(indgen(255)))

;;add vector field
;v_arr = fltarr(2,nx,ny)
;v_arr(0,*,*) = up(*,*,1,0)
;v_arr(1,*,*) = up(*,*,1,1)
;vert = fltarr(2,nx*ny/8)
;cnt = 0
;for i = 0,nx-1,8 do begin
;   for j = 0,ny-1,8 do begin
;      vert(0,cnt) = i
;      vert(1,cnt) = j
;      cnt = cnt+1
;   endfor
;endfor

;vector_field,v_arr,outverts,outconn,scale=0.04,vertices=vert
;ovect = OBJ_NEW('IDLgrPolyline',DATA=outverts,POLYLINES=outconn,COLOR=[255,255,255],THICK=4.0)


;add all objects here---------------------------------------------

oModel->Add, xtit
oModel->Add, ytit
oModel->Add, ztit
oModel->Add,xax
oModel->Add,yax
oModel->Add,zax
;oModel -> ADD, oPolygons  ;this is Io
oModel->Add, oStreamlines 
;oModel->Add, oI24
;oModel->Add, oI27
;oModel->Add, oI31
;oModel->Add, oJ0

;oModel->Add, otemp
;oModel->Add, obx
;oModel->Add, oby

oModel->Add, oCont
oModel->Add, oCont1
oModel->Add, oCont2
;oModel->Add, ovect
;oModel->Add, oContup
;oModel->Add, oContup1
oModel->Add, oden

;-----------------------------------------------------------------


oModel->Rotate, [1,0,0], -90
oModel->Rotate, [0,1,0], 30*10.5
oModel->Rotate, [1,0,0], 30

oLight = OBJ_NEW('IDLgrLight', loc=[100,-100,0], Inten=1.0, type=2)
oModel->Add,oLight

aLight = OBJ_NEW('IDLgrLight', Inten=1.0, type=0)
oModel->Add,aLight


view = obj_new('IDLgrView',color=[0,0,0],dimensions=[0,0]);,eye=9,$
; 	viewplane_rect=10*[-3.5,-3.5,7.0,7.0],zclip = 10*[7,-1000])
view->add, oModel

;win = obj_new('IDLgrWindow',dimensions=[1200,1200], graphics_tree=view)

;win->draw, view
;win->GetProperty, image_data=im

;write_image,'test','png',im,im(0,*,*),im(1,*,*),im(2,*,*)

;XINTERANIMATE, SET=[1200,1200,361], /SHOWLOAD


;for i = 0,360/4. do begin
;;for i = 0,n_elements(x31)-2 do begin

;   omodel->rotate,[0,1,0],1*4.

;;   dx = x31(i+1)-x31(i)
;;   dy = y31(i+1)-y31(i)
;;   dz = z31(i+1)-z31(i)
;;   oModel->translate,-dx,-dy,-dz
;;view->add, oModel
;   win->draw, view
;   win->GetProperty, image_data=im
;   xinteranimate, frame = i, image = im
   
;endfor
;xinteranimate,/keep_pixmaps

XOBJVIEW, oModel, scale=1.0,/BLOCK,xsize=1200,ysize=1200
;xOBJVIEW, oModel,scale=0.9,xsize=1200,ysize=1200
;xobjview_write_image,'test.png','png'


OBJ_DESTROY, [oModel, oPalette, oPalette2]

end















