import numpy as np
from scipy.io import FortranFile
import matplotlib.pyplot as plt
import mayavi
from mayavi import mlab
from hybrid_read import Hybrid_read

class plot_hybrid(Hybrid_read):
    def __init__(self):
        super().__init__(dir)
        self.read_para()
        self.read_coords()

    def s_animate_xz(self,file,nfrm,yslc):
        import matplotlib.animation as animation
        fig = plt.figure()
        plt.set_cmap('seismic')
        ims = []
        for i in range(nfrm):
            x = self.read_scalar(file,i+1)
            x[x>1] = 1.0
            im = plt.imshow(x[:,yslc,:].T,origin='lower',animated=True)
            ims.append([im])
        ani = animation.ArtistAnimation(fig,ims[1:],interval=200,blit=True,
                                        repeat_delay=100)
        plt.show()

    def s_animate_vec_xz(self,file,nfrm,yslc,v_comp):
        import matplotlib.animation as animation
        fig = plt.figure()
        plt.set_cmap('seismic')
        ims = []
        for i in range(nfrm):
            x = self.read_vector(file,i+1)
#            x[x>1] = 1.0
            im = plt.imshow(x[:,yslc,:,v_comp].T,origin='lower',animated=True)
            ims.append([im])
        ani = animation.ArtistAnimation(fig,ims[1:],interval=200,blit=True,
                                        repeat_delay=100)
        plt.show()

    def s_plot_iso(self,nfrm):
        
        fig = mlab.figure(size=(800,600),bgcolor=(1,1,1),
                          fgcolor=(0,0,0))
        clr = 'Spectral'
        #mlab.set_cmap('seismic')
        arr = self.read_scalar('c.mixed',nfrm)
        arr[arr>1] = 1.0
        b = self.read_vector('c.b1',nfrm)
        bx = b[:,:,:,0]
        by = b[:,:,:,1]
        bz = b[:,:,:,2]
        x,y,z = np.mgrid[0:np.max(self.x):self.dx,
                         np.min(self.y):np.max(self.y):self.dy,
                         np.min(self.z):np.max(self.z):self.delz]/(self.di/1e3)
        mlab.volume_slice(x,y,z,arr,plane_orientation='y_axes',
                          slice_index=0,colormap=clr)
        mlab.volume_slice(x,y,z,arr,plane_orientation='y_axes',
                          slice_index=np.int(self.ny/2),colormap=clr)
        mlab.volume_slice(x,y,z,arr,plane_orientation='y_axes',
                          slice_index=self.ny-1,colormap=clr)
        mlab.outline(color=(0,0,0))
        mlab.axes(xlabel='x (di)',ylabel='y (di)',zlabel='z (di)',color=(0,0,0))
        #mlab.contour3d(x,y,z,arr,contours=[0.5],
        #               opacity=0.6,transparent=True,colormap='gist_gray')
#        f = mlab.flow(x,y,z,bx,by,bz,scalars=arr,linetype='line',seedtype='plane',line_width=0.25,
#                      seed_resolution=200,seed_scale=30.0,seed_visible=False,
#                      integration_direction='both',colormap=clr,name='flow_plot')
        #mayavi.tools.pipeline.tube(figure = 'f',name='flow_plot',tube_radius=1.0)
        b_lines, topo = self.get_stream('c.b1',nfrm)
        #print(b_lines[0])
        for i in range(len(b_lines)):
            b = b_lines[i]
            xb = b[0,:]
            yb = b[1,:]
            zb = b[2,:]
            ii = np.round((xb*self.di/1e3)/self.dx).astype(int)-1
            jj = np.round((yb*self.di/1e3)/self.dy).astype(int)-1
            kk = np.round((zb*self.di/1e3)/self.delz).astype(int)-1
            s = arr[ii,jj,kk]
            mlab.plot3d(xb,yb,zb,s,tube_radius=0.5,tube_sides=12,colormap='Spectral')

        mlab.show()

    def get_stream(self,file,nfrm):
        import viscid
        from viscid.plot import vlab
        mix = self.read_scalar('c.mixed',nfrm)
        arr = self.read_vector(file,nfrm)
        b = viscid.zeros([self.x,self.y,self.z]/self.di*1e3,nr_comps=3,layout='interlaced')
        X,Y,Z = b.get_crds('xyz',shaped=False)
        b['x'] = arr[1:,1:,1:,0]
        b['y'] = arr[1:,1:,1:,1]
        b['z'] = arr[1:,1:,1:,2]
        obound0 = np.array([X.min(),Y.min(),Z.min()])
        obound1 = np.array([X.max(),Y.max(),Z.max()])
        seeds = viscid.Line((X[1], Y[int(self.ny/2)], Z[int(self.nz/2)]),
                            (X[-2],Y[int(self.ny/2)], Z[int(self.nz/2)]),10)
        b_lines, topo = viscid.calc_streamlines(b,seeds,obound0=obound0,obound1=obound1,
                                                stream_dir = viscid.DIR_BOTH, method=viscid.RK4,
                                                output=viscid.OUTPUT_BOTH)
        
#        vlab.plot_lines(b_lines,scalars=viscid.topology2color(topo))
#        mixed = viscid.project(mix,b)
        #mixarr = viscid.project(mix,b)
#        vlab.streamline(b,seed_resolution=10,seedtype='line',integration_direction='both')
#        print(np.shape(b_lines[0]))
#        help(viscid.calc_streamlines)
                               
#        vlab.show()
        
        return b_lines,topo

    
dir = '/data/hybrid/run_tel_10/'
h = Hybrid_read(dir)

p = plot_hybrid()
#print(p.nx,p.ny,p.nz)
file = 'c.mixed'
p.s_animate_xz(file,10,np.int(p.ny/2))
#file = 'c.up'
#p.s_animate_vec_xz(file,10,np.int(p.ny/2),2)
#print(p.di,np.max(p.x))
#p.s_plot_iso(17)
#p.get_stream('c.b1',17)
