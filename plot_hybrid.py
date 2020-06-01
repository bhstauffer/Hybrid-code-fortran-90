import numpy as np
from scipy.io import FortranFile
import matplotlib.pyplot as plt
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

    def s_plot_iso(self,file,nfrm):
        
        fig = mlab.figure(size=(800,600),bgcolor=(1,1,1),
                          fgcolor=(0,0,0))
        clr = 'Spectral'
        #mlab.set_cmap('seismic')
        for i in range(nfrm):
            arr = self.read_scalar(file,i+1)
            arr[arr>1] = 1.0
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
        mlab.contour3d(x,y,z,arr,contours=[0.5],
                       opacity=0.6,transparent=True,colormap='gist_gray')
        mlab.show()


    def get_stream(self):
        pass

    
dir = './run_va_0.8_beta_3_hires/'
h = Hybrid_read(dir)

p = plot_hybrid()
#print(p.nx,p.ny,p.nz)
file = 'c.mixed'
p.s_animate_xz(file,7,np.int(p.ny/2))
#print(p.di,np.max(p.x))
#p.s_plot_iso(file,17)
