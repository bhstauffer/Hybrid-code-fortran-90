import numpy as np
from scipy.io import FortranFile
import matplotlib.pyplot as plt

class Hybrid_read:

    def __init__(self,dir):
        self.dir = dir
        self.nx = 0
        self.ny = 0
        self.nz = 0
        self.dx = 0.0
        self.dy = 0.0
        self.delz = 0.0
        self.nt = 0
        self.dtsub_init = 0.0
        self.ntsub = 0
        self.dt = 0.0
        self.nout = 0
        self.vtop = 0.0
        self.vbottom = 0.0
        self.Ni_max = 0
        self.mproton = 0.0
        self.m_pu = 0.0
        self.m_heavy = 0.0
        self.np_top = 0.0
        self.np_bottom = 0.0
        self.b0_top = 0.0
        self.b0_bottom = 0.0
        self.vth_top = 0.0
        self.vth_bottom = 0.0
        self.alpha = 0.0
        self.beta = 0.0
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0
        self.va = 0.0
        self.di = 0.0
        
    def read_para(self):
        f = FortranFile(self.dir + 'para.dat','r')
        r = f.read_record('i4','i4','i4','f4','f4','f4')
        self.nx = r[0].item() 
        self.ny = r[1].item()
        self.nz = r[2].item()
        self.dx = r[3].item()
        self.dy = r[4].item()
        self.delz = r[5].item()

        r1 = f.read_record('i4','f4','i4','f4','i4','f4')
        self.nt = r1[0].item()
        self.dtsub_init = r1[1].item()
        self.ntsub = r1[2].item()
        self.dt = r1[3].item()
        self.nout = r1[4].item()

        r2 = f.read_record('a50')
        out_dir = r2[0].item()
        
        r3 = f.read_record('f4','f4')
        self.vtop = r3[0].item()
        self.vbottom = r3[1].item()

        r4 = f.read_record('i4')
        self.Ni_max = r4[0].item()

        r5 = f.read_record('f4','f4','f4')
        self.mproton = r5[0].item()
        self.m_pu = r5[1].item()
        self.m_heavy = r5[2].item()

        r6 = f.read_record('f4','f4')
        self.np_top = r6[0].item()
        self.np_bottom = r6[1].item()

        r7 = f.read_record('f4','f4')
        self.b0_top = r7[0].item()
        self.b0_bottom = r7[1].item()

        r8 = f.read_record('f4','f4')
        self.vth_top = r8[0].item()
        self.vth_bottom = r8[1].item()

        r9 = f.read_record('f4','f4')
        self.alpha = r9[0].item()
        self.beta = r9[1].item()

        q = 1.67e-19
        mu0 = np.pi*4e-7
        ep = 8.85e-12
        c = 3e8 #m/s
        
        self.va = self.b0_bottom/np.sqrt(mu0*self.mproton*self.np_bottom/1e9)
        wpi = np.sqrt(q*q*self.np_bottom/1e9/(ep*self.mproton))
        self.di = c/wpi
        
        f.close()

    def read_coords(self):
        f = FortranFile(self.dir + 'c.coord.dat','r')
        nnx = f.read_ints('i4').item()
        nny = f.read_ints('i4').item()
        nnz = f.read_ints('i4').item()
        #print(nnx,nny,nnz)
        self.x = f.read_reals('f4')
        self.y = f.read_reals('f4')
        self.z = f.read_reals('f4')
    
        f.close()

    def read_scalar(self,file,nfrm):
        file = self.dir+file+'.dat'
        f = FortranFile(file,'r')
        #x = np.zeros(self.nx,self.ny,self.nz)
        print(self.nx,self.ny,self.nz)
        for _ in range(nfrm):        
            frm = f.read_ints('i4').item()
            print(frm)
            x = f.read_reals('f4').reshape((self.nx,self.ny,self.nz),order = 'F')
        f.close()
        return x

    def read_vector(self,file,nfrm):
        file = self.dir+file+'.dat'
        f = FortranFile(file,'r')
        for _ in range(nfrm):        
            frm = f.read_ints('i4').item()
            print(frm)
            x = f.read_reals('f4').reshape((self.nx,self.ny,self.nz,3),order = 'F')
        f.close()
        return x

    def read_part(self,file,nfrm):
        file = self.dir+file+'.dat'
        f = FortranFile(file,'r')
        for _ in range(nfrm):
            frm = f.read_ints('i4').item()
            print(frm)
            x = f.read_reals('f4').reshape((self.Ni_max,3),order = 'F')
        f.close()
        return x

    def read_part_scalar(self,file,nfrm):
        file = self.dir+file+'.dat'
        f = FortranFile(file,'r')
        for _ in range(nfrm):
            frm = f.read_ints('i4').item()
            print(frm)
            x = f.read_reals('f4').reshape((self.Ni_max),order = 'F')
        f.close()
        return x
            
    def write_hdf5(self,filename,arr):
        import h5py
        self.read_coords()
        self.read_para()
        with h5py.File(filename+'.hdf5','w') as f:
            dt = f.create_dataset("dt", data = self.dt)
            nout = f.create_dataset("nout", data = self.nout)
            nx = f.create_dataset("nx", data = self.nx)
            ny = f.create_dataset("ny", data = self.ny)
            nz = f.create_dataset("nz", data = self.nz)
            x = f.create_dataset("x", data = self.x)
            y = f.create_dataset("y", data = self.y)
            z = f.create_dataset("z", data = self.z)
            dset = f.create_dataset("array",data = arr,compression="gzip")
    
#---------------------------------------------------------        
#dir = './run_va_0.8_beta_1/'
#h = Hybrid_read(dir)
#h.read_para()
#h.read_coords()
#mix = h.read_scalar('c.mixed',20)
#for i in range(1,25):
#    b1 = h.read_vector('c.E',i)
#    b1 = b1*h.mproton*1e3/1.6e-19
#    print(i,b1.max())
#    h.write_hdf5('./hdf5/electric_field_'+str(i),b1)

#mrat = h.read_part_scalar('c.mrat_ 1',1)
      
#plt.imshow(mix[:,np.int(h.ny/2),:].T,origin='lower')
#plt.show()

#plt.imshow(up[:,np.int(h.ny/2),:,2].T,origin='lower')
#plt.show()
