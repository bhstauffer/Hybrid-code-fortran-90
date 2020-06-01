import numpy as np
import matplotlib.pyplot as plt
from hybrid_read import Hybrid_read

dir = './run_va_0.8_beta_3/'
h = Hybrid_read(dir)

nfrm = 18
h.read_para()
h.read_coords()
b1 = h.read_vector('c.b1',nfrm)
up = h.read_vector('c.up',nfrm)
ni = h.read_scalar('c.np',nfrm)

j = np.int(h.ny/2+10)
k = np.int(h.nz/2-5)
q = 1.6e-19

b1 = b1*h.mproton/q
mu0 = np.pi*4e-7

va = b1[:,j,k,0]/np.sqrt(mu0*h.mproton*ni[:,j,k]/1e9)

fig,ax = plt.subplots()
ax.plot(h.x,(va/1e3),'r')
ax.set_xlabel('x')
ax.set_ylabel('vA (km/s)')
ax.set_ylim(-100,100)
ax2 = ax.twinx()
ax2.plot(h.x,(up[:,j,k,0]),'b')
ax2.set_ylabel('up (km/s)')
ax2.set_ylim(-100,100)
plt.show()
