import numpy as np
from structure_tensor import eig_special_3d, structure_tensor_3d

data = .... #Seismic cube here

sigma = 1.
rho = 5.

S = structure_tensor_3d(data, sigma, rho)
val, vec = eig_special_3d(S.astype('float64'), full=True)

dt = vec[-1,0,:,:,:]
dy = vec[-1,1,:,:,:]
dx = vec[-1,2,:,:,:]

mask = dt<0
dt[mask]*=-1
dy[mask]*=-1
dx[mask]*=-1

dip_x = dx/dt
dip_y = dy/dt
coherency = (val[2,:,:,:])/(val[2,:,:,:]+val[1,:,:,:]+val[0,:,:,:])
