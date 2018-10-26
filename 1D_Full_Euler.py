

import numpy as np
import scipy as sc

#Number of gridpoints
nu = 10
N = 2**nu


#parameters
dx = 0.01
dt = 0.001
test = 1
if test==1:
    gamma = 1.4
    v = 0.0
    p = 8.0/gamma
    rho = 8.0
    x_l = -0.5
    x_r = 0.5
    t_max =0.2
    rho_r=1
    p_r=1
elif test==2:
    # gamma = 1.4
    # v = 0.0
    # p = 8.0 / gamma
    # rho = 8.0
    x_l = 0.0
    x_r = 1.0
    t_max = 3
s=(3,3)
# A = np.zeros(s)
# #conservative variables U
# U = np.array([rho,v*rho,p/(gamma-1)+rho*v**2/2.0])
# A = np.matrix([[0,1,0],[(gamma-3)/2.0*(U[1]/U[0])**2,(3-gamma)*U[1]/U[0],gamma-1],[-gamma*U[1]*U[2]/U[0]**2+(gamma-1)*(U[1]/U[0])**3,gamma*U[2]/U[0]-3/2*(gamma-1)*(U[1]/U[0])**2,gamma*U[1]/U[0]]])
#
# w,r=np.linalg.eig(A)
# print(w)

#primitive variables V
U = np.array([rho,v,p])
A_p = np.matrix([[v,rho,0],[0,v,1.0/rho],[0,gamma*p,v]])
w_l,r_l=np.linalg.eig(A_p)
print(r_p)


#create mesh
grid = np.arange(x_l-dx/2,x_r+dx,dx)
interior_points = range(1,np.shape(grid)[0]-1)

#print(interior_points)
#values on grid n(dt) x m (dx)
t=0
while t<t_max:
    t+=dt
    #make new timestep in matrix
    #set boundary conditions
    #change interior points
    #


