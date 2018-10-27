

import numpy as np
import scipy as sc

#parameters
dx = 0.01
dt = 0.001

#Number of gridpoints
# nu = 10
# N = 2**nu

test = 1
if test==1:
    x_l = -0.5
    x_r = 0.5
    # create mesh
    grid = np.arange(x_l - dx / 2, x_r + dx, dx)
    interior_points = range(1, np.shape(grid)[0] - 1)
    n_grid = np.shape(grid)[0]
    #initial conditions primitive variables
    gamma = 1.4
    V = np.zeros((1,n_grid))
    P = 8/gamma*np.ones((1,n_grid/2))
    P = np.append(P,np.ones((1,n_grid/2)))
    Rho = np.ones((1, n_grid/2))*8
    Rho = np.append(Rho,np.ones((1, n_grid/2)))
    t_max =0.2

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

# #primitive variables V
# U = np.array([rho,v,p])
# A_p = np.matrix([[v,rho,0],[0,v,1.0/rho],[0,gamma*p,v]])
# w_l,r_l=np.linalg.eig(A_p)
# print(r_l)



#array with solutions
values = np.zeros((3,np.shape(grid)[0],1))
values[0,:,0]=Rho
values[1,:,0]=V
values[2,:,0]=P
print(values)
#print(interior_points)
#values on grid n(dt) x m (dx)
t=0
while t<t_max:
    t+=dt
    #make new timestep in matrix
    #set boundary conditions
    #change interior points
    #


def compute_eigenvalues(U,gamma):
    rho=U[0]
    v = U[1]
    p = U[2]
    A = np.matrix([[v,rho,0],[0,v,1.0/rho],[0,gamma*p,v]])
    e, R = np.linalg.eig(A)
    return {'e':e,'R':R}