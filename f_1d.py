from __future__ import print_function
#https://github.com/ibackus/sod-shocktube
import sod
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt


def main():
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
        gamma = 1.4
        # v = 0.0
        # p = 8.0 / gamma
        # rho = 8.0
        x_l = 0.0
        x_r = 1.0
        # create mesh
        grid = np.arange(x_l - dx / 2, x_r + dx, dx)
        interior_points = range(1, np.shape(grid)[0] - 1)
        n_grid = np.shape(grid)[0]
        V = np.zeros((1, n_grid))
        P = 0.1 + 0.001 * np.exp(-(grid - 0.5) ** 2 / 0.1 ** 2)
        Rho = gamma * P
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
    #print(values)
    #print(interior_points)
    #values on grid n(dt) x m (dx)
    t=0




    while t<t_max:
        #make new timestep in matrix
        #set boundary conditions
        #change interior points
        Rs = np.array([[[]]])
        values = np.append(values, np.zeros((3, n_grid, 1)), axis=2)

        for i in xrange(0,n_grid):
            if (i==0):
                if test==2:
                    decom = compute_eigenvalues(values[:, i, -2], gamma)
                    diagsi = convert_to_diag(values[:, i, -2], decom['L'])
                    # if (i>0):
                    diagsimin = convert_to_diag(values[:, i - 1, -2], decom['L'])
                    diagsiplus = convert_to_diag(values[:, i + 1, -2], decom['L'])

                    verschil = np.zeros((3, 1))
                    for j in xrange(0, 3):
                        if (decom['e'][j] >= 0):
                            verschil[j][0] = diagsi[j] - diagsimin[j]
                        elif (decom['e'][j] < 0):
                            verschil[j][0] = diagsiplus[j] - diagsi[j]
                    z = np.array([np.matmul(decom['L'], values[:, i, -2])]).T - dt / (dx) * np.multiply(
                        np.array([decom['e']]).T, verschil)
                    # print(decom['R'])
                    values[:, i, -1] = np.matmul(decom['R'], z[:, 0])

            elif(i==n_grid-1):
                if test == 2:
                    decom = compute_eigenvalues(values[:, i, -2], gamma)
                    diagsi = convert_to_diag(values[:, i, -2], decom['L'])
                    # if (i>0):
                    diagsimin = convert_to_diag(values[:, i - 1, -2], decom['L'])
                    diagsiplus = convert_to_diag(values[:, 0, -2], decom['L'])

                    verschil = np.zeros((3, 1))
                    for j in xrange(0, 3):
                        if (decom['e'][j] >= 0):
                            verschil[j][0] = diagsi[j] - diagsimin[j]
                        elif (decom['e'][j] < 0):
                            verschil[j][0] = diagsiplus[j] - diagsi[j]
                    z = np.array([np.matmul(decom['L'], values[:, i, -2])]).T - dt / (dx) * np.multiply(
                        np.array([decom['e']]).T, verschil)
                    # print(decom['R'])
                    values[:, i, -1] = np.matmul(decom['R'], z[:, 0])

            else:
                decom = compute_eigenvalues(values[:,i,-2],gamma)
                diagsi = convert_to_diag(values[:,i,-2],decom['L'])
                #if (i>0):
                diagsimin = convert_to_diag(values[:,i-1,-2],decom['L'])
                diagsiplus = convert_to_diag(values[:,i+1,-2],decom['L'])

                verschil = np.zeros((3,1))
                for j in xrange(0,3):
                    if (decom['e'][j]>=0):
                        verschil[j][0] = diagsi[j]-diagsimin[j]
                    elif (decom['e'][j]<0):
                        verschil[j][0] = diagsiplus[j]-diagsi[j]
                z = np.array([np.matmul(decom['L'],values[:,i,-2])]).T-dt/(dx)*np.multiply(np.array([decom['e']]).T,verschil)
                #print(decom['R'])
                values[:,i,-1] = np.matmul(decom['R'],z[:,0])
            # if (i==0):
            #     Rs = np.array([decom['R']])
            #
            #     diag = np.matrix([diagsi]).T
            # else:
            #     Rs = np.append(Rs, np.array([decom['R']]), axis=0)
            #     diag=np.append(diag,np.matrix([diagsi]).T,axis=1)
            # if (i == n_grid-1):
            #     z = diag[:,1:-1] -get_flux(decom['e'], diag, dx, dt)[:,0:-1]
        if test == 1:
            values[:, 0, -1] = values[:, 1, -1]
            values[:, -1, -1] = values[:, -2, -1]

        t += dt
        print(t/t_max,'%',end='\r')
        # for i in xrange(1,n_grid-1):
        #     if(i==1):
        #         prim = convert_to_prim(z[:,i-1],Rs[i,:,:])
        #     else:
        #         prim = np.append(prim,convert_to_prim(z[:,i-1],Rs[i,:,:]),axis=1)
        #
        # #keep boundary conditions
        #
        # values[:,1:-1,-1] = prim
        #print(decom['R']* decom['L'])
    #values.tofile('test2x00t001.dat')
    plt.figure(1)
    plt.plot(grid,values[0,:,-1])
    plt.ylabel('density')
    plt.xlabel('x')

    plt.figure(2)
    plt.plot(grid,values[1, :, -1])
    plt.ylabel('velocity')
    plt.xlabel('x')

    plt.figure(3)
    plt.plot(grid,values[2, :, -1])
    plt.ylabel('pressure')
    plt.xlabel('x')

    if test==1:
        analytic = analytic_solution(gamma,n_grid)
        p = analytic['p']
        rho = analytic['rho']
        u = analytic['u']
        x_analytic = analytic['x']
        plt.figure(1)
        plt.plot(x_analytic, rho, color='g')
        plt.ylabel('pressure')
        plt.xlabel('x')
        #plt.axis([0, 1, 0, 1.1])

        plt.figure(2)
        plt.plot(x_analytic, u, color='g')
        plt.ylabel('density')
        plt.xlabel('x')
        # plt.axis([0, 1, 0, 1.1])

        plt.figure(3)
        plt.plot(x_analytic,p, color='g')
        plt.ylabel('velocity')
        plt.xlabel('x')

    plt.show()
    plt.close()

def get_flux(_e, _values, _dx, _dt):
    f = np.multiply(np.matrix([_e]).T,(_values[:,1:]-_values[:,:-1])*_dt/_dx)
    return f


def compute_eigenvalues(_U,_gamma):

    rho=_U[0]
    v = _U[1]
    p = _U[2]
    A = np.array([[v,rho,0.],[0.,v,1.0/rho],[0.,_gamma*p,v]])

    e,R = np.linalg.eig(A)
    L = np.linalg.inv(R)
    return {'e':e,'L':L,'R':R}


def convert_to_diag(_values,_L):
    diags = np.matmul(_L,_values)
    return diags


def convert_to_prim(_diags,_R):
    prim = np.matmul(_R,_diags)
    return prim

def analytic_solution(gamma, npts):
    positions, regions, values = sod.solve(left_state=(8./gamma, 8, 0), right_state=(1, 1, 0.),
                                           geometry=(-0.5, 0.5, 0.), t=0.2, gamma=gamma, npts=npts)
    return values
#abs(max([3, 7, -10], key=abs)) EOC

main()

