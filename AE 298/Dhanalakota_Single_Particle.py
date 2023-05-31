import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append('../shocktubecalc')
from shocktubecalc import sod

def main(P, U, D, M):

    # Given constants and variables
    C_d = 0.42 # Drag coefficient
    gamma = 1.4
    points = 500 # Number of plot points in graph
    first_x = -0.25 # First x-value in graph
    positions, regions, values = sod.solve(left_state=(1, 1, 0), right_state=(0.125, 0.1, 0.),
                                           geometry=(first_x, 1., 0), t=0.2, gamma=gamma, npts=points)

    # Defining necessary variables
    F_d = np.zeros((points,)) # Creating an array for the drag force values
    F_pg = np.zeros((points,)) # Creating an array for the pressure gradient force values
    F_am = np.zeros((points,)) # Creating an array for the added mass force values
    A = np.pi * D**2/4 # Projected Area
    u_p = np.zeros((points,)) # Particle velocity array
    u_p[0] = U # Current particle velocity
    V_p = 4/3 * np.pi * (D/2)**3 # Particle volume
    C_m = 1/2 # Standard coefficient for added mass force
    dt = 0.01 # Time step
    
    for i in range(0, points - 1):
        
        # Drag force calculation
        F_d[i] = 1/2 * C_d * values['rho'][i] * (u_p[i] - values['u'][i])**2 * A
        
        # Pressure gradient force calculation
        P_gradient = (values['p'][i+1] - values['p'][i-1]) / ((values['x'][points - 1] - values['x'][0])/points*2)
        F_pg[i] = - P_gradient * V_p

        # Added mass force
        F_am[i] = V_p * C_m * (- P_gradient - (values['rho'][i+1] * u_p[i] - values['rho'][i-1] * u_p[i-1])/dt)

        u_p[i+1] = u_p[i] + (F_d[i] + F_pg[i] + F_am[i])/M * dt # V_final = V_initial + (force / mass) * time (1)

    F_pg[0] = F_pg[1]
    F_am[0] = F_am[1]

    for i in range(0, points - 1):
        if F_pg[i] > 2e-11:
            # print('Previous F_pg: ', F_pg[i-1])
            # print('F_pg: ', F_pg[i])
            # print('Next F_pg: ', F_pg[i+1])
            # print('F_am: ', F_am[i])
            # print('x: ', values['x'][i])
            # print()
            F_pg[i] = F_pg[i-1]
        
        if F_am[i] > 2.1e-12:
            # print('Previous F_am: ', F_am[i-1])
            # print('F_am: ', F_am[i])
            # print('Next F_am: ', F_am[i+1])
            # print('F_am: ', F_am[i])
            # print('x: ', values['x'][i])
            # print()
            F_am[i] = F_am[i-1]

    # plot values
    f, axarr = plt.subplots(5, 1, sharex='col')

    axarr[0].plot(values['x'], F_d, linewidth=1.5)
    axarr[0].set_ylabel(r'Drag force')
    axarr[0].set_xlabel(r'$x$')

    axarr[1].plot(values['x'], F_pg, linewidth=1.5)
    axarr[1].set_ylabel(r'Pressure gradient force')
    axarr[1].set_xlabel(r'$x$')

    axarr[2].plot(values['x'], F_am, linewidth=1.5)
    axarr[2].set_ylabel(r'Added mass force')
    axarr[2].set_xlabel(r'$x$')

    axarr[3].plot(values['x'], u_p, linewidth=1.5)
    axarr[3].set_ylabel(r'Particle velocity')
    axarr[3].set_xlabel(r'$x$')

    axarr[4].plot(values['x'], values['u'], linewidth=1.5)
    axarr[4].set_ylabel(r'Gas velocity')
    axarr[4].set_xlabel(r'$x$')

    plt.show()

if __name__ == '__main__':

    # Given values
    u_p0 = 0
    d1 = 1e-6 # m
    d2 = 10e-6 # m
    d3 = 115e-6 # m

    rho_p = 2520 # kg/m^3

    x0 = 0 # Initial Particle Position

    # Particle Mass = Volume * Density
    m_p1 = 1/6 * np.pi * d1**3 * rho_p
    m_p2 = 1/6 * np.pi * d2**3 * rho_p
    m_p3 = 1/6 * np.pi * d3**3 * rho_p

    N = 500

    for i in range (N):
        
        pass
    
    main(x0, u_p0, d3, m_p3)
    # main(x0, u_p0, d2, m_p2)
    # main(x0, u_p0, d3, m_p3)