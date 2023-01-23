#!/usr/bin/env python
# coding: utf-8

# BOUNDARY CONDITION
# 
# PARTICLES ARE IN A 3D CUBE. IF A PARTICLE TOUCHES THE BOTTOM PLANE, IT SETTLES (EXITS THE SIMULATION).
# 
# IF A PARTICLE TOUCHES ANY OF THE OTHER 5 PLANES OF THE CUBE, IT BOUNCES BACK INSIDE THE CUBE.
# 
# while t<time:
# 
#     initialize particle property array X (position, velocity, diameter) [init_particles] 
#     particle array X [8,N,2] 
#         8:properties
#         N:particles
#         2:refers to old and new timestep
#     
#     update diameter, velocity and position  [upd_particles]
# 
#     is new position outside the cube?
#         no: then continue the loop above
#         yes: boundary condition needed. is position out by >1 coordinate?
#                 yes: find which wall was hit first by the particle
#                 no: is particle settled? (is z<zmin?)
#                         yes: move particle to the end of the array, stop updating it
#                         no: apply elastic collision with the wall (find intercept of particle with the wall it went through,
#                             invert corresponding velocity component)
# visualize    

# In[103]:


import os

import cupy as cp
from rich.progress import track
from scipy.stats import gamma, pareto
from numba import njit


# In[104]:


"""CONSTANTS IN SI"""
RH, Tamb, TD = 40, 10, 1 # Relative humidity [%], Ambient temperature [oC], Turbulent diffusivity (empirical range, varies from 0 to 3)
SAMPL = 10  # [s], interval of seconds for printing output
VEMIT = 2  # [m/s], particle initial velocity magnitude
RHOL = 1000  # [kg/m3], density of liquid
RHOG = 1  # [kg/m3], density of air
MUG = 1.8e-5  # [kg/m/s], air viscosity
G = 9.8  # [m/s2], magnitude of acceleration of gravity
K = 1.38e-23  # [J/K], Boltzmann constant
DMIN = 1e-8  # [m], smallest particle diameter size (evaporation threshold)
MW = 0.01802  # [kg/mol], molecular weight of water
R = 8.3145  # [J/mol/K]  or [m3*Pa/K/mol], universal gas constant 
dt = 0.1  # [s], timestep 

INTSTEPS = int(SAMPL / dt)  # calculated number of steps to save output

### SOME USEFUL FUNCTIONS ###


# In[105]:


def init_vel(vel_magn, nd):
    """Generate a random initial velocity orientation for a given velocity magnitude
    
    Input: vel_magn: m/s, velocity magnitude, float
           nd: number of droplets, int
           
    Output: array of 
    """
    # initial particle velocity with random orientation
    Y = cp.random.randn(3, nd)
    norm = cp.sqrt(cp.sum((Y**2), axis=0))  # L2, not Frobenius
    Y /= norm
    return Y * vel_magn


# In[106]:


# drag calculation
def drag(Reyn, Rebl):  # Reynolds slip velocity, Reynolds blowing velocity
    """Calculate particle drag"""
    a = 0.09 + 0.077 * cp.exp(-0.4 * Reyn)
    b = 0.4 + 0.77 * cp.exp(-0.04 * Reyn)
    nom = 1 + 0.0545 * Reyn + 0.1 * (Reyn) ** (1 / 2) * (1 - 0.03 * Reyn)
    denom = 1 + (a * abs((Rebl)) ** b)

    return nom / denom


# gravity-buoyancy
""" Calculate gravity vector (acting only on z-axis)"""
grav = cp.array([0, 0, G * ((RHOG / RHOL) - 1)])
grav = grav.reshape(3, 1)


def fun(D, F):
    """Calculates two quantities needed for velocity update"""
    tau = (RHOL * D**2) / (18 * MUG)  # Stokes regime!

    m = (cp.pi * RHOL * D**3) / 6  # particle mass
    S = K * (Tamb + 273.15) / (m * tau)
    Sdif = cp.sqrt(2 * S)
    taust = tau / F
    return taust, Sdif


# evaporation
""" CONSTANTS FOR PARTICLE EVAPORATION 
Input values of RH: Relative humidity and Tamb: ambient temperature
are then used to find Tw, DV, VP that are needed for function evap"""
a0 = 5.1055
a1 = 0.4295
b0 = -0.04703
b1 = -0.005951
c0 = -4.005e-5
c1 = 1.66e-5

Tw = Tamb - ((a0 + a1 * Tamb) + (b0 + b1 * Tamb) * RH + (c0 + c1 * Tamb) * RH**2)
DV = 21.2e-6 * (1 + 0.0071 * Tw)
VP = 67 * (Tamb - Tw)

def evap(D, V, dt):
    """RETURNS UPDATED PARTICLE DIAMETER AFTER EVAPORATION
    Input: D diameter, float, from row 6 of particle property array, X[6,:,0]
           V velocity, float, [3:6 row of particle property array
           dt timestep
    Output:
           D diameter, float, from row 6 of particle property array, X[6,:,1]
    """
    alpha_evap = 4 * MW * DV * VP / (RHOL * R * (Tw + 273.15))
    beta_evap = 0.276 * (RHOG / (MUG * (DV) ** 2)) ** (1 / 6)
    # evaporation rate, D diameter, V velocity
    sqrttest = D**2 - 2 * dt * alpha_evap * (1 + beta_evap * cp.sqrt(D * V))
    Dnew = cp.empty_like(D)
    Dnew = cp.where(
        sqrttest < 0,
        DMIN,
        cp.sqrt(D**2 - 2 * dt * alpha_evap * (1 + beta_evap * cp.sqrt(D * V))),
    )

    # extracheck
    Dnew = cp.where(
        Dnew < 0, DMIN, Dnew
    )  # if negative value occurs (timestep dt is too large), set D to Dthres
    return Dnew


# In[107]:


### BOUNDARY ###

## define limits ##
DIM = 3
XLIM = cp.array([0, 1]) * DIM
YLIM = cp.array([0, 1]) * DIM
ZLIM = cp.array([0, 1]) * DIM


# In[108]:


def init_particles(D, restart, N=500000):
    """Initialize array X[:,:,0] with particle properties.
    X [8,N,2]
    8 : particle position (3 rows), particle velocity (3 rows), diameter (1 row), time elapsed in simulation per particle (1 row).
    N : one column per particle
    2 : two states of X, for old and new timestep (during integration we need quantities from previous timestep)
    """
    if restart:
        # checkpoint = np.load(restart)
        # checkpoint = cp.asarray(checkpoint)
        checkpoint = cp.load(restart)
        X = checkpoint["X"]
        t = int(checkpoint["t"])
        N = X.shape[1]
    else:
        t = 0
        X = cp.empty((9, N, 2))
        X[:3, :, 0] = 1.5
        #X[0,:,0]=2.65
        #X[1,:,0]=2.9
        #X[2,:,0]=2.95
        # initial particle velocity
        X[3:6, :, 0] = init_vel(
            VEMIT, N
        )  # initialize with some randomly oriented velocity
        # diameter sampling
        X[6, :, 0] = D
        # time elapsed in room
        X[7, :, 0] = 0
        # flag 0:active 1:inactive(settled)
        X[8, :, 0] = 0
    return t, X


# particle loop, n particles (only for particles that haven't settled)
def upd_particles(X):
    """Updates particle array X[:,:,1] based on X[:,:,0]
    Input:
        X, array: particle coordinates and property array
    Output:
        X, array: particle coordinates and property array with updated X[:,:,1]

    """

    n = X.shape[1]  # n no. of particles, X contains particle states (X, V, D, t)

    ### EVAPORATION ###
    # X[6,:,1] = X[6,:,0] #deactivate evaporation
    X[6, :, 1] = cp.where(
        X[6, :, 0] > DMIN,
        evap(X[6, :, 0], cp.linalg.norm(X[3:6, :, 0]), dt),
        X[6, :, 0],
    )

    ### DRAG ###

    # blowing velocity
    Vbl = (
        RHOL
        * (cp.power(X[6, :, 0], 3) - cp.power(X[6, :, 1], 3))
        / (6 * RHOG * cp.power(X[6, :, 1], 2))
    )

    # Reynolds numbers with slip and blowing velocity
    Re_s = RHOG * X[6, :, 0] * (cp.linalg.norm(X[3:6, :, 0], axis=0)) / MUG
    Re_bl = RHOG * X[6, :, 0] * Vbl / MUG

    f = drag(Re_s, Re_bl)

    ### TIME INTEGRATION ###
    # characteristic time and diffusion coeff
    taus, Sd = fun(X[6, :, 0], f)

    # velocity update
    X[3:6, :, 1] = (
        (cp.exp(-dt / taus)) * X[3:6, :, 0]
        + taus * (1 - cp.exp(-dt / taus)) * (grav.reshape(3, 1))
        + (Sd + TD)
        * (cp.sqrt((taus / 2) * (-cp.expm1(-2 * dt / taus))))
        * cp.random.randn(3, n)
    )

    # position update
    X[:3, :, 1] = X[:3, :, 0] + dt * X[3:6, :, 0]
    
    # time elapsed per particle
    X[7, :, 1] = X[7, :, 0] + dt


# In[109]:


#@njit
def boundary(activep,X,Sp): #active p: int, active particles, ind_set: list with indices of settled particles, Sp: settled particles
    """ LOOPS OVER ACTIVE PARTICLES AND APPLIES BOUNDARY CONDITION ON THEIR NEWLY-CALCULATED POSITION
    Input:
    activep: int, number of active particles
    X: particle array to be updated
    Sp: int, number of settled particles
    
    Output:
    X: updated particle array
    ind_set: list, indices of settled particles
    Sp: int, counter of settled particles
    
    """
    

    for ii in range(activep): #loop over active particles, starts at 0
            ## CONDITIONS ##
            cxmin=X[0,ii,1]-XLIM[0]<=0  #xpos < xmin
            cxmax=XLIM[1]-X[0,ii,1]<=0  #xpos > xmax 
            cymin=X[1,ii,1]-YLIM[0]<=0  #ypos < ymin
            cymax=YLIM[1]-X[1,ii,1]<=0  #ypos > ymax
            czmin=X[2,ii,1]-ZLIM[0]<=0  #zpos < zmin -- LEAVES SIMULATION DOMAIN
            czmax=ZLIM[1]-X[2,ii,1]<=0  #zpos > zmax
            
            flags = [cxmin, cxmax, cymin, cymax, czmin, czmax] 
            #p0=cp.array([[0,0,0], [ 3,0,0], [ 0,0,0], [ 0,3,0], [ 0,0,0], [ 0,0,3]])
            
            ind_set=[] #list to put indices of settled particles for this timestep
            
            # XOR condition??
            # values = [32, 16, 8, 4, 2, 1]
            # boundary_val = sum ([value for (flag, value) in list(zip(flags, values)) if flag])

            number_of_sides = sum([1 for flag in flags if flag ])  # counts number of sides by which particle is off
            
            if number_of_sides == 0: # particle is inside box
                continue            
            
            if number_of_sides == 1: # it went out of one side
                l=X[:3,ii,1]-X[:3,ii,0] #line connecting particle's old and new positions
                l0=X[:3,ii,0] #a point on that line (e.g. the old position)
                if cxmin: 
                    #find intercept
                    p0=cp.array([0,0,0]) #point on -x plane
                    nrm=cp.array([-1,0,0]) #normal (pointing outward) of +x plane                  
                    d=(cp.subtract(p0,l0)*nrm).sum(axis=0)/(l*nrm).sum(axis=0)
                    p=l0+l*d
                    #check p
                    X[:3,ii,1]=p
                    X[3:6,ii,1]=cp.multiply(X[3:6,ii,0],cp.array([-1,1,1])) #flip x-vel
                    
                elif cxmax:
                    p0=cp.array([3,0,0]) #point on +x plane
                    nrm=cp.array([1,0,0])
                    d=(cp.subtract(p0,l0)*nrm).sum(axis=0)/(l*nrm).sum(axis=0)
                    p=l0+l*d
                    X[:3,ii,1]=p
                    X[3:6,ii,1]=cp.multiply(X[3:6,ii,0],cp.array([-1,1,1])) #flip x-vel

                elif cymin: 
                    p0=cp.array([0,0,0]) #point on -y plane
                    nrm=cp.array([0,-1,0])
                    d=(cp.subtract(p0,l0)*nrm).sum(axis=0)/(l*nrm).sum(axis=0)
                    p=l0+l*d
                    X[:3,ii,1]=p
                    X[3:6,ii,1]=cp.multiply(X[3:6,ii,0],cp.array([1,-1,1])) #flip y-vel
                    
                elif cymax: 
                    p0=cp.array([0,3,0]) #point on +y plane
                    nrm=cp.array([0,1,0])
                    d=(cp.subtract(p0,l0)*nrm).sum(axis=0)/(l*nrm).sum(axis=0)
                    p=l0+l*d
                    X[:3,ii,1]=p
                    X[3:6,ii,1]=cp.multiply(X[3:6,ii,0],cp.array([1,-1,1])) #flip y-vel
                elif czmin: #particle exits simulation
                    #roll?? cp.roll(,-1, axis=1)
                    ind_set.append(ii) #add particle index to settled particles
                    X[8,ii,1]=1 #make particle inactive
                    Sp += 1
                    
                elif czmax: 
                    p0=cp.array([0,0,3]) #point on +z plane
                    nrm=cp.array([0,0,1])
                    d=(cp.subtract(p0,l0)*nrm).sum(axis=0)/(l*nrm).sum(axis=0)
                    p=l0+l*d
                    X[:3,ii,1]=p 
                    X[3:6,ii,1]=cp.multiply(X[3:6,ii,0],cp.array([1,1,-1])) #flip z-vel

            
            if number_of_sides > 1:
              # it went out of two or three sides.
              # time it took for it to reach the wall based on previous timestep's pos and velocity
                txmin=(XLIM[0]-X[0,ii,0])/X[3,ii,0] #Xmin-X<0, if -Vx then resulting t is >0 
                txmax=(XLIM[1]-X[0,ii,0])/X[3,ii,0] 
                tymin=(YLIM[0]-X[1,ii,0])/X[4,ii,0]
                tymax=(YLIM[1]-X[1,ii,0])/X[4,ii,0]
                tzmin=(ZLIM[0]-X[2,ii,0])/X[5,ii,0]
                tzmax=(ZLIM[1]-X[2,ii,0])/X[5,ii,0]

                times=[txmin, txmax, tymin, tymax, tzmin, tzmax]
                value = min([i for i in times if i >= 0])
                
                l=X[:3,ii,1]-X[:3,ii,0] #line connecting particle's old and new positions
                l0=X[:3,ii,0] #a point on that line (e.g. the old position)
                
                if value == txmin:
                    #find intercept
                    p0=cp.array([0,0,0]) #point on -x plane
                    nrm=cp.array([-1,0,0]) #normal (pointing outward) of +x plane                  
                    d=(cp.subtract(p0,l0)*nrm).sum(axis=0)/(l*nrm).sum(axis=0)
                    p=l0+l*d
                    X[:3,ii,1]=p
                    X[3:6,ii,1]=cp.multiply(X[3:6,ii,0],cp.array([-1,1,1])) #flip x-vel
                elif value == txmax:
                    p0=cp.array([3,0,0]) #point on +x plane
                    nrm=cp.array([1,0,0])
                    d=(cp.subtract(p0,l0)*nrm).sum(axis=0)/(l*nrm).sum(axis=0)
                    p=l0+l*d
                    X[:3,ii,1]=p
                    X[3:6,ii,1]=cp.multiply(X[3:6,ii,0],cp.array([-1,1,1])) #flip x-vel
                elif value == tymin:
                    p0=cp.array([0,0,0]) #point on -y plane
                    nrm=cp.array([0,-1,0])
                    d=(cp.subtract(p0,l0)*nrm).sum(axis=0)/(l*nrm).sum(axis=0)
                    p=l0+l*d
                    X[:3,ii,1]=p
                    X[3:6,ii,1]=cp.multiply(X[3:6,ii,0],cp.array([1,-1,1])) #flip y-vel
                elif value == tymax:
                    p0=cp.array([0,3,0]) #point on +y plane
                    nrm=cp.array([0,1,0])
                    d=(cp.subtract(p0,l0)*nrm).sum(axis=0)/(l*nrm).sum(axis=0)
                    p=l0+l*d
                    X[:3,ii,1]=p
                    X[3:6,ii,1]=cp.multiply(X[3:6,ii,0],cp.array([1,-1,1])) #flip y-vel
                elif value == tzmin:
                    ind_set.append(ii) #add particle index to settled particles
                    X[8,ii,1]=1 #make particle inactive
                    Sp += 1
                elif value == tzmax:
                    p0=cp.array([0,0,3]) #point on +z plane
                    nrm=cp.array([0,0,1])
                    d=(cp.subtract(p0,l0)*nrm).sum(axis=0)/(l*nrm).sum(axis=0)
                    p=l0+l*d
                    X[:3,ii,1]=p 
                    X[3:6,ii,1]=cp.multiply(X[3:6,ii,0],cp.array([1,1,-1])) #flip z-vel

    return X,ind_set,Sp


# In[110]:


def main(D, restart=None):
    """TIME LOOP
    FOR EVERY STEP, FIND IF PARTICLES HAVE SETTLED AND MOVE THEM TO END OF X PARTICLE ARRAY.
    UPDATE VELOCITY AND POSITION OF PARTICLES AT BEGINNING OF ARRAY
    Input: diameter, float in [m]
           restart, path to file that restarts simulation from a certain particle array
    """
    Np = 1000
    Sp = 0  # initialize exit counter (counts particles that settle)
    t, X = init_particles(D, restart, Np)
    
    for i in track(range(t, T + 1)):

        upd_particles(X[:, : (Np - Sp), :])  # apply only on non-settled particles
        
        activep=int(Np-Sp) #python indexing from 0, range function    
        
        #@jit(parallel=True)
        X,ind_set,Sp=boundary(activep,X,Sp)
        ### Move inactive particles to end of array, Sp = counter of inactive particles   
        ### OR JUST DELETE np.delete(X,ind_set,axis=1)
        all_id=list(range(X.shape[1])) #all indices of X
        active_first = [i for j, i in enumerate(all_id) if j not in ind_set] #delete inactive particles (keep only active indices)
        active_first.extend(ind_set) #add the inactive ones at trail (indices list: active, inactive)

        X[:,:,:] = X[:,active_first,:]
        
        if Sp == Np:  # if all particles have settled, break the loop
            print("all particles settled")
            PATH = os.getcwd()
            filename = "bin_D{}_RH_{}_T_{}_TD{}_dt{}_Np{}k_i_{}settled.npz".format(
                D, RH, Tamb, TD, dt, int(Np / 1000), int(i * dt)
            )
            cp.savez(
                os.path.join(PATH, filename), S=Sp, X=X[:, :, :], t=i, dt=dt, T=T, TD=TD
            )
            return X

        if i % INTSTEPS == 0:  # T intsteps
            PATH = os.getcwd()
            filename = "bin_D{}_RH_{}_T_{}_TD{}_dt{}_Np{}k_i_{}.npz".format(
                D, RH, Tamb, TD, dt, int(Np / 1000), int(i * dt)
            )  # _{}_of20_RH_{}_T_{}_Ct_{}_Np_{}k_i_{}.npz".format(kk,RH,Tamb,Ct,int(Np/1000),i)
            cp.savez(
                os.path.join(PATH, filename), S=Sp, X=X[:, :, :], t=i, dt=dt, T=T, TD=TD
            )

        X[:, :, 0] = X[:, :, 1]  # Xold = Xnew

    print("Done, Sp is", Sp)
    return X


# In[ ]:


time = 10 * 60 # [s] simulated time, e.g. 10 minutes * 60 seconds/min = 600 seconds
T = int(time / dt)  # steps needed    

X = main(150e-6, None) # diameter in [m], here: 150 um * 10-6 m/um , 

