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
#     calculate Delta_max = max - position, Delta_min = position - min (if <0, particle is outside)
#     if D_max > 0 or D_min > 0, set == 0 (particles are inside the domain, no change will take effect)
#     while there is D_max != 0 or D_min != 0 :
#         if particle settled (z<zmin), move it to end of property array and change its flag
#         update active particle counter
#         position = position - 2* |D_max| + 2 *|D_min|
#         where D_max != 0 or D_min != 0 , flip the velocity sign
# 
# visualize    

# In[44]:


import os
import time as timer

import cupy as cp
from numba import njit
from rich.progress import track
from scipy.stats import gamma, pareto


# In[45]:


"""CONSTANTS IN SI"""
RH, Tamb, TD = (
    40,
    10,
    1,
)  # Relative humidity [%], Ambient temperature [oC], Turbulent diffusivity (empirical range, varies from 0 to 3)
SAMPL = 1  # [s], interval of seconds for printing output
VEMIT = 0.3  # [m/s], particle initial velocity magnitude
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


# In[46]:


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


# In[58]:


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
Tf = Tw + 1/3*(Tamb-Tw)

def evap(D, V, dt):
    """RETURNS UPDATED PARTICLE DIAMETER AFTER EVAPORATION
    Input: D diameter, float, from row 6 of particle property array, X[6,:,0]
           V velocity, float, [3:6 row of particle property array
           dt timestep
    Output:
           D diameter, float, from row 6 of particle property array, X[6,:,1]
    """
    Tf = Tw + 1/3 *(Tamb - Tw) # oC
    alpha_evap = 4 * MW * DV * VP / (RHOL * R * (Tf + 273.15))
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


# In[59]:


### BOUNDARY ###

## define limits ##
box_min = cp.array([0,0,0])
box_max = cp.array([3,3,3])


# In[60]:


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

        ## initial particle velocity
        X[3:6, :, 0] = init_vel(
            VEMIT, N
        )  # initialize with some randomly oriented velocity

        # diameter sampling
        X[6, :, 0] = D
        # time elapsed in room
        X[7, :, 0] = 0
        # flag 1:active 0:inactive(settled)
        X[8, :, :] = 1
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
        * (cp.power(X[6, :, 0], 3) - cp.power(X[6, :, 1], 3))/dt 
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


# In[61]:


def settlingflag(X,activep,Np,mask):
    '''For active particles of the property array that have settled, this function changes their flag from 1 to 0.
    First creates trail of zeros for inactive particles. Then appends it to a boolean array on the active particles (True:settled, False:active). 
    
    INPUT: 
    X, array: particle property array
    activep, int: counter of active particles
    Np, int: total number of particles
    mask, boolean array: TRUE: settled (z<zmin, FALSE: active)

    OUTPUT:
    X, array: particle property array with updated settling flags          
    '''
    # check if particle settled and change flag
    trail=[False] * (Np-activep)      
    maskfull=mask.tolist() + trail
    X[8,maskfull,1]=0 #change flag for inactive particles

def settled_ids(mask):
    '''Returns indices of inactive particles from mask
    INPUT: mask, boolean array: TRUE: settled (z<zmin, FALSE: active)
    OUTPUT: settled_id, tuple: indices of inactive particles
    '''
    settled_id=cp.where(mask)[0].tolist() #[0] is for indices, [1] for values. Return index of inactive particles from previous mask
    return settled_id

def movetoend(X,activep,settled_id):
    '''Re-orders particles within the property array based on their flag. First creates list of all active indices and deletes the inactive ones. 
    Then appends inactive indices at trail. Re-arranges order of columns based on list.
    INPUT:
    X, array: particle property array
    activep, int: counter of active particles
    settled_id: tuple with indices of settled particles
    
    OUTPUT:
    X, array: particle property array with columns corresponding to inactive particles moved to end of array
    '''
    all_id = list(range(activep))  # all indices of activep (aka particles that were active at start of this timestep)
    active_first = [i for j, i in enumerate(all_id) if j not in settled_id]  # delete inactive particles (keep only active indices)
    active_first.extend(settled_id)  # add the inactive ones at trail (indices list: active, inactive)
    X[:,:activep,1]=X[:,active_first,1]

def delta_arrays(X,activep):
    '''Generates arrays that calculate difference of all positions from maxima and minima.
    INPUT:
    X, array: particle property array
    activep, int: counter of active particles
    
    OUTPUT:
    delta_max, delta_min: 3x(activep) arrays that contain (box_max - pos) and (pos - box_min)
    '''
    delta_max = box_max[:,None] - X[:3,:activep,1] # if <0, particle is outside
    delta_min = X[:3,:activep,1] - box_min[:,None] # if <0, particle is outside

    delta_max[delta_max>0] = 0 # if the delta_max > 0, particle is inside , no correction applied
    delta_min[delta_min>0] = 0 # if the delta_min > 0, particle is inside , no correction applied
    return delta_max, delta_min

def settling(X,activep,Np,mask):
    '''If particle has settled based on boolean mask, change its flag and move it to end of particle property array
    INPUT:
    X, array: particle property array
    activep, int: counter of active particles
    Np,int: total number of particles
    mask, boolean array: TRUE: settled (z<zmin, FALSE: active)
    '''
    #mask=X[2,:activep,1]<0 #boolean, True: particle became inactive, False: particle active
    
    if mask.any():
        settlingflag(X,activep,Np,mask) #change flag of settled particles
        settled_id=settled_ids(mask) #return indices of settled particles

        ### Move inactive particles to end of array - activep here refers to particles considered 'active' during previous timestep
        movetoend(X,activep,settled_id)
        
def updatecounters(activep,Sp,mask,Np):
    '''Updates the counters of settled and active particles based on TRUE values of boolean array that checks (z<zmin).
    INPUT:
    activep, int: counter of active particles
    Sp, int: counter of settled particles
    mask, boolean array: TRUE: settled (z<zmin, FALSE: active)
    Np, int: total number of particles
    '''
    settled_id=settled_ids(mask) #return indices of settled particles
    Sp += len(settled_id)
    activep = Np - Sp
    return activep,Sp


# In[65]:


def main(D, restart=None):
    """TIME LOOP
    FOR EVERY STEP, FIND IF PARTICLES HAVE SETTLED AND MOVE THEM TO END OF X PARTICLE ARRAY.
    UPDATE VELOCITY AND POSITION OF PARTICLES AT BEGINNING OF ARRAY
    Input: diameter, float in [m]
           restart, path to file that restarts simulation from a certain particle array
    """
    Np = int(1000)
    Sp = int(0)  # initialize exit counter (counts particles that settle)
    t, X = init_particles(D, restart, Np)

    for i in track(range(t, T + 1)):
        print(i,"timestep", file = sourceFile)
        activep = Np - Sp # python indexing from 0, range function
        print('active_parts',activep,file = sourceFile)
        upd_particles(X[:, : activep, :])  # apply only on non-settled particles
        print("new positions calculated", file=sourceFile)
        print(X[:,:,1], file = sourceFile)   
        
        # Boundary condition - only on active particles
        # generate delta_arrays
        
        delta_max,delta_min=delta_arrays(X,activep)        
        
        while (delta_max.any()) or (delta_min.any()):
            
            mask=X[2,:activep,1]<0
            settling(X,activep,Np,mask) #changes flag and moves to end
            print('settled particles: changed flag and moved to end', file=sourceFile)
            print(X[:,:,1], file=sourceFile)
            ### Update settled particle counter
            activep,Sp=updatecounters(activep,Sp,mask,Np)
            print('after counter update', file=sourceFile)
            print('active',activep,'settled',Sp, file=sourceFile)  
            
            if Sp == Np:  # if all particles have settled, break the loop
                print("all particles settled")
                PATH = os.path.join(os.getcwd(), "data_gen")
                filename = "simple_bin_D{}_RH_{}_T_{}_TD{}_dt{}_Np{}k_i_{}settled.npz".format(
                    D, RH, Tamb, TD, dt, int(Np / 1000), int(i * dt)
                )
                cp.savez(
                    os.path.join(PATH, filename), S=Sp, X=X[:, :, :], t=i, dt=dt, T=T, TD=TD
                )
                return X
            
            #re-generate delta arrays after having moved particles to end    
            delta_max,delta_min=delta_arrays(X,activep)
            print('delta_max',delta_max,file=sourceFile)
            print('delta_min',delta_min,file=sourceFile)
            print('condition check',(delta_max.any()) or (delta_min.any()),file=sourceFile)
            
            print('in the loop of reflection. Old pos:', file=sourceFile)        
            print(X[:3,:,1], file=sourceFile) 
            
            X[:3,:activep,1] += - 2*abs(delta_max) + 2*abs(delta_min)
            
            print('in the loop of reflection. New pos:', file=sourceFile)        
            print(X[:3,:,1], file=sourceFile) 
            
            # update velocity
            # make mask. wherever non-zero element in delta arrays is TRUE, factor is (-1), where FALSE, factor is 1 (no flip of sign)
            mask_vel = cp.where((delta_max!=0) | (delta_min!=0), -1, 1)
            X[3:6,:activep,1] = mask_vel*X[3:6,:activep,1]
            
            print('in the loop of reflection. New vel:', file=sourceFile)        
            print(X[3:6,:,1], file=sourceFile)
            
            #update the delta arrays, check condition again
            delta_max,delta_min=delta_arrays(X,activep)

            print('cond',(delta_max.any()) or (delta_min.any()),file=sourceFile)


        if i % INTSTEPS == 0:  # T intsteps
            PATH = os.path.join(os.getcwd(), "data_gen")
            filename = "simple_bin_D{}_RH_{}_T_{}_TD{}_dt{}_Np{}k_i_{}.npz".format(
                D, RH, Tamb, TD, dt, int(Np / 1000), int(i * dt)
            )  # _{}_of20_RH_{}_T_{}_Ct_{}_Np_{}k_i_{}.npz".format(kk,RH,Tamb,Ct,int(Np/1000),i)
            cp.savez(
                os.path.join(PATH, filename), S=Sp, X=X[:, :activep, :], t=i, dt=dt, T=T, TD=TD
            )

        X[:, :, 0] = X[:, :, 1]  # Xold = Xnew

    print("Done, Sp is", Sp)
    return X


# In[66]:


time = 10*60 # [s] simulated time, e.g. 10 minutes * 60 seconds/min = 600 seconds
T = int(time / dt)  # steps needed


# In[67]:
# get_ipython().run_cell_magic('capture', 'cap --no-stderr', 'sourceFile = open("outputsimple.txt", "w")\nX = main(100e-6, None)  # diameter in [m], here: 150 um * 10-6 m/um\nsourceFile.close()\n# with open(\'output.txt\', \'w\') as f:\n#    f.write(cap.stdout)')

