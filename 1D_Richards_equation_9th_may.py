
"""

@author: Saurabh Kumar 
Code based on Dogan and Mortz 2002a,b 
# Dogan, A., & Motz, L. H. (2005). Saturated-unsaturated 3D groundwater model. I: Development. 
# Dogan, A., & Motz, L. H. (2005). Saturated-unsaturated 3D groundwater model. II: Verification and application

Trial version of uniform grid 


"""
import numpy as np
import matplotlib.pyplot as plt
from VanGenuchten_function_files import f1,f2, Van_Genuchten_moisure,Van_Genuchten_K,Van_Genuchten_specific_moisure_capacity
import pandas as pd
import scipy.io




df = pd.read_excel('meteo_data.xlsx')  # For .xlsx files
rain = np.array(df["rainfall_mm_per_d"])  # mm/day
PE = np.array(df["PE_mm_per_d"])   # mm/day
rain = rain/1000 # m/day
Ep = PE/1000 # m/day potential evaporation


zmin = 0
zmax = 3  # in m
dz = 0.1 # in m
z = np.arange(zmin + dz/2,zmax+dz/2,dz)

# # spacing is non uniform 
# z_1 = 0
# z_2 = 1
# z_3 = 3
# dz_2 = 0.1818
# dz_i = 0.0022

# # z_2 = np.arange(3,1,-0.1818)
# # z_2_nodes = len(z_2)
# # z_1_nodes = z_nodes - z_2_nodes 

# grid_data = scipy.io.loadmat("grid_spacing.mat")
# dz_all = grid_data.get("dz_all")[0]
# z = grid_data.get("z_req")[0]



tmin = 1; # starting time in day
tmax = len(rain); # max time in day
dt = 1; #in day
# dt_min = 0.001 
# dt_max = 1
# inc_factor = 0.2  # increasing factor in the code
# dec_fact = 0.1 # decreasing factor in the code



t = np.arange(tmin,tmax+dt,dt) 

Ksat = 3.220   # m/day
thetas = 0.3769 
thetar = 0.0515
alpha = 3.321    # air entry presure in m
N = 2.503 # van Genucten Paramter 
Ss = 0.00 #  specific storage 
n_eta =  -0.8653 # Relative permeability exponen

q = - rain  # flux in m/day at the top boundary 
tnodes = len(rain)
znodes = len(z)
# preallocate h and theta 
h = np.empty([znodes,tnodes])
theta=np.empty([znodes,tnodes])
# storage water level within the system
water_storage = np.empty([tnodes])
Ea = np.empty([znodes,tnodes])


z_wt = 8 # m. the water table is situated at depth greater than z_max


rho_pu = 600.6 # storage capacity at root water uptake, mm
rho_w = 174.0  # water_level at wiltng point,  mm

# root water uptake function 

psi_a = -0.05 # critical pressure heads associated with anaerobiosis,
psi_d = -4 # critical pressure heads associated with soilwater-limited evapotranspiration
psi_w = -150 # # critical pressure head associated with plant wilting

Lr = 1 # m # depth of root zone
a = 2 # m # measuremnet of how fast root zone declines with depth
depth = zmax- z


# calculation of 1st part of root water uptake term by Feddes 
f1_vals = f1(depth,a,Lr,zmax)

# initial presure condition 
h[:,0] = depth - z_wt 
# initial moisture condition
theta[:,0] = Van_Genuchten_moisure( h[:,0] , thetar, thetas , alpha, N)


# defining the maximum iteration condition and the convergence criteria
threshold = .001 # minimum convergence between sucessive iteration
MaxIterations = 50 # each time step have max number of iteration 
Error = np.zeros([tnodes,MaxIterations]) # preallocation of error calculated at each iteration and time step


# n refers to time step and m refers to itertion levels 

for n in range(len(rain)-1):
    h[:,n+1] = h[:,n] # assumption of new time step
    theta [:,n+1] = theta[:,n]
    print(n) # time step no
    
    for m in range(MaxIterations):
        # calculation of A and RHS for linear algebgra formation AX = RHS
        RHS = np.empty([znodes])
        A = np.zeros([znodes,znodes])
        K = Van_Genuchten_K(h[:,n+1],Ksat, thetar, thetas , alpha, N, n_eta) # hydraulic condutivity
        C = Van_Genuchten_specific_moisure_capacity(h[:,n+1], thetar, thetas, alpha, N ) # specific moisture capacity
        
        f2_vals = f2(h[:,n+1], psi_a, psi_d, psi_w)  # 2nd term of Feddes root water uptake model 
        RWU = f1_vals*f2_vals*Ep[n] # root water uptakte at each time step
        Ea[:,n] = RWU # storing the RWU as actual Evapotranspiration  # to correct it

        # calculation at each nodes          
        for i in range(0,znodes): # i refers to the node number, z is taken upwards
            
            h2  = h[i,n+1]
            K2  = K[i]
            theta_current = theta[i,n+1]
            theta_old = theta[i,n]
            # calculation of s
            s =(theta_current-theta_old)/dt          
          # calculation of p1 
            p1 = C[i]/dt
            # calculation of p2
            p2 = 0; # Sw*Ss/dt; since Sw=0
            
            if i ==0 :  # free drainage at the bottom layer
                K3Z = K[i+1]

                Kszplus = ( K2 + K3Z )/ 2
                CNZplusmean = -Kszplus/ 2
                g = -CNZplusmean/dz 
                d = -(  g + p1 + p2)
                A[i,i] = d      # coeffiecient of H(k)
                A[i,i+1] = g  # coeffiecient of H(k+1)
                RHS[i] = s - p1 *( h[i,n+1] + z[i] ) - p2*( h[i,n] + z[i]) - K2/dz - RWU[i]/dz  # RHS at i

           
            elif i < znodes-1: # central zone
                K1Z = K[i-1] 
                Kszminus = ( K2 + K1Z ) / 2 
                CNZminusmean = -Kszminus/ 2
                c = -CNZminusmean/dz;             
                # calculation of p1 
    #            C = Van_Genuchten_and_Nelson_C( h2 , thetar, thetas , hs , ho , N, Ss )
                p1 = C[i]/dt           
                # calculation of g
                K3Z = K[i+1]
                Kszplus = ( K2+ K3Z )/2 
                CNZplusmean = -Kszplus/dz
                g = -CNZplusmean/dz  
                # calculation of d          
                d = -(  c + g + p1 + p2)
                # forming the A matrix and RHS vector
                A[i,i-1] = c   # coeffiecient of H(k-1)
                A[i,i] = d      # coeffiecient of H(k)
                A[i,i+1] = g  # coeffiecient of H(k+1)
                RHS[i] = s - p1 *( h[i,n+1] + z[i] ) - p2*( h[i,n] + z[i]) - RWU[i]/dz # RHS at i
            elif i == znodes-1:   # surface boundary condition, only flux boundary condition
                K1Z = K[i-1] 
                Kszminus = ( K2 + K1Z )  / 2 
                CNZminusmean = -Kszminus/2
                c = -CNZminusmean/dz;             
                # calculation of p1 
    #            C = Van_Genuchten_and_Nelson_C( h2 , thetar, thetas , hs , ho , N, Ss )
                p1 = C[i]/dt 
                d = -(  c +  p1 + p2)
                A[i,i-1] = c   # coeffiecient of H(k-1)
                A[i,i] = d      # coeffiecient of H(k)
                RHS[i] = s - p1 *( h[i,n+1] + z[i] ) - p2*( h[i,n] + z[i]) + q[n+1]/dz - RWU[i]/dz
           
            
        X_input = h[:,n+1]+z
        # X = np.linalg.inv(A).dot(RHS)
        # X = np.linalg.solve(A,RHS)
        X = np.linalg.lstsq(A,RHS, rcond=None)[0] # due to formation of singular maxtrix, this one is being used 
        Error[n,m] =  np.sqrt(np.mean(pow(X-X_input,2))) # RMSE calculation  
        # updating the presure head and theta
        h[:,n+1] = X-z # calculation of current pressur head 
        theta[:,n+1] = Van_Genuchten_moisure( h[:,n+1], thetar, thetas , alpha, N) # calculation of current moisture head 
        water_storage[n+1] = np.trapz(theta[:,n+1], x=z, dx=dz)*1000

        if Error[n,m] <= threshold:
            break
            
# -------figures-----------------------------
# storage water level within the system
# integration of theta from 0 to L wrt to dz

# code break at 46th timestep

time_req = np.arange(tmin,tnodes,dt) 
nodes_req = np.arange(9861, 10000)
# nodes_req = np.arange(1, 1000)


plt.figure(figsize=(18, 6))
plt.plot(time_req[nodes_req],theta[29,nodes_req],label='depth = 5 cm')
plt.plot(time_req[nodes_req],theta[27,nodes_req],label='depth = 25 cm')
plt.plot(time_req[nodes_req],theta[25,nodes_req],label='depth = 45 cm') 
plt.ylim([0,0.5])      
plt.xlim([min(nodes_req),max(nodes_req)])      
plt.xlabel('Time (day)')
plt.ylabel(r"$\theta$")
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig("soil_moisture_variation.png",dpi = 200)


plt.figure(figsize=(18, 6))
plt.plot(time_req[nodes_req],water_storage[nodes_req])
plt.ylim([0,600])      
plt.xlim([min(nodes_req),max(nodes_req)])      
plt.xlabel('Time (day)')
plt.ylabel(r"$\Theta$ (mm)")
plt.title("Storage water level (mm)")
plt.grid()
plt.tight_layout()
plt.savefig("water_storage.png",dpi = 200)

