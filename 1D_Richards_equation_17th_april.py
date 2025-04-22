# -*- coding: utf-8 -*-
"""
Created on Sat Mar 27 15:06:44 2021

@author: Saurabh Kumar 
Code based on Dogan and Mortz 2002a,b 
# Dogan, A., & Motz, L. H. (2005). Saturated-unsaturated 3D groundwater model. I: Development. 
# Dogan, A., & Motz, L. H. (2005). Saturated-unsaturated 3D groundwater model. II: Verification and application

"""
import numpy as np
import matplotlib.pyplot as plt
from VanGenuchten_function_files import f1,f2, Van_Genuchten_moisure,Van_Genuchten_K,Van_Genuchten_specific_moisure_capacity
import pandas as pd
import scipy.io




df = pd.read_excel('meteo_data.xlsx')  # For .xlsx files
rain = np.array(df["rainfall(mm/d)"])  # mm/day
PE = np.array(df["PE(mm/d)"])   # mm/day
rain = rain/1000 # m/day
Ep = PE/1000 # m/day potential evaporation


zmin = 0
zmax = 3  # in m
# dz = 0.1 # in m
# z = np.arange(zmin + dz/2,zmax+dz/2,dz)

# spacing is non uniform 
z_1 = 0
z_2 = 1
z_3 = 3
dz_2 = 0.1818
dz_i = 0.0022

# z_2 = np.arange(3,1,-0.1818)
# z_2_nodes = len(z_2)
# z_1_nodes = z_nodes - z_2_nodes 

grid_data = scipy.io.loadmat("grid_spacing.mat")
dz_all = grid_data.get("dz_all")[0]
z = grid_data.get("z_req")[0]



tmin = 0; # starting time in day
tmax = len(rain); # max time in day
dt = 1; #in day
t = np.arange(tmin,tmax+dt,dt) 

Ksat = 3.220   # m/day
thetas = 0.3769 
thetar = 0.0515
alpha = 3.321    # air entry presure in m
N = 2.503 # van Genucten Paramter 
Ss = 0.00 #  specific storage 
n_eta =  -0.8653 # Relativepermeabilityexponen

q = Ep - rain  # flux in m/day at the top boundary 
tnodes = len(t)
znodes = len(z)
# preallocate h and theta 
h = np.empty([znodes,tnodes])
theta=np.empty([znodes,tnodes])
# storage water level within the system
water_storage = np.empty([tnodes])
z_wt = 5 # m. the water table is situated at depth greater than z_max


rho_pu = 600.6 # storage capacity at root water uptake, mm
rho_w = 174.0  # water_level at wiltng point,  mm

# root water uptake function 

psi_a = -0.05 # critical pressure heads associated with anaerobiosis,
psi_d = -4 # critical pressure heads associated with soilwater-limited evapotranspiration
psi_w = -150 # # critical pressure head associated with plant wilting

Lr = 1 # m # depth of root zone
a = 2 # m # measuremnet of how fast root zone declines with depth
depth = zmax- z

f1_vals = f1(depth,a,Lr,zmax)

# initial presure condition 
h[:,0] = depth - z_wt 
# initial moisture condition
theta[:,0] = Van_Genuchten_moisure( h[:,0] , thetar, thetas , alpha, N)



threshold = .001
MaxIterations = 50
Error = np.zeros([tnodes,MaxIterations])


RHS = np.empty([znodes])
A = np.empty([znodes,znodes])

#K = Van_Genuchten_K(h[:,0], Ksat, hs, N)


# def free_drainage_boundary(A, RHS, ):
    
#     return A, RHS

# tnodes = 750





for n in range(tnodes-1):
    h[:,n+1] = h[:,n]
    theta [:,n+1] = theta[:,n]
    
    for m in range(MaxIterations):
        RHS = np.empty([znodes])
        A = np.zeros([znodes,znodes])
        K = Van_Genuchten_K(h[:,n+1],Ksat, thetar, thetas , alpha, N, n_eta)
        C = Van_Genuchten_specific_moisure_capacity(h[:,n+1], thetar, thetas, alpha, N )
        
        f2_vals = f2(h[:,n+1], psi_a, psi_d, psi_w)
        RWU = f1_vals*f2_vals*Ep[n+1]

         
        for i in range(0,znodes):

            d2z = dz_all[i]

            
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
            
            if i ==0 :  # free drainage 
                K3Z = K[i+1]
                d3z = dz_all[i+1]
                Kszplus = ( K2*d2z + K3Z*d3z )/( (d2z+d3z)/2)
                CNZplusmean = -Kszplus/( (d2z+d3z)/2)
                g = -CNZplusmean/d2z 
                d = -(  g + p1 + p2)
                A[i,i] = d      # coeffiecient of H(k)
                A[i,i+1] = g  # coeffiecient of H(k+1)
                RHS[i] = s - p1 *( h[i,n+1] + z[i] ) - p2*( h[i,n] + z[i]) - K2/d2z - RWU[i]/d2z  # RHS at i

           
            elif i < znodes-1: # central zone
                d1z = dz_all[i-1]
                d3z = dz_all[i+1]
                K1Z = K[i-1] 
                Kszminus = ( K2*d2z + K1Z*d1z ) / ( (d1z+d2z)/2 )
                CNZminusmean = -Kszminus/( (d2z+d1z)/2)
                c = -CNZminusmean/d2z;             
                # calculation of p1 
    #            C = Van_Genuchten_and_Nelson_C( h2 , thetar, thetas , hs , ho , N, Ss )
                p1 = C[i]/dt           
                # calculation of g
                K3Z = K[i+1]
                Kszplus = ( K2*d2z + K3Z*d3z )/( (d3z+d2z)/2 )
                CNZplusmean = -Kszplus/( (d2z+d3z)/2)
                g = -CNZplusmean/d2z  
                # calculation of d          
                d = -(  c + g + p1 + p2)
                # forming the A matrix and RHS vector
                A[i,i-1] = c   # coeffiecient of H(k-1)
                A[i,i] = d      # coeffiecient of H(k)
                A[i,i+1] = g  # coeffiecient of H(k+1)
                RHS[i] = s - p1 *( h[i,n+1] + z[i] ) - p2*( h[i,n] + z[i])  - RWU[i]/d2z # RHS at i
            elif i == znodes-1:   # surface boundary condition, only flux boundary condition
                K1Z = K[i-1] 
                Kszminus = ( K2*d2z + K1Z*d1z )  / ( (d1z+d2z)/2 )
                CNZminusmean = -Kszminus/( (d2z+d1z)/2)
                c = -CNZminusmean/d2z;             
                # calculation of p1 
    #            C = Van_Genuchten_and_Nelson_C( h2 , thetar, thetas , hs , ho , N, Ss )
                p1 = C[i]/dt
  
                d = -(  c +  p1 + p2)
                A[i,i-1] = c   # coeffiecient of H(k-1)
                A[i,i] = d      # coeffiecient of H(k)
                RHS[i] = s - p1 *( h[i,n+1] + z[i] ) - p2*( h[i,n] + z[i]) + q[n+1]/d2z - RWU[i]/d2z 
           
            
        X_input = h[:,n+1]+z
        # X = np.linalg.inv(A).dot(RHS)
        # X = np.linalg.solve(A,RHS)
        # X = np.linalg.lstsq(A,RHS, rcond=None)[0] 
        X = np.linalg.lstsq(A,RHS,rcond=None)[0] 
        Error[n,m] =  np.sqrt(np.sum(pow(X-X_input,2)))
        # updating the presure head and theta
        h[:,n+1] = X-z
        theta[:,n+1] = Van_Genuchten_moisure( h[:,n+1], thetar, thetas , alpha, N)
        water_storage[n+1] = np.trapz(theta[:,n+1], x=z, dx=dz_all)*1000

        if Error[n,m] <= threshold:
            break
            
# -------figures-----------------------------
# storage water level within the system
# integration of theta from 0 to L wrt to dz

# code break at 46th timestep

time_req = np.arange(tmin,tnodes,dt) 
nodes_req = np.arange(9861, tnodes)

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



plt.figure(figsize=(18, 6))
plt.plot(time_req[nodes_req],water_storage[nodes_req])
plt.ylim([0,1000])      
plt.xlim([min(nodes_req),max(nodes_req)])      
plt.xlabel('Time (day)')
plt.ylabel(r"$\Theta$ (mm)")
plt.title("Storage water level (mm)")
plt.legend()
plt.grid()
plt.tight_layout()

plt.figure()
plt.plot(rain*1000)
plt.ylabel('Rain (mm)')

plt.figure()
plt.plot(Ep*1000)
plt.ylabel('E_p (mm/day)')

plt.figure()
plt.plot(Ep*1000)


