
"""
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
from scipy.interpolate import CubicSpline

import time
t0= time.perf_counter()

# importing meterological data
df = pd.read_excel('meteo_data.xlsx')  # For .xlsx files

df = df.iloc[10592:len(df)] # reading from year 1990 to 1997 
rain = np.array(df["rainfall_mm_per_d"])  # mm/day
PE = np.array(df["PE_mm_per_d"])   # mm/day
rain = rain/1000 # m/day
Ep = PE/1000 # m/day potential evaporation


# non uniform grid spacing, we have generated the data and directly importing it
zmin = 0  # in m 
zmax = 3  # in m
grid_data = scipy.io.loadmat("grid_spacing.mat")
dz_all = grid_data.get("dz_all")[0]
z = grid_data.get("z_req")[0]
znodes = len(z) # nos of nodes 

# run time details and interval  
tmin = 1; # starting time in day
tmax = len(rain); # max time in day
dt = 1; #in day  
dt_min = 0.000000001 # min time step 
dt_max = 2 # max time step 
plus_factor = 0.12  # factor with which time step dt will either increase or decrease
dec_factor = 0.2 # decrease factor for time step
time_given = np.arange(tmin,tmax+dt,dt)
# tnodes = len(rain)

tnodes = 3*round( 2*tmax/( dt_min + dt_max )  + 0.5*tmax/dt_max )  #  max avialabe timesteps 
t = np.empty([tnodes]);
# t = np.NaN
t[0] = tmin 
 

Ksat = 3.220   # m/day
thetas = 0.3769 
thetar = 0.0515
alpha = 3.321    # air entry presure in m
N = 2.503 # van Genucten Paramter 
Ss = 0.00 #  specific storage 
n_eta =  -0.8653 # Relative permeability exponen

q = rain  # flux in m/day at the top boundary 


# preallocate h and theta 
h = np.empty([znodes,tnodes])
theta=np.empty([znodes,tnodes])
# storage water level within the system
water_storage = np.empty([tnodes])
Ea = np.empty([znodes,tnodes])


# root water uptake function 

psi_a = -0.05 # critical pressure heads associated with anaerobiosis,
psi_d = -4 # critical pressure heads associated with soilwater-limited evapotranspiration
psi_w = -150 # # critical pressure head associated with plant wilting

Lr = 1 # m # depth of root zone
a = 2 # m # measuremnet of how fast root zone declines with depth
depth = zmax- z # z is taken upwards positive, depth is being measured from the surface
f1_vals = f1(depth,a,Lr,zmax) # 1st term of the root water uptake model by Feddes

# initial condition 
z_wt = 10 # m. assumption the water table is situated at depth greater than z_max
h[:,0] = depth - z_wt # given 
# initial moisture condition
theta[:,0] = Van_Genuchten_moisure( h[:,0] , thetar, thetas , alpha, N)
water_storage[0] = np.trapz(theta[:,0], x=z, dx=dz_all)*1000

# stopping criteria
threshold = .001
MaxIterations = 20
Error = np.zeros([tnodes,MaxIterations]) # storing RMSE error at each time step and iteration 

# begin simulation 
n = 0  # initialize the timestep counter explicitly    

while n < tnodes:
    
    Repeat_itr = False  # timestep required to repeated due to max iteration condition 
    # Repeat_itr = False

    current_time = t[n] + dt # unknown time for which prediction will be done
    print(current_time)
    # estimating the values of Ea and q 
    Ep_current = np.interp(current_time, time_given, Ep)  # linear interpolation
    q_current = np.interp(current_time, time_given, q)
    # using cubic spline to interpolate the values 
    # rain_spline = CubicSpline(time_given, q, bc_type='natural')
    # ET_spline = CubicSpline(time_given, Ep, bc_type='natural')
        
#    # Interpolated rainfall and ET values
    # q_current = rain_spline(current_time)
    # Ep_current = ET_spline(current_time)

# # Avoid negative interpolated rainfall

    h[:,n+1] = h[:,n]
    theta [:,n+1] = theta[:,n]
    # print(n)
    Error[n,:] = np.NaN
    for m in range(MaxIterations):
        RHS = np.empty([znodes])
        A = np.zeros([znodes,znodes])
        K = Van_Genuchten_K(h[:,n+1],Ksat, thetar, thetas , alpha, N, n_eta)
        C = Van_Genuchten_specific_moisure_capacity(h[:,n+1], thetar, thetas, alpha, N )
        
        f2_vals = f2(h[:,n+1], psi_a, psi_d, psi_w)
        RWU = f1_vals*f2_vals*Ep_current
        Ea[:,n+1] = RWU
         
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
            
            if i ==0 :  # free drainage at the bottom 
                K3Z = K[i+1]
                d3z = dz_all[i+1]
                Kszplus = ( K2*d2z + K3Z*d3z )/( (d2z+d3z)/2)
                CNZplusmean = -Kszplus/( (d2z+d3z)/2)
                g = -CNZplusmean/d2z 
                d = -(  g + p1 + p2)
                A[i,i] = d      # coeffiecient of H(k)
                A[i,i+1] = g  # coeffiecient of H(k+1)
                RHS[i] = s - p1 *( h[i,n+1] + z[i] ) - p2*( h[i,n] + z[i]) - K2/d2z + RWU[i]/d2z  # RHS at i

           
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
                RHS[i] = s - p1 *( h[i,n+1] + z[i] ) - p2*( h[i,n] + z[i])  + RWU[i]/d2z # RHS at i
            elif i == znodes-1:   # surface boundary condition, only flux boundary condition
                K1Z = K[i-1] 
                Kszminus = ( K2*d2z + K1Z*d1z )  / ( (d1z+d2z)/2 )
                CNZminusmean = -Kszminus/( (d2z+d1z)/2)
                c = -CNZminusmean/d2z;             
                # calculation of p1 
                p1 = C[i]/dt  
                d = -(  c +  p1 + p2)
                A[i,i-1] = c   # coeffiecient of H(k-1)
                A[i,i] = d      # coeffiecient of H(k)
                RHS[i] = s - p1 *( h[i,n+1] + z[i] ) - p2*( h[i,n] + z[i]) - q_current/d2z + RWU[i]/d2z
           
            
        X_input = h[:,n+1]+z
        # X = np.linalg.inv(A).dot(RHS)
        # X = np.linalg.solve(A,RHS)
        # X = np.linalg.lstsq(A,RHS, rcond=None)[0] 
        try:
            # X = np.linalg.solve(A, RHS)
            X = np.linalg.inv(A).dot(RHS)
        except np.linalg.LinAlgError:
    # print(f"Singular matrix encountered at timestep {n}, iteration {m}. Reducing dt and retrying...")
            # dt = max(dt / 3, dt_min)
            Repeat_itr = True
            dt = dt*(1-plus_factor)

            n = n -1# time step will not change 
            break  # Exit the iteration loop to repeat the timestep with smaller dt

        Error[n,m] =  np.sqrt(np.mean(pow(X-X_input,2)))
        
        # updating the presure head and theta
        h[:,n+1] = X-z
        theta[:,n+1] = Van_Genuchten_moisure( h[:,n+1], thetar, thetas , alpha, N)
        water_storage[n+1] = np.trapz(theta[:,n+1], x=z, dx=dz_all)*1000
        
        # stopping criteria for the iteration, if threshold is reached, move to the next iteration 
        if Error[n,m] <= threshold:
            break
        
    else: # repeating the loop if only iteration is not repeated
        Repeat_itr = True
        # break 
        n = n -1# time step will not change 
        dt = dt / 3 # current value of dt is reduced by a third
        print("repeat")
    

    t[n + 1] = t[n] + dt # moving to next time step
    # n = n + 1  # moving to next time node
    time_elasped = t[n + 1]
    n = n + 1 # moving to the next iteration     
        
        # adaptive time settings [Mls, 1982; Šimunek et al., 1992]:
    if m <= 3 and dt >= dt_min and dt <= dt_max:
        dt = dt * (1 + plus_factor)  # for next time step
    elif m >= 7 and dt >= dt_min and dt <= dt_max:
        dt = dt * (1 - dec_factor)  # for next time step
    # limiting dt range
    dt = np.clip(dt, dt_min, dt_max)
    
    # ensuring time_elasped does not increase t_max
    time_left = tmax - time_elasped
    if time_left < dt:
        dt = time_left
        
    if time_elasped == tmax:
        break


t1 = time.perf_counter() - t0
print("Time elapsed: ", t1) # CPU seconds elapsed (floating point)





plt.figure()
plt.plot(t[0:n], water_storage[0:n])
plt.show()





        
# -------figures-----------------------------
time_req = t # np.arange(tmin,tnodes,dt) 
# nodes_req = np.arange(9861, 10000)
# nodes_req = np.arange(1, 1000)
fig, (ax1, ax2,ax3) = plt.subplots(3, 1, figsize=(12, 14))
ax1.plot(time_req,theta[29,:],label='depth = 5 cm')
ax1.plot(time_req,theta[27,:],label='depth = 25 cm')
ax1.plot(time_req,theta[25,:],label='depth = 45 cm') 
ax1.set_ylim([0,0.5])      
# ax1.set_xlim([min(nodes_req),max(nodes_req)])      
ax1.set_xlabel('Time (day)')
ax1.set_ylabel(r"$\theta$")
ax1.set_title("Soil Moisture Variation betweeen ")
ax1.legend()



ax2.plot(time_req,water_storage)
ax2.set_ylim([0,600])      
# ax2.set_xlim([min(nodes_req),max(nodes_req)])      
ax2.set_xlabel('Time (day)')
ax2.set_ylabel(r"$\Theta$ (mm)")
ax2.set_title("Storage water level (mm)")
ax2.grid()

ax3.plot(rain[0:round(max(t))])
# ax3.set_ylim([0,60])      
# ax2.set_xlim([min(nodes_req),max(nodes_req)])      
ax3.set_xlabel('Time (day)')
ax3.set_ylabel("rain (mm/day)")
ax3.set_title("rainfall data")
ax3.grid()


plt.tight_layout()
plt.show()
# plt.savefig("non_uniform_grid_soil_moisture_variation.png", dpi= 300)
# plt.savefig("non_uniform_grid_soil_moisture_variation.svg", dpi= 300)
