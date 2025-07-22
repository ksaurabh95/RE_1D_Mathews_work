# -*- coding: utf-8 -*-
"""
Created on Fri Jul 11 13:10:36 2025

@author: saurabh kumar
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Jul  9 10:51:13 2025

@author: saurabh kumar
"""


"""
@author: Saurabh Kumar 
Code based on Dogan and Mortz 2002a,b 
# Dogan, A., & Motz, L. H. (2005). Saturated-unsaturated 3D groundwater model. I: Development. 
# Dogan, A., & Motz, L. H. (2005). Saturated-unsaturated 3D groundwater model. II: Verification and application

"""
import numpy as np
import matplotlib.pyplot as plt
from RE_codes_depth_variation import f1,f2
import pandas as pd
import scipy.io
from scipy.interpolate import interp1d, CubicSpline
from RE_codes_depth_variation import central_zone,top_flux_boundary,bottom_free_drainage,ponding_flux,top_ponding_head_boundary
from RE_codes_depth_variation import cal_Van_Genuchten_moisure,Van_Genuchten_K,Van_Genuchten_moisure,Van_Genuchten_specific_moisure_capacity
from scipy.optimize import fsolve

import time


t0= time.perf_counter()
np.set_printoptions(precision=17)

# importing meterological data
df = pd.read_excel('meteo_data.xlsx')  # For .xlsx files
# df = df.iloc[9000:12782] # reading from year 1990 to 1997 
df = df.iloc[9496:12782] # from 1st jan 1987 to 31st dec 1995 

df = df.reset_index()
df = df.drop(["index"], axis=1)
df['datetime'] = pd.to_datetime(df[['year', 'month', 'day']])

rain = np.array(df["rainfall_mm_per_d"])  # mm/day
PET = np.array(df["PE_mm_per_d"])   # mm/day
rain = rain/1000 # m/day
PET = PET/1000 # m/day potential evaporation



# non uniform grid spacing, we have generated the data and directly importing it
zmin = 0  # in m 
zmax = 3  # in m
grid_data = scipy.io.loadmat("grid_spacing.mat")
dz_all = grid_data.get("dz_all")
z = grid_data.get("z_act")
z = z[:,0]
dz_all = dz_all[:,0]
znodes = len(z) # nos of nodes 

# run time details and interval  
tmin = 1; # starting time in day
tmax = len(rain); # max time in day
dt = 1; #in day  
dt_min = 0.0000000000001 # min time step 
dt_max = 1 # max time step 
plus_factor = 0.12  # factor with which time step dt will either increase or decrease
dec_factor = 0.12 # decrease factor for time step
time_given = np.arange(tmin,tmax+dt,dt)
# tnodes = len(rain)

tnodes = 100*round( 2*tmax/( dt_min + dt_max )  + 0.5*tmax/dt_max )  #  max avialabe timesteps 
t = np.empty([tnodes]);
# t = np.nan
t[0] = tmin 


# to import SHPs 
Ksat = 108.5/1000  *np.ones(znodes)
thetas = 0.4616 *np.ones(znodes)
thetar = 0.0961 *np.ones(znodes)
alpha = 2.711*np.ones(znodes)
N = 1.149*np.ones(znodes)
Ss = 0.00 #  specific storage 
n_eta = -5.153*np.ones(znodes)
q = rain  # flux in m/day at the top boundary 


# preallocate h and theta 
h = np.empty([znodes,tnodes])
theta=np.empty([znodes,tnodes])
# storage water level within the system
water_storage = np.empty([tnodes])
water_storage_100_cm = np.empty([tnodes])  # water content in top 1 m depth
Ea = np.empty([tnodes])
q_ponding_flux = np.empty([tnodes])
q_rain = np.empty([tnodes])
q_simulation  = np.empty([tnodes])

# root water uptake function 

psi_a = -0.05 # critical pressure heads associated with anaerobiosis,
psi_d = -4 # critical pressure heads associated with soilwater-limited evapotranspiration
psi_w = -150 # # critical pressure head associated with plant wilting

Lr = 1 # m # depth of root zone
a = 2 # m # measuremnet of how fast root zone declines with depth
depth = zmax- z # z is taken upwards positive, depth is being measured from the surface
f1_vals = f1(depth,a,Lr,zmax) # 1st term of the root water uptake model by Feddes

# initial condition 
z_wt = 3*zmax # m. assumption the water table is situated at depth 3 times of z_max
h[:,0] = depth - z_wt # given 
# initial moisture condition
theta[:,0] = cal_Van_Genuchten_moisure( h[:,0] , thetar, thetas , alpha, N)
water_storage[0] = np.trapezoid(theta[:,0], x=z, dx=dz_all)*1000

nodes_1m_depth = np.arange(10,znodes)
water_storage_100_cm[0] = np.trapezoid(theta[nodes_1m_depth,0], x=z[nodes_1m_depth])*1000

# stopping criteria
threshold = .001
MaxIterations = 30
Error = np.zeros([tnodes,MaxIterations]) # storing RMSE error at each time step and iteration 
residual = np.zeros([tnodes])
# begin simulation 
n = 0  # initialize the timestep counter explicitly    

# Creating the interpolation function with step-wise behavior
rain_interp_func = interp1d(time_given, q, kind='next', bounds_error=False, fill_value='extrapolate')
Ep_interp_func = interp1d(time_given, PET, kind='next', bounds_error=False, fill_value='extrapolate')

# rain_interp_func = CubicSpline(time_given, q, bc_type='natural')
# Ep_interp_func = CubicSpline(time_given, Ep, bc_type='natural')

while n < tnodes:
    
    h[:,n+1] = h[:,n]
    theta [:,n+1] = theta[:,n]
    Error[n,:] = np.nan
    current_time = t[n] + dt# unknown time for which prediction will be done

    # # Interpolated rainfall and ET values
    q_rain[n+1] = rain_interp_func(current_time)
    Ep_current = Ep_interp_func(current_time)
    
     
    for m in range(MaxIterations):
        RHS = np.empty([znodes])
        A = np.zeros([znodes,znodes])        
        f2_vals = f2(h[:,n+1], psi_a, psi_d, psi_w)
        RWU = f1_vals*f2_vals*Ep_current
        Ea[n+1] =  np.trapezoid(RWU[nodes_1m_depth], x=z[nodes_1m_depth])*1000 
        
        q_ponding_flux[n+1] = ponding_flux(h[:,n+1],Ksat, thetar, thetas , alpha, N, n_eta,dz_all,z)
        # check for the ponding condition 
        if q_rain[n+1] > q_ponding_flux[n+1]:
            q_i = q_ponding_flux[n+1]
            # print('ponding flux')
        else: q_i = q_rain[n+1]
            
        q_simulation[n+1] = q_i

             
        A, RHS = bottom_free_drainage(A, RHS,h,thetas, thetar, alpha,Ksat,N,n_eta, RWU,dz_all, n,dt,z) 
        A, RHS = central_zone(A, RHS, h, thetas, thetar, alpha,Ksat,N,n_eta, RWU,dz_all, n,dt,z)
        A, RHS = top_flux_boundary(A, RHS,h, thetas, thetar, alpha,Ksat,N,n_eta, RWU,dz_all, n,dt,z,q_i) 
        
            
        X_input = h[:,n+1]+z

        try:
            X = np.linalg.solve(A, RHS)
            # X = np.linalg.inv(A).dot(RHS)
        except np.linalg.LinAlgError:  # to supress error and simulation break
            Repeat_itr = True
            dt = dt*(1-plus_factor)

            n = n -1# time step will not change 
            break  # Exit the iteration loop to repeat the timestep with smaller dt

        Error[n,m] =  np.sqrt(np.sum(pow(X-X_input,2)))  # RMSE
        residual [n] = Error[n,m]
        
        # updating the presure head and theta
        h[:,n+1] = X-z
        theta[:,n+1] = cal_Van_Genuchten_moisure( h[:,n+1], thetar, thetas , alpha, N)
        water_storage[n+1] = np.trapezoid(theta[:,n+1], x=z)*1000
        
        water_storage_100_cm [n+1] = np.trapezoid(theta[nodes_1m_depth,n+1], x=z[nodes_1m_depth])*1000 
        
        # stopping criteria for the iteration, if threshold is reached, move to the next iteration 
        if Error[n,m] <= threshold:             
            break
        

            
        
        
    else: # repeating the loop if only iteration is not repeated
        Repeat_itr = True
        # break 
        n = n -1# time step will not change 
        dt = dt / 3 # current value of dt is reduced by a third
        # print("repeat")
    

    t[n + 1] = t[n] + dt # moving to next time step
    # n = n + 1  # moving to next time node
    time_elasped = t[n + 1]
    
    if n % 300 == 0:    
        print("simulation run time =",time_elasped)
        
    
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
print("Time elapsed: ", t1) # CPU seconds elapsed 
# removing the extra nodes from the required variables 
t = t[0:n+1]
residual = residual[0:n+1]
q_rain = q_rain[0:n+1]
q_ponding_flux = q_ponding_flux[0:n+1]
water_storage = water_storage[0:n+1] 
water_storage_100_cm = water_storage_100_cm[0:n+1] 

h = h[:, 0:n+1] 
theta = theta[:, 0:n+1] 
RWU = RWU[0:n+1] 
Error=  Error[0:n+1,:]



        


# -------figures-----------------------------
time_req = t[0:n+1] # np.arange(tmin,tnodes,dt) 

time_req = time_req
water_storage_req = water_storage[0:n+1]
# plt.figure()
# plt.plot(time_req,water_storage_req)
# plt.ylabel(r"$\Theta$ (mm)")
# plt.title("Daily water storage level from 2000 to 2025")
# # plt.ylim([510,1050])
# plt.show()

water_storage_100_cm_spline= CubicSpline(t, water_storage_100_cm, bc_type='natural')
water_storage_100_cm_final = water_storage_100_cm_spline(time_given)

water_storage_spline= CubicSpline(t, water_storage, bc_type='natural')
water_storage_final = water_storage_spline(time_given)

df['water_storage_final'] =water_storage_final
df['water_storage_100_cm_final'] =water_storage_100_cm_final


# comparison of predcited values with Mathis result 
import datetime as dt

df['datetime'] = pd.to_datetime(df[['year', 'month', 'day']])
df['water_storage_final'] =water_storage_final

start_date = "1988-01-01"
end_date = "1995-01-01"  # Adjust if you get data beyond 1995
df_filtered = df[(df["datetime"] >= start_date) & (df["datetime"] <= end_date)]

# comparison of predcited values with Mathis result 

# Read Excel file (assuming column name is 'Time')
Mathias = pd.read_excel("Mathias_result.xlsx",sheet_name="clay")

# Convert fractional years to datetime
years = Mathias["Time"].astype(int)  # Extract integer year
fractional_part = Mathias["Time"] - years  # Get the fractional component

# Compute actual dates
start_of_year = pd.to_datetime(years.astype(str) + "-01-01")  # Start of the year
days_in_year = (pd.to_datetime((years + 1).astype(str) + "-01-01") - start_of_year).dt.days
Mathias["Date"] = start_of_year + pd.to_timedelta(fractional_part * days_in_year, unit="D")

# Format to show only date
Mathias["Date"] = Mathias["Date"].dt.strftime('%Y-%m-%d')
Mathias["Date"] = pd.to_datetime(Mathias["Date"])  # Convert if



# Plot the daily rainfall
plt.figure(figsize=(15,3))
plt.plot(df_filtered["datetime"], df_filtered["water_storage_final"],label="RE Model")
plt.plot(Mathias["Date"], Mathias["Water_Storage"],"o" ,label="Mathias et al")
# plt.xlabel("Date")
plt.ylabel(r"$\Theta$ (mm)")
plt.title("Daily water storage level from 1988 to 1994")
plt.ylim([900,1400])
plt.yticks(range(1000, 1410, 200))
plt.xticks(pd.date_range(start=start_date, end=end_date, freq="YS"), rotation=0)  # 'YS' stands for Year Start
plt.xlim(pd.Timestamp(start_date), pd.Timestamp(end_date))  # Set x-axis limits

plt.grid()
plt.tight_layout()
plt.legend()
# plt.savefig("water_storage_level_clay_soil.png",dpi = 300)
# plt.savefig("water_storage_level_clay_soil.svg",dpi = 300)

plt.show()

# df_filtered.to_excel("filtered_data.xlsx", index=False)

# Mathias.to_excel("Mathias_soil.xlsx", index=False)










# df.to_excel("filtered_data_depth_variation.xlsx", index=False)

# scipy.io.savemat('filtered_data_depth_variation.mat', {'water_storage_100_cm_final': water_storage_100_cm_final, 
#                                                        'water_storage_final': water_storage_final,'residual': residual,'q_rain': q_rain,
#                                                        'q_ponding_flux': q_ponding_flux,'h': h,'theta': theta,'RWU': RWU,
#                                                        'Error': Error,'water_storage': water_storage,
#                                                        'water_storage_100_cm': water_storage_100_cm,'time_given': time_given })


plt.figure()

plt.plot(t[0:n+1],q_ponding_flux[0:n+1],label ='ponding flux')
plt.plot(t[0:n+1],q_rain[0:n+1],label ='rain')
plt.plot(t[0:n+1],q_simulation[0:n+1],label ='act')
plt.legend()
plt.show()

# plt.figure()

# plt.plot(t[0:n+1],theta[ 35, 0:n+1])
# plt.show()


q_ponding_flux[ q_ponding_flux <=q_rain]










