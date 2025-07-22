# -*- coding: utf-8 -*-
"""
Created on Tue Jul  8 14:31:56 2025

@author: saurabh kumar
"""
import numpy as np
from scipy.optimize import fsolve


def cal_Van_Genuchten_moisure( h , thetar, thetas , alpha, N):
    """
    calculation of soil moisture based on  Van Genuchten relationship (1981)
    thetar = residual moisture content, thetas = saturated soil moisture , h = presure head, 
    hs = air entry presure head , N = Van Genuchten paramter , Ss =  specific storage
    """
    y = np.empty([len(h)])
    m = 1- 1/N
    
    for i in range(len(thetas)):
        if h[i] <= 0 :
            beta = pow( abs(alpha[i]*h[i]), N[i] )
            y[i] = thetar[i] + (thetas[i] - thetar[i])* pow( 1+ beta, -m[i] )  
        else: 
            y[i] = thetas[i]

    
    return y



def Van_Genuchten_moisure( h , thetar, thetas , alpha, N):
    """
    calculation of soil moisture based on  Van Genuchten relationship (1981)
    thetar = residual moisture content, thetas = saturated soil moisture , h = presure head, 
    hs = air entry presure head , N = Van Genuchten paramter , Ss =  specific storage
    """
    m = 1- 1/N
    if h <= 0 :
        beta = pow( abs(alpha*h), N )
        y = thetar + (thetas - thetar)* pow( 1+ beta, -m )  
    else: 
        y = thetas
    
    return y


def Van_Genuchten_K( h , Ksat, thetar, thetas , alpha, N, n_eta):
    """
    calculation of soil moisture based on modified Van Genuchten and  Nelson (1985)
    thetar = residual moisture content, thetas = saturated soil moisture , h = presure head, 
    hs = air entry presure head , N = Van Genuchten paramter , Ss =  specific storage
    n_eta =  Relativepermeabilityexponen
    """
    # y = np.empty([len(h)])
    m = 1 - 1/N
    beta = pow( abs(alpha*h), N )
    Se = pow( 1+ beta,-m)
    y =  Ksat*pow(Se,n_eta)*pow( ( 1 - pow((1 - pow(Se,1/m)),m) ) , 2) 
    #if h < 0:    
     #   y =  Ksat*pow(Se,n_eta)*pow( ( 1 - pow((1 - pow(Se,1/m)),m) ) , 2) 
    # else: y = Ksat 
    
    return y



def Van_Genuchten_specific_moisure_capacity( h , thetar, thetas, alpha, N ):
    """
    calculation of soil specific moisture capacity based on 
    modified Van Genuchten and Nelson (1985)
    thetar = residual moisture content, thetas = saturated soil moisture , h = presure head, 
    hs = air entry presure head , N = Van Genuchten paramter , % Ss =  specific storage
    """
    # tempI = h > 0
    m = 1 - 1/N 
    hs = 1./alpha
    if h <= 0 :
        beta = pow( abs(alpha*h), N )
        y = (N - 1)*(thetas - thetar)*pow(abs(h),N-1 )/( pow(abs(hs),N)*pow(1 + beta, m +1) ) 
    else: 
        y = 0
        
    return y




def central_zone(A, RHS, h, thetas, thetar, alpha,Ksat,N,n_eta, RWU,dz_all, n,dt,z):
    znodes = len(dz_all)
    
    for i in range(1,znodes-1):
        
        thetas_current = thetas[i]
        thetar_current = thetar[i]
        alpha_current = alpha[i]
        N_current = N[i]
        Ksat_current = Ksat[i]
        n_eta_current = n_eta[i]
        
        d1z = dz_all[i-1]
        d2z = dz_all[i]
        d3z = dz_all[i+1]
        
        h1Z  = h[i-1,n+1]
        h2  = h[i,n+1]
        h3Z  = h[i+1,n+1]
        
        h_old = h[i,n]

        
        K1Z = Van_Genuchten_K( h1Z , Ksat_current, thetar_current, thetas_current , alpha_current, N_current, n_eta_current)
        K2  = Van_Genuchten_K( h2 , Ksat_current, thetar_current, thetas_current , alpha_current, N_current, n_eta_current)
        K3Z = Van_Genuchten_K( h3Z , Ksat_current, thetar_current, thetas_current , alpha_current, N_current, n_eta_current)
        
        
        theta_current = Van_Genuchten_moisure( h2 , thetar_current, thetas_current , alpha_current, N_current)
        theta_old = Van_Genuchten_moisure( h_old , thetar_current, thetas_current , alpha_current, N_current)
        
        C = Van_Genuchten_specific_moisure_capacity( h2 , thetar_current, thetas_current, alpha_current, N_current )
        
        # calculation of s
        s =(theta_current-theta_old)/dt          
      # calculation of p1 
        p1 = C/dt
        # calculation of p2
        p2 = 0; # Sw*Ss/dt; since Sw=0
            
        Kszminus = ( K2*d2z + K1Z*d1z ) / ( d1z+d2z )
        CNZminusmean = -Kszminus/( (d2z+d1z)/2)
        c = -CNZminusmean/d2z;             
      
        # calculation of g
        
        Kszplus = ( K2*d2z + K3Z*d3z )/ (d3z+d2z)
        CNZplusmean = -Kszplus/( (d2z+d3z)/2)
        g = -CNZplusmean/d2z  
        # calculation of d          
        d = -(  c + g + p1 + p2)
        # forming the A matrix and RHS vector
        A[i,i-1] = c   # coeffiecient of H(k-1)
        A[i,i] = d      # coeffiecient of H(k)
        A[i,i+1] = g  # coeffiecient of H(k+1)
        RHS[i] = s - p1 *( h[i,n+1] + z[i] ) - p2*( h[i,n] + z[i])  + RWU[i] # RHS at i
    
    return A,RHS


def top_flux_boundary(A, RHS,h, thetas, thetar, alpha,Ksat,N,n_eta, RWU,dz_all, n,dt,z,q_current):
    znodes = len(dz_all)
    # RHS = np.empty([znodes])   # leaving the top and boundary 
    # A = np.zeros([znodes,znodes]) # leaving the top and boundary 
    
    i = znodes-1
    thetas_current = thetas[i]
    thetar_current = thetar[i]
    alpha_current = alpha[i]
    N_current = N[i]
    Ksat_current = Ksat[i]
    n_eta_current = n_eta[i]
    
    d1z = dz_all[i-1]
    d2z = dz_all[i]
    
    h1Z  = h[i-1,n+1]
    h2  = h[i,n+1]
    
    h_old = h[i,n]

    
    K1Z = Van_Genuchten_K( h1Z , Ksat_current, thetar_current, thetas_current , alpha_current, N_current, n_eta_current)
    K2  = Van_Genuchten_K( h2 , Ksat_current, thetar_current, thetas_current , alpha_current, N_current, n_eta_current)
    
    
    theta_current = Van_Genuchten_moisure( h2 , thetar_current, thetas_current , alpha_current, N_current)
    theta_old = Van_Genuchten_moisure( h_old , thetar_current, thetas_current , alpha_current, N_current)
    
    C = Van_Genuchten_specific_moisure_capacity( h2 , thetar_current, thetas_current, alpha_current, N_current )

    # calculation of s
    s =(theta_current-theta_old)/dt          
    # calculation of p1 
    p1 = C/dt  
    # calculation of p2
    p2 = 0;         
    
    Kszminus = ( K2*d2z + K1Z*d1z ) / ( d1z+d2z )
    CNZminusmean = -Kszminus/( (d2z+d1z)/2)
    c = -CNZminusmean/d2z;             

    d = -(  c +  p1 + p2)
    
    A[i,i-1] = c   # coeffiecient of H(k-1)
    A[i,i] = d      # coeffiecient of H(k)
    RHS[i] = s - p1 *( h[i,n+1] + z[i] ) - p2*( h[i,n] + z[i]) - q_current/d2z + RWU[i]

    
    return A,RHS


def bottom_free_drainage(A, RHS, h, thetas, thetar, alpha,Ksat,N,n_eta, RWU,dz_all, n,dt,z):
    # znodes = len(dz_all)
    i = 0   # free drainage at the bottom 

    thetas_current = thetas[i]
    thetar_current = thetar[i]
    alpha_current = alpha[i]
    N_current = N[i]
    Ksat_current = Ksat[i]
    n_eta_current = n_eta[i]
    

    d2z = dz_all[i]
    d3z = dz_all[i+1]
    

    h2  = h[i,n+1]
    h3Z  = h[i+1,n+1]
    
    h_old = h[i,n]

    
    K2  = Van_Genuchten_K( h2 , Ksat_current, thetar_current, thetas_current , alpha_current, N_current, n_eta_current)
    K3Z = Van_Genuchten_K( h3Z , Ksat_current, thetar_current, thetas_current , alpha_current, N_current, n_eta_current)
    
    
    theta_current = Van_Genuchten_moisure( h2 , thetar_current, thetas_current , alpha_current, N_current)
    theta_old = Van_Genuchten_moisure( h_old , thetar_current, thetas_current , alpha_current, N_current)
    
    C = Van_Genuchten_specific_moisure_capacity( h2 , thetar_current, thetas_current, alpha_current, N_current )

    # calculation of s
    s =(theta_current-theta_old)/dt          
  # calculation of p1 
    p1 = C/dt
    # calculation of p2
    p2 = 0; 

    Kszplus = ( K2*d2z + K3Z*d3z )/ (d3z+d2z)
    CNZplusmean = -Kszplus/( (d2z+d3z)/2)
    g = -CNZplusmean/d2z 
    d = -(  g + p1 + p2)
    A[i,i] = d      # coeffiecient of H(k)
    A[i,i+1] = g  # coeffiecient of H(k+1)
    RHS[i] = s - p1 *( h[i,n+1] + z[i] ) - p2*( h[i,n] + z[i]) + K2/d2z + RWU[i]  # RHS at i
  
    return A,RHS



def bottom_no_flow(A, RHS, h, thetas, thetar, alpha,Ksat,N,n_eta, RWU,dz_all, n,dt,z):
    # znodes = len(dz_all)
    i = 0   # free drainage at the bottom 

    thetas_current = thetas[i]
    thetar_current = thetar[i]
    alpha_current = alpha[i]
    N_current = N[i]
    Ksat_current = Ksat[i]
    n_eta_current = n_eta[i]
    

    d2z = dz_all[i]
    d3z = dz_all[i+1]
    

    h2  = h[i,n+1]
    h3Z  = h[i+1,n+1]
    
    h_old = h[i,n]

    
    K2  = Van_Genuchten_K( h2 , Ksat_current, thetar_current, thetas_current , alpha_current, N_current, n_eta_current)
    K3Z = Van_Genuchten_K( h3Z , Ksat_current, thetar_current, thetas_current , alpha_current, N_current, n_eta_current)
    
    
    theta_current = Van_Genuchten_moisure( h2 , thetar_current, thetas_current , alpha_current, N_current)
    theta_old = Van_Genuchten_moisure( h_old , thetar_current, thetas_current , alpha_current, N_current)
    
    C = Van_Genuchten_specific_moisure_capacity( h2 , thetar_current, thetas_current, alpha_current, N_current )

    # calculation of s
    s =(theta_current-theta_old)/dt          
  # calculation of p1 
    p1 = C/dt
    # calculation of p2
    p2 = 0; 

    Kszplus = ( K2*d2z + K3Z*d3z )/ (d3z+d2z)
    CNZplusmean = -Kszplus/( (d2z+d3z)/2)
    g = -CNZplusmean/d2z 
    d = -(  g + p1 + p2)
    A[i,i] = d      # coeffiecient of H(k)
    A[i,i+1] = g  # coeffiecient of H(k+1)
    RHS[i] = s - p1 *( h[i,n+1] + z[i] ) - p2*( h[i,n] + z[i])  + RWU[i]  # RHS at i
  
    return A,RHS





def ponding_flux(h,Ksat, thetar, thetas , alpha, N, n_eta,dz_all,z):
    # estimation of flux at the top when ponding condition has been reached, i.e., Se >= 0.999 
    Se_top = 0.999
    
    # code to estimate the pressure head the point 
    znodes = len(dz_all)
    
    i = znodes-1
    thetas_current = thetas[i]
    thetar_current = thetar[i]
    alpha_current = alpha[i]
    N_current = N[i]
    Ksat_current = Ksat[i]
    n_eta_current = n_eta[i]

  
    theta_top = Se_top*(thetas_current - thetar_current) + thetar_current
    m = 1-1/N_current
    
    K_top = Ksat_current*(Se_top**n_eta_current)*( ( 1 - (1 - Se_top**(1/m))**m  )**2)

    f = lambda h: Van_Genuchten_moisure(h, thetar_current, thetas_current, alpha_current, N_current) - theta_top

    # Initial guess (h must be negative in unsaturated soil)
    # h_top  = fsolve(f, x0=0.2)[0]
    h_top  = 0.01

    # K_top = Van_Genuchten_K(h_top, Ksat_current, thetar_current, thetas_current , alpha_current, N_current, n_eta_current)
  
    h2  = h[i]
    K2 = Van_Genuchten_K( h2 , Ksat_current, thetar_current, thetas_current , alpha_current, N_current, n_eta_current)
    
    d2z = dz_all[i]
    Kszplus = ( K2 + K_top ) / 2 
    CNZplusmean = -Kszplus/d2z
    
  
    # Kszminus = ( K2 + K_top )  / 2 
    # Kszminus = ( K2 )
    
    H2 = h2 + z[i]
    H_top = h_top + 3 # 
    
    q = -CNZplusmean*(H_top - H2)
    # q = -CNZplusmean*(H_top - H2)

    # q = -Kszplus*( (  h_top - h2)/(d2z/2) -1  )
        
    
    return q

def f2(psi, psi_a, psi_d, psi_w):
    """
    # plant stress function
    Computes f2(psi) based on a piecewise function.

    Parameters:
    psi   : input pressure head (can be scalar or numpy array)
    psi_a : air entry pressure head
    psi_d : critical pressure head (dry threshold)
    psi_w : wilting point pressure head

    Returns:
    f2_val : result of f2(psi)
    """
    psi = np.array(psi, dtype=float)
    f2_val = np.zeros_like(psi)

    # Case: psi_a > psi > psi_d --> f2 = 1
    mask1 = (psi < psi_a) & (psi > psi_d)
    f2_val[mask1] = 1

    # Case: psi_d >= psi >= psi_w --> f2 = 1 - (psi - psi_d)/(psi_w - psi_d)
    mask2 = (psi <= psi_d) & (psi >= psi_w)
    f2_val[mask2] = 1 - (psi[mask2] - psi_d) / (psi_w - psi_d)

    # psi >= psi_a and psi < psi_w are already 0 by initialization

    return f2_val



def f1(depth,a,Lr,zmax):
    """
    exponential root distribution function
    
    """
    depth = np.array(depth, dtype=float)
    f1 = np.zeros_like(depth)
    mask = depth > Lr
    f1[~mask] = a/Lr *( np.exp(-a)-np.exp(-a*depth[~mask]/Lr) ) / ( (1+a)*np.exp(-a) -1) 
                                                                   
    f1[mask] = 0
    
    
    return f1


def top_ponding_head_boundary(A, RHS,h_pond, z):
    znodes = len(z)
    # RHS = np.empty([znodes])   # leaving the top and boundary 
    # A = np.zeros([znodes,znodes]) # leaving the top and boundary 
    
    i = znodes-1
    
    A[i,i] = 1      # coeffiecient of H(k)
    RHS[i] =  h_pond + z[i]  

    
    return A,RHS


