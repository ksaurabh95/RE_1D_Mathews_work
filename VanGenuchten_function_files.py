# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 17:14:47 2021

@author: Richa_Pc01
"""

"""
----------------------------defining required Van genuchten functions
"""

import numpy as np


def Van_Genuchten_moisure( h , thetar, thetas , alpha, N):
    """
    calculation of soil moisture based on modified Van Genuchten and  Nelson (1985)
    thetar = residual moisture content, thetas = saturated soil moisture , h = presure head, 
    hs = air entry presure head , N = Van Genuchten paramter , Ss =  specific storage
    """
    # y = np.empty([len(h)])
    m = 1 - 1/N
    beta = pow( abs(alpha*h), N )
    y =  thetar + ( thetas - thetar )* pow( 1+ beta,-m) 
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
    return y


def specific_moisture_capacity(h, thetas, thetar, alpha, N):
    """
    Calculate the specific moisture capacity C(h) using van Genuchten parameters.

    Parameters:
    h      : Pressure head (negative for unsaturated zone), scalar or numpy array
    thetas : Saturated water content
    thetar : Residual water content
    alpha  : van Genuchten alpha parameter (1/cm)
    n      : van Genuchten n parameter (>1)

    Returns:
    C(h)   : Specific moisture capacity
    """
    m = 1 - 1/ N

    # Only compute for h < 0 (unsaturated conditions)

    term = (alpha * h) ** N
    denom = (1 + term) ** (m + 1)
    C = (thetas - thetar) * alpha * N * m * (alpha * h) ** (N - 1) / denom

    return C



def Van_Genuchten_specific_moisure_capacity( h , thetar, thetas, alpha, N ):
    """
    calculation of soil specific moisture capacity based on 
    modified Van Genuchten and Nelson (1985)
    thetar = residual moisture content, thetas = saturated soil moisture , h = presure head, 
    hs = air entry presure head , N = Van Genuchten paramter , % Ss =  specific storage
    """
    y = np.empty([len(h)])
    tempI = h > 0
    m = 1 - 1/N 
    hs = 1/alpha
    beta = pow( abs(alpha*h[~tempI]), N )

    y[~tempI] = (N-1)*(thetas - thetar)* pow(abs(h[~tempI]),N-1)/( pow(abs(hs),N)*pow(1 + beta, m+1) )   
    y[tempI] = 0
            
    return y


# def Feddes_function( h , psi_a, psi_d, psi_w ):
#     """
#     calculation of soil specific moisture capacity based on 
#     modified Van Genuchten and Nelson (1985)
#     thetar = residual moisture content, thetas = saturated soil moisture , h = presure head, 
#     hs = air entry presure head , N = Van Genuchten paramter , % Ss =  specific storage
#     """
    
            
#     return y



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










def Van_Genuchten_and_Nelson_moisure( x , thetar, thetas , hs , ho , N, Ss):
    """
    calculation of soil moisture based on modified Van Genuchten and  Nelson (1985)
    thetar = residual moisture content, thetas = saturated soil moisture , h = presure head, 
    hs = air entry presure head , N = Van Genuchten paramter , Ss =  specific storage
    """
    y = np.empty([len(x)])
    tempI = x > ho
    m = 1 - 1/N
    beta = pow( abs(x[~tempI]/hs), N )
    betao = pow(ho/hs,N)
    y[~tempI] =  thetar + ( thetas - thetar )* pow( 1+ beta,-m) 
    y[tempI] =  thetar + ( thetas - thetar )*pow( 1 + betao, -m) + Ss*( x[tempI] - ho )   
    return y


def Van_Genuchten_and_Nelson_hydraulic_conductivity(x, Ksat, hs, N):
    """
    calculation of soil moisture based on modified Van Genuchten and  Nelson (1985)
    h = presure head, 
    hs = air entry presure head , N = Van Genuchten paramter , 
    Ksat = saturated hyraulic conductivity ,
    """
    y = np.empty([len(x)])
    tempI = x < 0
    m = 1 - 1/N
    beta = pow( abs(x[tempI]/hs), N )
    y[tempI] =Ksat * (  pow(1+beta,-5*m/2))* pow( pow(1 + beta, m) -  pow(beta,m) , 2 ) 
    y[~tempI] =  Ksat  
    return y

def Van_Genuchten_and_Nelson_specific_moisure_capacity( x , thetar, thetas, hs, ho, N, Ss ):
    """
    calculation of soil specific moisture capacity based on 
    modified Van Genuchten and Nelson (1985)
    thetar = residual moisture content, thetas = saturated soil moisture , h = presure head, 
    hs = air entry presure head , N = Van Genuchten paramter , % Ss =  specific storage
    """
    y = np.empty([len(x)])
    tempI = x > ho 
    m = 1 - 1/N 
    beta = pow( abs(x[~tempI]/hs), N )
    y[tempI] = Ss
    y[~tempI] = (N-1)*(thetas - thetar)* pow(abs(x[~tempI]),N-1)/( pow(abs(hs),N)*pow(1 + beta, m+1) )   
    return y

def Van_Genuchten_and_Nelson_C( x , thetar, thetas , hs , ho , N, Ss ):
    """
    calculation of soil specific moisture capacity based on 
    modified Van Genuchten and Nelson (1985)
    thetar = residual moisture content, thetas = saturated soil moisture , h = presure head, 
    hs = air entry presure head , N = Van Genuchten paramter , % Ss =  specific storage
    """
    if x > ho:
        y = Ss
    else:
        m = 1 - 1/N 
        beta = pow( abs(x/hs), N )
        y = (N-1)*(thetas - thetar)* pow(abs(x),N-1)/( pow(abs(hs),N)*pow(1+ beta, m+1) )   
    return y

