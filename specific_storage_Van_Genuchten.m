function [ Ss ] = specific_storage_Van_Genuchten( data )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
thetas = data.thetas  ; 
thetar = data.thetar  ;
hs = 1/data.alpha   ;   % air entry presure in m
ho = data.ho ; % in cm/h
N = data.N  ; % van Genucten Paramter 
m = 1 - 1/N ; 
beta = (abs(ho/hs))^N ; 
Ss = (N-1)*(thetas-thetar)*abs(ho)/((abs(hs)^N)*(1+beta)^(m+1));
end

