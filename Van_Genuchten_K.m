function [ K ] = Van_Genuchten_K( h,data )
% calculation of hydraulic permeablity based on Van Genuchten (1985)
% Ksat = saturated hyraulic conductivity , h = presure head, 
% hs = air entry presure head , N = Van Genuchten paramter 
Ksat = data.Ksat ; 
alpha = data.alpha  ;   % air entry presure in m 
N = data.N ;
n_eta = data.n_eta ;
m = 1 - 1/N ; % van Genucten Paramter 

% tempI = h < 0;
Se = (1 + (abs(alpha*h)).^N).^(-m);
K = Ksat.*(Se.^n_eta).*(1 - (1 - Se.^(1/m) ).^m ).^2 ;

% K(tempI) =  Ksat * ((1+beta).^(-5*m/2)).* ( ( (1 + beta).^m - beta.^m).^2 ) ;
% K(~tempI) = Ksat ;
end

