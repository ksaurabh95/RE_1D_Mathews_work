function [ theta ] = Van_Genuchten_theta( h , data)
% calculation of soil moisture based on modified Van Genuchten and  Nelson (1985)
% thetar = residual moisture content, thetas = saturated soil moisture , h = presure head, 
% hs = air entry presure head , N = Van Genuchten paramter , Ss =  specific storage
thetas = data.thetas ; 
thetar = data.thetar ;
alpha = data.alpha  ;   % air entry presure in cm 
N = data.N ;
m = 1 - 1/N ; % van Genucten Paramter 
Se = (1 + (abs(alpha*h)).^N).^(-m);
theta =  thetar + ( thetas - thetar ).*Se ;  

end

