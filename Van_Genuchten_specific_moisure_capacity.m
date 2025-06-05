function [ C ] = Van_Genuchten_specific_moisure_capacity( h, data)
% calculation of soil specific moisture capacity based on modified Van Genuchten and  Nelson (1985)
% thetar = residual moisture content, thetas = saturated soil moisture , h = presure head, 
% hs = air entry presure head , N = Van Genuchten paramter , % Ss =  specific storage
thetas = data.thetas ;
thetar = data.thetar ;

hs = 1/data.alpha  ;    % air entry presure in m
ho =  0 ; % in cm/h
N = data.N   ; % van Genucten Paramter 
Ss = 0 ; % specific storage in van genuchten model 
m = 1 - 1/N ; 
beta = (abs(h/hs)).^N;
 
 if h > ho
   C = Ss ;
 else
   C = ( N-1 )*( thetas - thetar )*(abs(h).^(N-1)) / ( (abs(hs).^N)*(( 1 + beta).^( m + 1)) ) ;
end

