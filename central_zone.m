function [ A,RHS ] = central_zone( coordinates, data,A,RHS,h,K,theta,dt,RWU)
% Detailed explanation goes here
dz_all=coordinates.dz;
n = coordinates.n;
znodes = coordinates.znodes;
z = data.z;
% forming the A and RHS matrix
 for k=2:znodes-1
     h2  = h(k,n+1) ;
     K2  = K(k) ;
     K1Z = K(k-1) ;
     K3Z = K(k+1) ;
     theta_current = theta(k,n+1) ;
     theta_old = theta(k,n) ;
     C= Van_Genuchten_specific_moisure_capacity( h2, data) ; 
     d1z = dz_all(k-1) ;
     d2z = dz_all(k) ;
     d3z = dz_all(k+1) ;

     
     
  % calculation of paramters  
  % calculation of c
  Kszminus = ( K2*d2z + K1Z*d1z ) / (d2z+d1z);
  CNZminusmean = -Kszminus/( (d2z + d1z)/2 ) ;
  c = -CNZminusmean/d2z;
  % calculation of g
  Kszplus = ( K2*d2z + K3Z*d3z )/(d2z + d3z);
  CNZplusmean = -Kszplus/( (d2z + d3z)/2 );
  g = -CNZplusmean/d2z;   
  % calculation of p1
  p1 = C/dt;
  % calculation of p2
  p2 = 0; % Sw*Ss/dt; since Ss=0
  % calculation of d
  d = -(  c + g + p1 + p2);
  % calculation of s
  s =(theta_current-theta_old)/dt;
  % forming the A matrix and RHS vector
  A( k,k-1) = c;   % coeffiecient of H(k-1)
  A( k,k) = d;          % coeffiecient of H(k)
  A( k,k+1) = g;   % coeffiecient of H(k+1)
  RHS(k) = s - p1 *( h(k,n+1) + z(k) ) - p2*( h(k,n) + z(k)) + RWU(k,n+1)   ;       % RHS at k     
end
end

