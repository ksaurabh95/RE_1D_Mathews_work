function [ A,RHS ] = upper_flux_boundary( coordinates, data,A,RHS,h,K,theta, dt,RWU,dz_all)
%Upper flux boundary condition using the dogan and motz(2005) paper  
z_initial = coordinates.z_initial;
n = coordinates.n;
% extracting the important data
q = data.q ;
z = data.z;
    k = z_initial ;
    h2  = h(k,n+1) ;
     K2  = K(k) ; 
     K1Z = K(k-1);
     theta_current = theta(k,n+1) ;
     theta_old = theta(k,n) ;
     C= Van_Genuchten_specific_moisure_capacity( h2, data ) ; 
     
     d1z = dz_all(k-1) ;
     d2z = dz_all(k) ;
     
% calculation of c
  Kszminus = ( K2*d2z + K1Z*d1z ) / (d2z+d1z);
  CNZminusmean = -Kszminus/( (d2z + d1z)/2 ) ;
  c = -CNZminusmean/d2z; 
 % calculation of p1
   p1 = C/dt;
% calculation of Sw
   % Sw = theta_current/porosity;
% calculation of p2
   p2 =0 ;  % Sw*Ss/dt;
% calculation of d
   d = -(  c + p1 + p2);
% calculation of s
   s =(theta_current-theta_old)/dt;   
% forming the A matrix and RHS vector
   A( k,k-1) = c;   % coeffiecient of H(i,j,k-1)
   A(k,k)=d;                 % coeffiecient of H(k)   
   RHS(k) = s - p1 *( h(k,n+1) + z(k)) - p2*(h(k,n)+z(k)) - q(n+1)/d2z + RWU(k,n+1);   % RHS at k    
end 

