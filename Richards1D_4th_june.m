% 1D richards equation based on Dogan and Motz,2005 a,b 
tic
% reading the meteological data from the excel file

meteodata = readtable("meteo_data.xlsx");
rain = meteodata(9262:end,6); % rainfall mm/day
PET = meteodata(9262:end,7); % potential ET in mm/day
rain = table2array(rain)/1000 ; % m /day
PET = table2array(PET)/1000 ; % m/day

% soil parameters
data.Ksat = 3.220;   % m/day
data.thetas = 0.3769; % saturated soil moisture content 
data.thetar = 0.0515; % residual soil moisture content
data.alpha = 3.321 ;   % air entry presure in m
data.ho = -0.19105548; % in m/day
data.N = 2.503; % van Genucten dimensionless Paramter,n  
data.Ss = 0.00; % specific storage 
data.n_eta = -0.8653 ; % relative permeability exponent, empirical exponent
% Input parameters
zmin = 0;
zmax = 3;  % in m
% dz = 0.005; % in m
% non uniform grids 
grid_data = load('grid_spacing.mat');
z = grid_data.z_act;
dz_all = grid_data.dz_all;
depth = zmax - z ;
% z is measured upwards 
% calculation of znodes
znodes = length(z);
data.z = z;
% root water uptake funtion details 
psi_a  = -0.05 ; % m
psi_d = -4 ; % m 
psi_w = -150 ; % m 

a = 2 ; % parameter describing how fast root grows
Lr = 1 ; % depth of the root zone in m 
f1 = exponential_root_distribution_function(a,Lr,depth) ; % exponential root distribution function
% time information
tmin = 1; % starting time in day
tmax = length(rain); % max time in day
time_given = (tmin:tmax)' ; 
% time step information
dt = 1 ; % in day
dt_min = 0.0000000000000000001 ; % hr
dt_max = 2 ; % hr
factor = 0.2 ;
tnodes = round( tmax/( dt_min + dt_max )  + 0.5*tmax/dt_max );  % max avialabe timesteps 

data.t = nan(tnodes,1);
data.t(1) = tmin ;
% preallocation for storing information the nodes
h=nan(znodes,tnodes);
theta=nan(znodes,tnodes);
storage_water_level = nan(tnodes,1);
RWU = nan(znodes,tnodes);
% initial condition 
depth_water_level = 2*zmax ; % assuminng water table is at 5 m depth 
h(:,1) = depth - depth_water_level ;
theta (:,1) = Van_Genuchten_theta( h(:,1), data);
storage_water_level (1) = trapz(z, theta (:,1))*1000; % estimating in mm
% Specifying the error threshold
threshold = .001;
MaxIterations = 30;
Error = nan(tnodes,MaxIterations);
residual_error = nan(tnodes,1);
time_elasped = tmin ;
n = 1;
% Solution 
while time_elasped <= tmax
    
    Repeat_itr = 0 ; % timestep required to repeated 
    Error(n,:)= nan ;
    t_req = data.t(n) + dt ;
    % use of spline interpolation 
%     q_interpolated = spline(time_given, rain, t_req); % using cubic spline to find the in between value of rain
%     PET_interpolated = spline(time_given, PET, t_req); % using cubic spline to find the in between value of ET
    
    % using step wise interpolation 
    q_interpolated = interp1(time_given, rain, t_req, 'previous', 'extrap'); % using cubic spline to find the in between value of rain
    PET_interpolated = interp1(time_given, PET, t_req, 'previous', 'extrap'); % using cubic spline to find the in between value of ET
    
    data.q(n+1) =  q_interpolated ; % to correct 
    h(:,n+1) = h(:,n);
    theta (:,n+1) = Van_Genuchten_theta(h(:,n+1),data);
    % calculaion of plant stress function 
    f2 = plant_stress_function(psi_a, psi_d, psi_w, h(:,n+1) ) ;
    % plant uptake term for the current time step 
    RWU(:,n+1) = f1.*f2.*PET_interpolated; % plant uptake term 
for m = 1:MaxIterations
        A=zeros(znodes,znodes);
        RHS=zeros(znodes,1);  
% Calculation of K for h(:,n+1)
  K = Van_Genuchten_K( h(:,n+1),data ) ;
%% for the central zone excluding upper, lower, 
 coordinates.dz = dz_all;
 coordinates.n = n;
 coordinates.znodes=znodes;
 [A,RHS] = central_zone(coordinates,data,A,RHS,h,K,theta,dt,RWU);
%% -------------lower boundary condition---------------------
% free drainage boundary condition
coordinates.z_initial = 1;
[A,RHS] = free_drainage( coordinates, data,A,RHS,h,K,theta,dt,RWU,dz_all);

%% ----------------------upper boundary condition --------
%------------- flux boundary condition ---------
coordinates.z_initial = znodes;
[A,RHS] = upper_flux_boundary(coordinates,data,A,RHS,h,K,theta,dt,RWU,dz_all);
% Reshaping the H(:,n+1) matrix such into the same shape as X1 so that
X_input = h(:,n+1)+z;
% % solution using conjugate gradient method
 tol=0.001;            % tol  : Break condition for the relative residuum
 m_max=50;             % m_max: Maximum number of iterations to perform
% [X1]=solvePCG(A, RHS,X_input,tol, m_max);
X1 = A\RHS;   
Error(n,m)=sqrt(sum((X1-X_input).^2));
residual_error(n) = Error(n,m) ;
%% updating the total head
% updating the presure head 
h(:,n+1) = X1-z;
% Updating the theata
 theta (:,n+1) = Van_Genuchten_theta(h(:,n+1),data);
 storage_water_level (n+1) = trapz(z, theta (:,n+1))*1000; % estimating in mm
 
            if Error(n,m) <= threshold 
                break
            end
                        
end

           if m >=MaxIterations 
               Repeat_itr = 1 ;
               dt = dt/3 ; 
           end
       if dt > dt_max 
              dt = dt_max ;
        elseif dt < dt_min 
              dt = dt_min ;
       end

 
   if Repeat_itr == 1
       % Repeating the Currrent timestep 
%       disp('Current Time Step is Repeated due to Reaching Max Iteration Condition ')       
     
   else
       % moving to next time step
        data.t(n+1) = data.t(n) + dt ;
        time_elasped = data.t(n+1) ;
%         clc
        if mod(n, 100) == 0
        disp(['Current time step: ', num2str(time_elasped)]);
        end

        
       % adaptive time settings [Mls, 1982; Šimunek et al., 1992]:
        if m <= 3 && dt >= dt_min && dt <= dt_max
            dt =  dt*(1+factor);  % for next time step
        elseif m >= 7  && dt >= dt_min && dt <= dt_max
             dt =  dt*(1-factor); % for next time step
        end
        % limiting dt range
        if dt > dt_max 
              dt = dt_max ;
        elseif dt < dt_min 
              dt = dt_min ;
        end
    % ensuring time_elasped doesnot increse t_max    
        time_left = tmax-time_elasped ;
        if time_left <dt 
            dt = tmax-time_elasped ;
        end   
         n = n+1 ; % moving to next timenode 
  end
  
             if time_elasped == tmax 
                 break 
             end

end
toc

storage_level_new = interp1(data.t,storage_water_level,time_given);

storage_level_new_spline = spline(data.t,storage_water_level,time_given);

day_start = 9862-9262+1 ; 
day_end = 12418 -9262+1;
idx = day_start:day_end;

startDate = datetime(1988, 1, 1);
endDate = datetime(1994, 12, 31);
dates = startDate:endDate; % Create time series

% comparison with Mathias results case 1 sandy soill 
Mathias_result = readtable("Mathias_result.xlsx");
water_storage_obs = table2array(Mathias_result(:,2)); % rainfall mm/day
Time_obs = table2array(Mathias_result(:,1)); % potential ET in mm/day
% Convert fractional years to datetime
startOfYear = datetime(floor(Time_obs), 1, 1); % Get the first day of each year
daysInYear = days(datetime(floor(Time_obs) + 1, 1, 1) - startOfYear); % Number of days in the year
fractionalPart = Time_obs - floor(Time_obs); % Extract the fractional part
dateValues = startOfYear + fractionalPart .* daysInYear; % Compute actual dates
dateValues.Format = 'yyyy-MM-dd'; % Display only date

% Display the results
% disp(dateValues);

figure()
% plot(time_given(idx), storage_level_new(idx))
plot(dates, storage_level_new_spline(idx),'-','LineWidth',2)
hold on 
plot(dateValues, water_storage_obs,'o')
ylabel("\Theta (mm)")
xlabel("Time (years)")
title("Water storage between 1988 to 1994")

% ylim([150, 310])
grid
years = year(startDate):year(endDate);
xtickDates = datetime(years, 1, 1);
xticks(xtickDates);
xticklabels(string(years));
legend("1D code","Mathias et al")

% % min_values = min(Error, [], 2);
% figure()
% plot(residual_error(7651:9701))
% xlabel("Days")

