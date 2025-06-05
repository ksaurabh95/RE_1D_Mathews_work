load('grid_spacing_original.mat')
zmax = 3 ; %m 
% issue is z has be measured from upwards, .i.e at surface height is 3 m. 
z_req = z_req';
dz_all = dz_all' ;
z_act = zmax-z_req ;
z_act = flip(z_act);
dz_all = flip(dz_all);
save("grid_spacing.mat","dz_all","z_act")


