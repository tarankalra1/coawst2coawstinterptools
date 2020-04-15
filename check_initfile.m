%check init file
clear all ; close all ; clc ;
% Written by Tarandeep S Kalra to check the refined grid initial file that 
% was created by COAWST coarse grid and initial file 
% using d_init_coawst2coawst.m 

%% Outer/Coarse grid and initial file
grid_name_outer = '/media/taran/DATADRIVE2/marsh_result/barnegat_bay/all_other_folders/runfiles_bbleh_zd_taran/grid/bbleh_grid_073d.nc';
netcdf_load(grid_name_outer); 
lon_rho_outer=lon_rho; 
lat_rho_outer=lat_rho; 
h_outer=h ; 

init_outer='../init_taran_bbleh_2.nc';
netcdf_load(init_outer); 
zeta_outer=zeta; 
h_outer=h ;

%% Refined grid 
grd_fname = '/media/taran/DATADRIVE2/marsh_result/barnegat_bay/easygrid/bbleh_reedy_grd.nc';
netcdf_load(grd_fname); 
h_small=h ; 
lon_rho_small=lon_rho; 
lat_rho_small=lat_rho; 

%Refined grid initial solution
netcdf_load('reedy_init.nc')
zeta_small=zeta; 
%netcdf_load('reedy_bry.nc'); 

figure(1) 
subplot(1,2,1)
pcolor(lon_rho_outer, lat_rho_outer, h_outer)
grid off
hold on 
pcolorjw(lon_rho_small, lat_rho_small, h_small)
colorbar
caxis([-3 3])
ylim([40 40.06])
xlim([-74.1 -74.05])
title('coarse + refined grid bathy')
 
subplot(1,2,2)
pcolorjw(lon_rho_outer, lat_rho_outer, h_outer)
%grid of
%hold on 
%pcolorjw(lon_rho_small, lat_rho_small, h_small)
colorbar
caxis([-3 3])
ylim([40 40.06])
xlim([-74.1 -74.05])
title('coarse grid bathy')

figure(2) 
subplot(1,2,1)
pcolor(lon_rho_outer, lat_rho_outer, zeta_outer)
grid off
hold on 
pcolorjw(lon_rho_small, lat_rho_small, zeta_small)
colorbar
caxis([-1 1])
ylim([40 40.06])
xlim([-74.1 -74.05])
title('coarse + refined grid zeta')
 
subplot(1,2,2)
pcolorjw(lon_rho_outer, lat_rho_outer, zeta_outer)
%grid off
%hold on 
%pcolorjw(lon_rho_small, lat_rho_small, h_small)
colorbar
caxis([-1 1])
ylim([40 40.06])
xlim([-74.1 -74.05])
title('coarse grid zeta')

