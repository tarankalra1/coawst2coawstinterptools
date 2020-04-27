clear all ; close all ;clc; 

grd_fname = '/media/taran/DATADRIVE2/marsh_result/barnegat_bay/easygrid/bbleh_reedy_grd.nc';
netcdf_load(grd_fname); 
h_small=h ; 
lon_rho_small=lon_rho; 
lat_rho_small=lat_rho; 

%Refined grid initial solution
netcdf_load('reedy_bry.nc')

plot(ubar_east(:,5))