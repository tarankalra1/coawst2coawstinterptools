%check init file
clear all ; close all ; clc ;
grd_fname = '/media/taran/DATADRIVE2/marsh_result/barnegat_bay/easygrid/bbleh_reedy_grd.nc';

netcdf_load(grd_fname); 

netcdf_load('reedy_init.nc')
pcolorjw(lon_rho, lat_rho, zeta)
shading flat