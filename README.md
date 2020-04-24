Functions to use a coarse/outer grid and its solution to create boundary forcing file 
for a refined grid.

Both are COAWST solution and grids. 

1. d_obc_coawst2coawst.m --> This is the main file that calls for the functions to 
create boundary forcing file
Required primary functions
* create_roms_bry_from_coawst 
* create_roms_netcdf_bndry_mwUL 


#. sw_dist.m --> from COAWST mtools. Define distance between two lat,lon coordinates

#. add_angle_to_grid.m --> will add the variable "angle" to grid 
