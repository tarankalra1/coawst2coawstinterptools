Functions to use a coarse/outer grid and its solution to create boundary forcing file 
for a refined grid.

Both are COAWST solution and grids. 

1.Boundary forcing file generation codes: 
* d_bndry_coawst2coawst.m --> This is the main file that calls for the functions to 
create boundary forcing file
Required primary functions
* create_roms_bry_from_coawst 
* create_roms_netcdf_bndry_mwUL 
* rho2u_2d_mai.m
* rho2v_2d_mai.m
* u2rho_2d_mai.m
* v2rho_2d_mai.m
* u2rho_3d_mai.m
* v2rho_3d_mai.m

2. Additional files folder
* 1. sw_dist.m --> from COAWST mtools. Define distance between two lat,lon coordinates

* 2. add_angle_to_grid.m --> will add the variable "angle" to grid 
