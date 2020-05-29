Functions to use a coarse/outer grid and its solution to create boundary forcing file 
for a refined grid.

Both are COAWST solution and grids. 

Boundary forcing file generation codes: 

##### 1. step1_scatinterp.m --> first call this file to create a 3d and 4d variable (in this
case zeta and temp) that can be replaced using scatinterp later. 


##### 2. step2_bndry_use_scatinterp_funcs.m --> Step 2 to use the scatinterp function
and replace the 3d and 4d vars with the desired parameters required for interpolation
(in this case, ubar, vbar, u, v, sand01-05, temp, salt)

Required primary functions
* save_4dbc_tsk.m 
*  maplev_4dvar.m
* interp3d_insert_tsk.m
* interp4d_insert_tsk.m
* create_roms_bry_from_coawst 
* create_roms_netcdf_bndry_mwUL 
* rho2u_2d_mai.m
* rho2v_2d_mai.m
* u2rho_2d_mai.m
* v2rho_2d_mai.m
* u2rho_3d_mai.m
* v2rho_3d_mai.m

* d_bndry_coawst2coawst_griddata.m --> Legacy code to create boundary forcing with grid
data (From Christie H. and Maitane O).

2. swan_forcing 
creating TPAR files from an existing COAWST solution
* create_TPAR_coawst.m --> Main file 

Required primary functions
*  find_nearest_point.m

3. Additional files folder
* 1. sw_dist.m --> from COAWST mtools. Define distance between two lat,lon coordinates

* 2. add_angle_to_grid.m --> will add the variable "angle" to grid 
