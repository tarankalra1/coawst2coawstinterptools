The main codebase in the repository is "create_ROMS_forcing". It helps to create a forcing
file every hour or at any other frequency for forcing the bigger (coarser) model grid solution
to a smaller/inner/refined model grid solution. If the forcing is imposed every hour, the files
get larger in size. Therefore, we follow two strategies:
a) first we save the structure of 3d and 4d vars in Step 1 that
   gets replaced using scatinterp.
b) Break down the forcing to have a file generated separately for northern, eastern, western bc. 

In this case, step 2 first creates the forcing only for northern bc. That can be easily changed. 
Step 3, Step 4 append to the forcing file by incorporating eastern and western boundaries. 

--> Hacks and running tips: 
In this code, temperature and salinity were separately added
because of the nature of the problem. There were 6 sand classes. All of that can be easily modified. 

Always run the code for "1 day"  of forcing file generation and test. 

The code could also be improved to concatenate two different contiguous periods of time if forcing
files start to become "huge". For that the "time" array needs to be modified.    


Functions to use a coarse/outer grid and its solution to create boundary forcing file 
for a refined grid. Both are COAWST solution and grids. 

Boundary forcing file generation codes: 

##### 1. step1_scatinterp.m --> first call this file to create a 3d and 4d variable (in this
case zeta and temp) that can be replaced using scatinterp later. 


##### 2. step2_bndry_north.m --> Step 2 modify northern bc to use the scatinterp function
northern bc addition
--and replace the 3d and 4d vars with the desired parameters required for interpolation
(in this case, ubar, vbar, u, v, sand01-05, temp, salt)


##### 3. step2_bndry_east.m  --> Step 3 modify eastern bc to use the scatinterp function
eastern bc addition 
--and replace the 3d and 4d vars with the desired parameters required for interpolation
(in this case, ubar, vbar, u, v, sand01-05, temp, salt)

##### 4. step2_bndry_west.m --> Step 4 modify western bc to use the scatinterp function
western bc addition 
--and replace the 3d and 4d vars with the desired parameters required for interpolation
(in this case, ubar, vbar, u, v, sand01-05, temp, salt)

--------------------------
Supporting files:
--------------------------
1.  interp3d_insert_tsk.m
2.  save_westbc.m
3.  interp4d_insert_tsk.m
4.  rho2v_2d_mai.m
5.  v2rho_2d_mai.m
5.  v2rho_3d_mai.m
6.  save_4dbc_tsk.m
7.  save_eastbc.m
8.  save_westbc.m
9.  create_append_roms_east.m
10. create_append_roms_west.m
11. maplev_4dvar_tsk.m
12. u2rho_2d_mai.m
13. u2rho_3d_mai.m
14. create_roms_bry_from_coawst.m
15. rho2u_2d_mai.m
16. rho2u_3d_mai.m


* d_bndry_coawst2coawst_griddata.m --> Legacy code to create boundary forcing with grid
data (From Christie H. and Maitane O).


Items to be modifed later and would not make sense right now.

2. swan_forcing 
creating TPAR files from an existing COAWST solution
* create_TPAR_coawst.m --> Main file 

Required primary functions
*  find_nearest_point.m

3. Additional files folder
* 1. sw_dist.m --> from COAWST mtools. Define distance between two lat,lon coordinates

* 2. add_angle_to_grid.m --> will add the variable "angle" to grid 

* 3. find_veg_mask_reedy.m --> Zafers contour lines for veg used to map on a nested grid


4. Common source of error 
*  wind forcing file may have NaN's double check, rerun the wind forcing code..
 sometimes it works

* NARR data commonly used in the past was not found after 2015. But double check

* I had to use this particular mexcdf for easy grid
C:\Users\tkalra\Desktop\COAWST\COAWST_v3_3_local-master\Tools\mfiles\matlab_tools\netcdf\ncutility\mexcdf.m

