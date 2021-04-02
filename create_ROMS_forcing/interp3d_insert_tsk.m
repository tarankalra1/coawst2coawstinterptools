function [F]=interp3d_insert_tsk(lon_rho_coarse_col, lat_rho_coarse_col, ....
                             v_3d, F_3d_coarse, .....
                             lon_rho_ref_col, lat_rho_ref_col, .....
                             nx_ref, ny_ref )


% v_3d=input coarse grid parameter (zeta, u , v)
% grd_size= coarse grid 
%F_3d_coarse= loaded from step 1 

grd_size=length(lon_rho_coarse_col);
%Reshape all the data in column vector 

V_3d_col=reshape(v_3d,[grd_size 1]) ;

F_3d_coarse.Values = V_3d_col;
    
F_interp=F_3d_coarse(lon_rho_ref_col, lat_rho_ref_col);
    
F=reshape(F_interp, nx_ref, ny_ref);
    
clear V_3d_col F_interp 
 

