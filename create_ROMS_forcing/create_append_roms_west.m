function create_append_roms_west(x_psi, x_rho, x_u, x_v, e_psi, e_rho, e_u, ....
e_v, s_rho, zeta_time, zeta_west, ubar_west, vbar_west, u_west, v_west, temp_west, salt_west, .....
sand_west_01, sand_west_02, sand_west_03, sand_west_04, sand_west_05, nc_bndry)
 
%, ......
%zeta_west, ubar_west, vbar_west, u_west, v_west, temp_west, salt_west, .....
%sand_west_01, sand_west_02, sand_west_03, sand_west_04, sand_west_05, obc) 
%
% jcw April 18, 2009
% updated 01Sep2015 to have all 4 sides
%    for all the vars.
%
 L=x_psi;
LP=x_rho;
L=x_u;
LP=x_v;

M=e_psi;
MP=e_rho;
MP=e_u;
M=e_v;

s=s_rho;

erhodimID=e_rho;
eudimID=e_u;
evdimID=e_v;
s_rhodimID=s;

zttdimID=zeta_time;
v2tdimID=zeta_time;
v3tdimID=zeta_time;
tptdimID=zeta_time;
sltdimID=zeta_time ;
sanddimID=zeta_time ;
 
nc_bndry, LP

% %psidimID = netcdf.defDim(nc_bndry,'x_psi',L);
% xrhodimID = netcdf.defDim(nc_bndry,'x_rho',LP);
% xudimID = netcdf.defDim(nc_bndry,'x_u',L);
% xvdimID = netcdf.defDim(nc_bndry,'x_v',LP);

% %psidimID = netcdf.defDim(nc_bndry,'e_psi',M);
% erhodimID = netcdf.defDim(nc_bndry,'e_rho',MP);
% eudimID = netcdf.defDim(nc_bndry,'e_u',MP);
% evdimID = netcdf.defDim(nc_bndry,'e_v',M);
% s_rhodimID = netcdf.defDim(nc_bndry,'s_rho',s);
% 
% %brydimID = netcdf.defDim(nc_bndry,'bry_time',t_clim);
% 
% % tsk
% zttdimID = netcdf.defDim(nc_bndry,'zeta_time',t_clim);
% v2tdimID = netcdf.defDim(nc_bndry,'v2d_time',t_clim);
% v3tdimID = netcdf.defDim(nc_bndry,'v3d_time',t_clim);
% sltdimID = netcdf.defDim(nc_bndry,'salt_time',t_clim);
% tptdimID = netcdf.defDim(nc_bndry,'temp_time',t_clim);
% sanddimID = netcdf.defDim(nc_bndry,'sand_time',t_clim);

% tsk

%% Variables and attributes:

nccreate(nc_bndry,'zeta_west','Dimensions',{'e_rho' erhodimID 'zeta_time' zttdimID },'Format','classic');
ncwrite(nc_bndry,'zeta_west',zeta_west);
ncwriteatt(nc_bndry,'zeta_west','long_name','free-surface western boundary condition');
ncwriteatt(nc_bndry,'zeta_west','units','meter');
ncwriteatt(nc_bndry,'zeta_west','field','zeta_west, scalar, series');

disp('after here')
nccreate(nc_bndry,'ubar_west','Dimensions',{'e_u' eudimID 'v2d_time' v2tdimID },'Format','classic');
ncwrite(nc_bndry,'ubar_west',ubar_west);
ncwriteatt(nc_bndry,'ubar_west','long_name','2D u-momentum western boundary condition');
ncwriteatt(nc_bndry,'ubar_west','units','meter second-1');
ncwriteatt(nc_bndry,'ubar_west','field','ubar_west, scalar, series');

nccreate(nc_bndry,'vbar_west','Dimensions',{'e_v' evdimID 'v2d_time' v2tdimID },'Format','classic');
ncwrite(nc_bndry,'vbar_west',vbar_west);
ncwriteatt(nc_bndry,'vbar_west','long_name','2D v-momentum western boundary condition');
ncwriteatt(nc_bndry,'vbar_west','units','meter second-1');
ncwriteatt(nc_bndry,'vbar_west','field','vbar_west, scalar, series');

nccreate(nc_bndry,'u_west','Dimensions',{'e_u' eudimID 's_rho' s_rhodimID .....
					    'v3d_time' v3tdimID },'Format','classic');
ncwrite(nc_bndry,'u_west', u_west);
ncwriteatt(nc_bndry,'u_west','long_name','3D u-momentum western boundary condition');
ncwriteatt(nc_bndry,'u_west','units','meter second-1');
ncwriteatt(nc_bndry,'u_west','field','u_west, scalar, series');

nccreate(nc_bndry,'v_west','Dimensions',{'e_v' evdimID 's_rho' s_rhodimID .....
					    'v3d_time' v3tdimID },'Format','classic');
ncwrite(nc_bndry,'v_west', v_west);
ncwriteatt(nc_bndry,'v_west','long_name','3D v-momentum western boundary condition');
ncwriteatt(nc_bndry,'v_west','units','meter second-1');
ncwriteatt(nc_bndry,'v_west','field','v_west, scalar, series');

nccreate(nc_bndry,'temp_west','Dimensions',{'e_rho' erhodimID 's_rho' s_rhodimID .....
					    'temp_time' tptdimID },'Format','classic');
ncwrite(nc_bndry,'temp_west', temp_west);
ncwriteatt(nc_bndry,'temp_west','long_name','3D temperature western boundary condition');
ncwriteatt(nc_bndry,'temp_west','units','C');
ncwriteatt(nc_bndry,'temp_west','field','temp_west, scalar, series');

nccreate(nc_bndry,'salt_west','Dimensions',{'e_rho' erhodimID 's_rho' s_rhodimID .....
					    'salt_time' sltdimID },'Format','classic');
ncwrite(nc_bndry,'salt_west', salt_west);
ncwriteatt(nc_bndry,'salt_west','long_name','3D salinity western boundary condition');
ncwriteatt(nc_bndry,'salt_west','units','psu');
ncwriteatt(nc_bndry,'salt_west','field','salt_west, scalar, series');

nccreate(nc_bndry,'sand_west_01','Dimensions',{'e_rho' erhodimID 's_rho' s_rhodimID .....
					    'sand_time' sanddimID },'Format','classic');
ncwrite(nc_bndry,'sand_west_01', sand_west_01);
ncwriteatt(nc_bndry,'sand_west_01','long_name','suspended noncohesive sediment western boundary condition');
ncwriteatt(nc_bndry,'sand_west_01','units','kilogram meter-3');
ncwriteatt(nc_bndry,'sand_west_01','field','sand_west_01, scalar, series');

nccreate(nc_bndry,'sand_west_02','Dimensions',{'e_rho' erhodimID 's_rho' s_rhodimID .....
					    'sand_time' sanddimID },'Format','classic');
ncwrite(nc_bndry,'sand_west_02', sand_west_02);
ncwriteatt(nc_bndry,'sand_west_02','long_name','suspended noncohesive sediment western boundary condition');
ncwriteatt(nc_bndry,'sand_west_02','units','kilogram meter-3');
ncwriteatt(nc_bndry,'sand_west_02','field','sand_west_02, scalar, series');

nccreate(nc_bndry,'sand_west_03','Dimensions',{'e_rho' erhodimID 's_rho' s_rhodimID .....
					    'sand_time' sanddimID },'Format','classic');
ncwrite(nc_bndry,'sand_west_03', sand_west_03);
ncwriteatt(nc_bndry,'sand_west_03','long_name','suspended noncohesive sediment western boundary condition');
ncwriteatt(nc_bndry,'sand_west_03','units','kilogram meter-3');
ncwriteatt(nc_bndry,'sand_west_03','field','sand_west_03, scalar, series');

nccreate(nc_bndry,'sand_west_04','Dimensions',{'e_rho' erhodimID 's_rho' s_rhodimID .....
					    'sand_time' sanddimID },'Format','classic');
ncwrite(nc_bndry,'sand_west_04', sand_west_04);
ncwriteatt(nc_bndry,'sand_west_04','long_name','suspended noncohesive sediment western boundary condition');
ncwriteatt(nc_bndry,'sand_west_04','units','kilogram meter-3');
ncwriteatt(nc_bndry,'sand_west_04','field','sand_west_04, scalar, series');

nccreate(nc_bndry,'sand_west_05','Dimensions',{'e_rho' erhodimID 's_rho' s_rhodimID .....
					    'sand_time' sanddimID },'Format','classic');
ncwrite(nc_bndry,'sand_west_05', sand_west_05);
ncwriteatt(nc_bndry,'sand_west_05','long_name','suspended noncohesive sediment western boundary condition');
ncwriteatt(nc_bndry,'sand_west_05','units','kilogram meter-3');
ncwriteatt(nc_bndry,'sand_west_05','field','sand_west_05, scalar, series');

