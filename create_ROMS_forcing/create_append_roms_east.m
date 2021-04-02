function create_append_roms_east(x_psi, x_rho, x_u, x_v, e_psi, e_rho, e_u, ....
e_v, s_rho, zeta_time, zeta_east, ubar_east, vbar_east, u_east, v_east, temp_east, salt_east, .....
sand_east_01, sand_east_02, sand_east_03, sand_east_04, sand_east_05, nc_bndry)
 
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

%if(obc.east==1)
 % EAST

nccreate(nc_bndry,'zeta_east','Dimensions',{'e_rho' erhodimID 'zeta_time' zttdimID },'Format','classic');
ncwrite(nc_bndry,'zeta_east',zeta_east);
ncwriteatt(nc_bndry,'zeta_east','long_name','free-surface eastern boundary condition');
ncwriteatt(nc_bndry,'zeta_east','units','meter');
ncwriteatt(nc_bndry,'zeta_east','field','zeta_east, scalar, series');

disp('after here')
nccreate(nc_bndry,'ubar_east','Dimensions',{'e_u' eudimID 'v2d_time' v2tdimID },'Format','classic');
ncwrite(nc_bndry,'ubar_east',ubar_east);
ncwriteatt(nc_bndry,'ubar_east','long_name','2D u-momentum eastern boundary condition');
ncwriteatt(nc_bndry,'ubar_east','units','meter second-1');
ncwriteatt(nc_bndry,'ubar_east','field','ubar_east, scalar, series');

nccreate(nc_bndry,'vbar_east','Dimensions',{'e_v' evdimID 'v2d_time' v2tdimID },'Format','classic');
ncwrite(nc_bndry,'vbar_east',vbar_east);
ncwriteatt(nc_bndry,'vbar_east','long_name','2D v-momentum eastern boundary condition');
ncwriteatt(nc_bndry,'vbar_east','units','meter second-1');
ncwriteatt(nc_bndry,'vbar_east','field','vbar_east, scalar, series');

nccreate(nc_bndry,'u_east','Dimensions',{'e_u' eudimID 's_rho' s_rhodimID .....
					    'v3d_time' v3tdimID },'Format','classic');
ncwrite(nc_bndry,'u_east', u_east);
ncwriteatt(nc_bndry,'u_east','long_name','3D u-momentum eastern boundary condition');
ncwriteatt(nc_bndry,'u_east','units','meter second-1');
ncwriteatt(nc_bndry,'u_east','field','u_east, scalar, series');

nccreate(nc_bndry,'v_east','Dimensions',{'e_v' evdimID 's_rho' s_rhodimID .....
					    'v3d_time' v3tdimID },'Format','classic');
ncwrite(nc_bndry,'v_east', v_east);
ncwriteatt(nc_bndry,'v_east','long_name','3D v-momentum eastern boundary condition');
ncwriteatt(nc_bndry,'v_east','units','meter second-1');
ncwriteatt(nc_bndry,'v_east','field','v_east, scalar, series');

nccreate(nc_bndry,'temp_east','Dimensions',{'e_rho' erhodimID 's_rho' s_rhodimID .....
					    'temp_time' tptdimID },'Format','classic');
ncwrite(nc_bndry,'temp_east', temp_east);
ncwriteatt(nc_bndry,'temp_east','long_name','3D temperature eastern boundary condition');
ncwriteatt(nc_bndry,'temp_east','units','C');
ncwriteatt(nc_bndry,'temp_east','field','temp_east, scalar, series');

nccreate(nc_bndry,'salt_east','Dimensions',{'e_rho' erhodimID 's_rho' s_rhodimID .....
					    'salt_time' sltdimID },'Format','classic');
ncwrite(nc_bndry,'salt_east', salt_east);
ncwriteatt(nc_bndry,'salt_east','long_name','3D salinity eastern boundary condition');
ncwriteatt(nc_bndry,'salt_east','units','psu');
ncwriteatt(nc_bndry,'salt_east','field','salt_east, scalar, series');

nccreate(nc_bndry,'sand_east_01','Dimensions',{'e_rho' erhodimID 's_rho' s_rhodimID .....
					    'sand_time' sanddimID },'Format','classic');
ncwrite(nc_bndry,'sand_east_01', sand_east_01);
ncwriteatt(nc_bndry,'sand_east_01','long_name','suspended noncohesive sediment eastern boundary condition');
ncwriteatt(nc_bndry,'sand_east_01','units','kilogram meter-3');
ncwriteatt(nc_bndry,'sand_east_01','field','sand_east_01, scalar, series');

nccreate(nc_bndry,'sand_east_02','Dimensions',{'e_rho' erhodimID 's_rho' s_rhodimID .....
					    'sand_time' sanddimID },'Format','classic');
ncwrite(nc_bndry,'sand_east_02', sand_east_02);
ncwriteatt(nc_bndry,'sand_east_02','long_name','suspended noncohesive sediment eastern boundary condition');
ncwriteatt(nc_bndry,'sand_east_02','units','kilogram meter-3');
ncwriteatt(nc_bndry,'sand_east_02','field','sand_east_02, scalar, series');

nccreate(nc_bndry,'sand_east_03','Dimensions',{'e_rho' erhodimID 's_rho' s_rhodimID .....
					    'sand_time' sanddimID },'Format','classic');
ncwrite(nc_bndry,'sand_east_03', sand_east_03);
ncwriteatt(nc_bndry,'sand_east_03','long_name','suspended noncohesive sediment eastern boundary condition');
ncwriteatt(nc_bndry,'sand_east_03','units','kilogram meter-3');
ncwriteatt(nc_bndry,'sand_east_03','field','sand_east_03, scalar, series');

nccreate(nc_bndry,'sand_east_04','Dimensions',{'e_rho' erhodimID 's_rho' s_rhodimID .....
					    'sand_time' sanddimID },'Format','classic');
ncwrite(nc_bndry,'sand_east_04', sand_east_04);
ncwriteatt(nc_bndry,'sand_east_04','long_name','suspended noncohesive sediment eastern boundary condition');
ncwriteatt(nc_bndry,'sand_east_04','units','kilogram meter-3');
ncwriteatt(nc_bndry,'sand_east_04','field','sand_east_04, scalar, series');

nccreate(nc_bndry,'sand_east_05','Dimensions',{'e_rho' erhodimID 's_rho' s_rhodimID .....
					    'sand_time' sanddimID },'Format','classic');
ncwrite(nc_bndry,'sand_east_05', sand_east_05);
ncwriteatt(nc_bndry,'sand_east_05','long_name','suspended noncohesive sediment eastern boundary condition');
ncwriteatt(nc_bndry,'sand_east_05','units','kilogram meter-3');
ncwriteatt(nc_bndry,'sand_east_05','field','sand_east_05, scalar, series');

%end
% 
% if(obc.west==1)
% % WEST
% nccreate('myfile.nc','zeta_west','Dimensions',{'erhodimId' erhodimID 'zttdimID' zttdimID },'Format','classic');
% ncwrite('myfile.nc','zeta_west',zeta_east);
% ncwriteatt('myfile.nc','zeta_west','long_name','free-surface eastern boundary condition');
% ncwriteatt('myfile.nc','zeta_west','units','meter');
% ncwriteatt('myfile.nc','zeta_west','field','zeta_east, scalar, series');
% 
% nccreate('myfile.nc','ubar_west','Dimensions',{'eudimId' eudimID 'zttdimID' v2tdimID },'Format','classic');
% ncwrite('myfile.nc','ubar_west',ubar_west);
% ncwriteatt('myfile.nc','ubar_west','long_name','2D u-momentum eastern boundary condition');
% ncwriteatt('myfile.nc','ubar_west','units','meter second-1');
% ncwriteatt('myfile.nc','ubar_west','field','ubar_west, scalar, series');
% 
% nccreate('myfile.nc','vbar_west','Dimensions',{'evdimId' evdimID 'zttdimID' v2tdimID },'Format','classic');
% ncwrite('myfile.nc','vbar_west',vbar_west);
% ncwriteatt('myfile.nc','vbar_west','long_name','2D v-momentum eastern boundary condition');
% ncwriteatt('myfile.nc','vbar_west','units','meter second-1');
% ncwriteatt('myfile.nc','vbar_west','field','vbar_west, scalar, series');
% 
% nccreate('myfile.nc','u_west','Dimensions',{'eudimId' eudimID 'srhodimID' s_rhodimID .....
% 					    'v3tdimID' v3tdimID },'Format','classic');
% ncwrite('myfile.nc','u_west', u_west);
% ncwriteatt('myfile.nc','u_west','long_name','3D u-momentum eastern boundary condition');
% ncwriteatt('myfile.nc','u_west','units','meter second-1');
% ncwriteatt('myfile.nc','u_west','field','u_west, scalar, series');
% 
% nccreate('myfile.nc','v_west','Dimensions',{'evdimId' evdimID 'srhodimID' s_rhodimID .....
% 					    'v3tdimID' v3tdimID },'Format','classic');
% ncwrite('myfile.nc','v_west', v_west);
% ncwriteatt('myfile.nc','v_west','long_name','3D v-momentum eastern boundary condition');
% ncwriteatt('myfile.nc','v_west','units','meter second-1');
% ncwriteatt('myfile.nc','v_west','field','v_west, scalar, series');
% 
% nccreate('myfile.nc','temp_west','Dimensions',{'erhodimId' erhodimID 'srhodimID' s_rhodimID .....
% 					    'tptdimID' tptdimID },'Format','classic');
% ncwrite('myfile.nc','temp_west', temp_west);
% ncwriteatt('myfile.nc','temp_west','long_name','3D temperature eastern boundary condition');
% ncwriteatt('myfile.nc','temp_west','units','C');
% ncwriteatt('myfile.nc','temp_west','field','temp_west, scalar, series');
% 
% nccreate('myfile.nc','salt_west','Dimensions',{'erhodimId' erhodimID 'srhodimID' s_rhodimID .....
% 					    'sltdimID' sltdimID },'Format','classic');
% ncwrite('myfile.nc','salt_west', salt_west);
% ncwriteatt('myfile.nc','salt_west','long_name','3D salinity eastern boundary condition');
% ncwriteatt('myfile.nc','salt_west','units','psu');
% ncwriteatt('myfile.nc','salt_west','field','salt_west, scalar, series');
% 
% nccreate('myfile.nc','sand_west_01','Dimensions',{'erhodimId' erhodimID 'srhodimID' s_rhodimID .....
% 					    'sanddimID' sanddimID },'Format','classic');
% ncwrite('myfile.nc','sand_west_01', sand_west_01);
% ncwriteatt('myfile.nc','sand_west_01','long_name','3D salinity eastern boundary condition');
% ncwriteatt('myfile.nc','sand_west_01','units','kilogram meter-3');
% ncwriteatt('myfile.nc','sand_west_01','field','sand_west_01, scalar, series');
% 
% nccreate('myfile.nc','sand_west_02','Dimensions',{'erhodimId' erhodimID 'srhodimID' s_rhodimID .....
% 					    'sanddimID' sanddimID },'Format','classic');
% ncwrite('myfile.nc','sand_west_02', sand_west_02);
% ncwriteatt('myfile.nc','sand_west_02','long_name','3D salinity eastern boundary condition');
% ncwriteatt('myfile.nc','sand_west_02','units','kilogram meter-3');
% ncwriteatt('myfile.nc','sand_west_02','field','sand_west_02, scalar, series');
% 
% nccreate('myfile.nc','sand_west_03','Dimensions',{'erhodimId' erhodimID 'srhodimID' s_rhodimID .....
% 					    'sanddimID' sanddimID },'Format','classic');
% ncwrite('myfile.nc','sand_west_03', sand_west_03);
% ncwriteatt('myfile.nc','sand_west_03','long_name','3D salinity eastern boundary condition');
% ncwriteatt('myfile.nc','sand_west_03','units','kilogram meter-3');
% ncwriteatt('myfile.nc','sand_west_03','field','sand_west_03, scalar, series');
% 
% nccreate('myfile.nc','sand_west_04','Dimensions',{'erhodimId' erhodimID 'srhodimID' s_rhodimID .....
% 					    'sanddimID' sanddimID },'Format','classic');
% ncwrite('myfile.nc','sand_west_04', sand_west_04);
% ncwriteatt('myfile.nc','sand_west_04','long_name','3D salinity eastern boundary condition');
% ncwriteatt('myfile.nc','sand_west_04','units','kilogram meter-3');
% ncwriteatt('myfile.nc','sand_west_04','field','sand_west_04, scalar, series');
% 
% nccreate('myfile.nc','sand_west_05','Dimensions',{'erhodimId' erhodimID 'srhodimID' s_rhodimID .....
% 					    'sanddimID' sanddimID },'Format','classic');
% ncwrite('myfile.nc','sand_west_05', sand_west_05);
% ncwriteatt('myfile.nc','sand_west_05','long_name','3D salinity eastern boundary condition');
% ncwriteatt('myfile.nc','sand_west_05','units','kilogram meter-3');
% ncwriteatt('myfile.nc','sand_west_05','field','sand_west_05, scalar, series');
% 
% end
% 
