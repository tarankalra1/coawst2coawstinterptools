% This routine :
%  - creates boundary for a COAWST grid that is refined
%    based on an existing COAWST coarse grid solution 
%    The refined grid is nested within coarse grid. 

% input -> 
% coarse grid -> grid_name_outer
% solution    -> ncml link
% refined grid -> grid_name_inner
%
% Solution 
% output -> boundary forcing file 
% 
%
% bry file includes zeta, ubar, vbar, sand, salt, temp
%
% This is currently set up to use opendap and nctoolbox.
%
% written by Maitane Olabarrieta 01/14/2018
% Christie Heggermiller edits
% Tarandeep Kalra adding sediment bc's 

clear all; close all;
 clf ;
 clc ; 
% Define filenames.
% model_url : outer grid simulations
%%%%%%%%%%%%%OUTER/COARSE GRID%%%%%%%%%%%%%%%%%%%%%%%%%%
%mdl_fname = 'http://geoport.whoi.edu/thredds/dodsC/sand/usgs/users/tkalra/bbleh/bbleh.ncml';

%           http://geoport.whoi.edu/thredds/dodsC/sand/usgs/users/tkalra/bbleh_run_2/61_91/bbleh.ncml
mdl_fname='http://geoport.whoi.edu/thredds/dodsC/sand/usgs/users/tkalra/bbleh_run_2/61_91/bbleh.ncml';
grid_name_outer = '/media/taran/DATADRIVE2/marsh_result/barnegat_bay/all_other_folders/runfiles_bbleh_zd_taran/grid/bbleh_grid_073d.nc';

setup_nctoolbox;
nc = ncgeodataset(mdl_fname);
time1 = nj_time(nc,'zeta');

%ocean_time=nc{'ocean_time'}(:); 
dstart=datenum(2015,04,15);

%time=nc{'ocean_time'}(:); 
t_insec(:)=nc{'ocean_time'}(:);
t_indays(:)=t_insec(:)/(3600*24);
% Converting time in days to julian time units;
t(:)=t_indays(:) +dstart;

% The time period for which the TPAR files need to be generated
t1=datenum(2015,06,17);
t2=datenum(2015,06,17,1,0,0);
init_1=near(t,t1);
init_end=near(t,t2);

time=t(init_1:init_end);


% These are grid parameters for the outer grid/bigger grid
% READ GRID PARAMETERS
masku = nc{'mask_u'}(:);
maskv = nc{'mask_v'}(:);
maskr = nc{'mask_rho'}(:);
angler = ncread(grid_name_outer,'angle');
angler = angler'; % Take the transpose of angler to be consistent for nc;
coawst_h = nc{'h'}(:);
coawst_Vtransform = nc{'Vtransform'}(:);
coawst_Vstretching = nc{'Vstretching'}(:);
coawst_hc = nc{'hc'}(:);
coawst_theta_s = nc{'theta_s'}(:);
coawst_theta_b = nc{'theta_b'}(:);
coawst_N = length(nc{'s_rho'}(:));
%
% This is lon and lat of the outer grid/biger grid saved 
%  
g = nc{'zeta'}(1,:,:).grid;
lon_rho = g.lon;
lat_rho = g.lat;
%
% This is lon and lat of the outer grid/biger grid saved 
%  
lon_rho_z = repmat(lon_rho,1,1,coawst_N);
lat_rho_z = repmat(lat_rho,1,1,coawst_N);

close(nc);

%%%%%%%%%%%%%INNER/REFINED GRID%%%%%%%%%%%%%%%%%%%%%%%%%%
% These are grid parameters for the refined grid/inner grid
% gridname : refined grid
grid_name_inner = '/media/taran/DATADRIVE2/marsh_result/barnegat_bay/easygrid/bbleh_reedy_grd.nc';
% initial name for refined grid that is user input 
init_fname=  'reedy_init.nc'
% bathymetry name for refined grid that is user input 
bry_fname=  'reedy_bry.nc'

% These parameters are for the refined/inner grid 
% Enter grid vertical coordinate parameters
% These need to be consistent with the refined ROMS setup.
theta_s = 0.0;
theta_b = 0.0;
Tcline  = 0;
N       = 7;
Vtransform  = 1;
Vstretching = 1;

disp('getting roms grid dimensions ...');

Sinp.N           = N;            % number of vertical levels
Sinp.Vtransform  = Vtransform;   % vertical transformation equation
Sinp.Vstretching = Vstretching;  % vertical stretching function
Sinp.theta_s     = theta_s;      % surface control parameter
Sinp.theta_b     = theta_b;      % bottom  control parameter
Sinp.Tcline      = Tcline;       % surface/bottom stretching width

if Vtransform == 1
    h = ncread(grid_name_inner,'h');
    hmin = min(h(:));
    hc = min(max(hmin,0),Tcline);
elseif Vtransform == 2
    h = ncread(grid_name_inner,'h');
    hmin = max(0.1,min(h(:)));
    hc = Tcline;
end

 % inputs pertaining to refined grids 
Sinp.hc = hc;                          % stretching width used in ROMS
gn = get_roms_grid(grid_name_inner,Sinp);
[nxr,nyr] = size(gn.lon_rho);
[nxu,nyu] = size(gn.lon_u);
[nxv,nyv] = size(gn.lon_v);

%%%%%%%%%%%%%%%%%Choose time steps for forcing file %%%%%%%%%%%%%%%%%%%%%%%5
%T =3;% length(time);
 
u = zeros(nxu,nyu,N);
v = zeros(nxv,nyv,N);
temp = zeros(nxr,nyr,N);
salt = zeros(nxr,nyr,N);
ubar = zeros(nxu,nyu,init_end);
vbar = zeros(nxv,nyv,init_end);
zeta = zeros(nxr,nyr,init_end);

% Save the lon and lat of smaller/refined grid for interpolation 
lon_rho_romsz = repmat(gn.lon_rho,1,1,N);
lat_rho_romsz = repmat(gn.lat_rho,1,1,N);


% Save 3D variables for interpolation from outer/bigger in temporary arrays
% READ VARIABLES

% this is required for calculation of water level 
report=0;

%% Initialize the variables of interest
% these are temporary variables created for the interpolation to the
% smaller/refined grid ... 

%% READ 3D VARIABLES FOR BRY
% 
for mm = init_1:init_end
    
    zeta_coawst(:,:)   = double(squeeze(nc{'zeta'}(mm,:,:)));
    aa = squeeze(zeta_coawst(:,:));
    size(aa)
    % Compute vertical elevations of the grid, this is time dependent
    % These depths are negative
 %   if mod(mm-1,24) == 0;
    zr = set_depth(coawst_Vtransform, coawst_Vstretching, ....
                   coawst_theta_s,    coawst_theta_b, .....
                   coawst_hc,               coawst_N,.........
                   5,                 coawst_h, aa, report);
               
%         zr = set_depth(coawst_Vtransform,coawst_Vstretching,coawst_theta_s,.......
%                        coawst_theta_b,coawst_hc,coawst_N, ...
%                        1,coawst_h,aa);
%    end
    
% interpolate using griddata
    aa(~maskr) = NaN;
    aa = maplev(aa);
    zz = griddata(lon_rho(:),lat_rho(:),aa(:),gn.lon_rho,gn.lat_rho);
% This will be zeta for refined mesh
    zeta(:,:,mm-init_1+1) = zz;
    
% PROCEED TO 4D variables
% These are grid(nc temp) parameters for the outer grid/bigger grid
    sand01_coawst=double(squeeze(nc{'sand_01'}(mm,:,:,:)));
    sand02_coawst=double(squeeze(nc{'sand_02'}(mm,:,:,:)));

% 4d vars 
    for zz = 1:coawst_N
        aa = squeeze(sand01_coawst(zz,:,:));
        aa(~maskr) = NaN;
        aa = maplev(aa);
        sand01_tmp_coarse(:,:,zz) = aa;
        clear aa
        
        aa = squeeze(sand02_coawst(zz,:,:));
        aa(~maskr) = NaN;
        aa = maplev(aa);
        sand02_tmp_coarse(:,:,zz) = aa;
        clear aa
    end 
    
    for zz = 1:coawst_N
        zr(:,:,zz) = maplev(squeeze(zr(:,:,zz)));
    end
    
    clear ur_coarse vr_coarse u2_tmp_coarse v2_tmp_coarse
    clear temp_coawst salt_coawst 
    clear sand01_coawst sand02_coawst sand03_coawst
    clear sand04_coawst sand05_coawst
    clear u_coawst v_coawst  
    
    % These are grid(nc temp) parameters for the outer grid/bigger grid
%  temp is the inner/refined grid
    sand01_tmp_ref= griddata(lon_rho_z,         lat_rho_z,      zr(:,:,1:7),.....
                             sand01_tmp_coarse, lon_rho_romsz,  lat_rho_romsz, gn.z_r);

    sand02_tmp_ref= griddata(lon_rho_z,         lat_rho_z,      zr(:,:,1:7),.....
                             sand02_tmp_coarse, lon_rho_romsz,  lat_rho_romsz, gn.z_r);
   
    clear ur_ref vr_ref ur2_ref  vr2_ref ur_tmp_coarse vr_tmp_coarse
    clear temp_coarse salt_coarse
    clear sand01_tmp_coarse sand02_tmp_coarse clear sand03_tmp_coarse 
    clear sand04_tmp_coarse sand05_tmp_coarse

       % Remove possible NaN values
    for zz = 1:N
        aa = squeeze(sand01_tmp_ref(:,:,zz));
        aa = maplev(aa);
        sand01_tmp_ref(:,:,zz) = aa;
        clear aa
        
        aa = squeeze(sand02_tmp_ref(:,:,zz));
        aa = maplev(aa);
        sand02_tmp_ref(:,:,zz) = aa;
        clear aa
    end 
% temp_4d is the refined/inner grid
    sand01_4d(:,:,:,mm-init_1+1)=sand01_tmp_ref(:,:,:);
    sand02_4d(:,:,:,mm-init_1+1)=sand02_tmp_ref(:,:,:); 
    
end
toc
clear temp_ref salt_ref sand01_tmp_ref sand02_tmp_ref
clear sand03_tmp_ref sand04_temp_ref sand05_tmp_ref
clear u_tmp_ref v_tmp_ref

% store local counters
count_1=1;
count_end=init_end-init_1+1 ; 

% %% CREATE BOUNDARY CONDITION
% 
%time=t(init_1:init_end); 

obc.north=false;
obc.south=false;
obc.east=true;
obc.west=false;

%if(inorth=='true')
%if(obc.north=true)
  zeta_north = squeeze(zeta(:, end,      count_1:count_end));
  ubar_north = squeeze(ubar(:, end,      count_1:count_end));
  vbar_north = squeeze(vbar(:, end,      count_1:count_end));
  
  temp_north = squeeze(temp_4d(:, end, :,count_1:count_end));
  salt_north = squeeze(salt_4d(:, end, :,count_1:count_end));
  u_north=squeeze(u_4d(:, end,:,         count_1:count_end));
  v_north=squeeze(v_4d(:, end,:,         count_1:count_end));
  
  sand_north_01=squeeze(sand01_4d(:,end,:,count_1:count_end));
  sand_north_02=squeeze(sand02_4d(:,end,:,count_1:count_end));
  sand_north_03=squeeze(sand03_4d(:,end,:,count_1:count_end));
  sand_north_04=squeeze(sand04_4d(:,end,:,count_1:count_end));
  sand_north_05=squeeze(sand05_4d(:,end,:,count_1:count_end));

%elseif(isouth='true')
  zeta_south = squeeze(zeta(:,  1,        count_1:count_end));
  ubar_south = squeeze(ubar(:,  1,        count_1:count_end));
  vbar_south = squeeze(vbar(:,  1,        count_1:count_end));
  
  temp_south = squeeze(temp_4d(:,  1,  :, count_1:count_end));
  salt_south = squeeze(salt_4d(:,  1,  :, count_1:count_end));
  u_south=squeeze(u_4d(:,  1, :,          count_1:count_end));
  v_south=squeeze(v_4d(:,  1, :,          count_1:count_end));

  sand_south_01=squeeze(sand01_4d(:,1,:, count_1:count_end));
  sand_south_02=squeeze(sand02_4d(:,1,:, count_1:count_end));
  sand_south_03=squeeze(sand03_4d(:,1,:, count_1:count_end));
  sand_south_04=squeeze(sand04_4d(:,1,:, count_1:count_end));
  sand_south_05=squeeze(sand05_4d(:,1,:, count_1:count_end));

%elseif(ieast='true')  
  zeta_east = squeeze(zeta(end, :,       count_1:count_end));
  ubar_east = squeeze(ubar(end, :,       count_1:count_end));
  vbar_east = squeeze(vbar(end, :,       count_1:count_end));
  temp_east  = squeeze(temp_4d(end,:,  :,count_1:count_end));
  salt_east  = squeeze(salt_4d(end,:,  :,count_1:count_end));
  u_east=squeeze (u_4d(end, :, :,        count_1:count_end));
  v_east=squeeze (v_4d(end, :, :,        count_1:count_end));

  sand_east_01=squeeze(sand01_4d(end,:,:,count_1:count_end));
  sand_east_02=squeeze(sand02_4d(end,:,:,count_1:count_end));
  sand_east_03=squeeze(sand03_4d(end,:,:,count_1:count_end));
  sand_east_04=squeeze(sand04_4d(end,:,:,count_1:count_end));
  sand_east_05=squeeze(sand05_4d(end,:,:,count_1:count_end));

%elseif(iwest='true')  
  zeta_west = squeeze(zeta(1,   :,       count_1:count_end));
  ubar_west = squeeze(ubar(1,   :,       count_1:count_end)); 
  vbar_west = squeeze(vbar(1,   :,       count_1:count_end)); 
  temp_west  = squeeze(temp_4d(1,  :,  :,count_1:count_end)); 
  salt_west  = squeeze(salt_4d(1,  :,  :,count_1:count_end));
  u_west=squeeze (u_4d(1,   :, :,        count_1:count_end));
  v_west=squeeze (v_4d(1,   :, :,        count_1:count_end));

  
  sand_west_01=squeeze(sand01_4d(1,:,:,count_1:count_end));
  sand_west_02=squeeze(sand02_4d(1,:,:,count_1:count_end));
  sand_west_03=squeeze(sand03_4d(1,:,:,count_1:count_end));
  sand_west_04=squeeze(sand04_4d(1,:,:,count_1:count_end));
  sand_west_05=squeeze(sand05_4d(1,:,:,count_1:count_end));
%end 

create_roms_bry_from_coawst(grid_name_inner, bry_fname, time,...
    Sinp.theta_s,Sinp.theta_b,Sinp.Tcline,Sinp.Vtransform,Sinp.Vstretching,Sinp.N,...
     zeta_north, ubar_north, vbar_north,...
    temp_north, salt_north, u_north, v_north, .....
    sand_north_01, sand_north_02, sand_north_03, sand_north_04, .....
    sand_north_05,......
     zeta_south, ubar_south, vbar_south,...
    temp_south, salt_south, u_south, v_south, .....
    sand_south_01, sand_south_02, sand_south_03, sand_south_04, .....
    sand_south_05,......
    zeta_east, ubar_east, vbar_east,...
    temp_east, salt_east, u_east, v_east, .....
    sand_east_01, sand_east_02, sand_east_03, sand_east_04, .....
    sand_east_05,......
    zeta_west, ubar_west, vbar_west,...
    temp_west, salt_west, u_west, v_west, .....
    sand_west_01, sand_west_02, sand_west_03, sand_west_04, .....
    sand_west_05, obc ) 

%  elseif(isouth='true')    

%  elseif(isouth='true')    
    
%  elseif(isouth='true')    
    
    
%    temp_south, temp_east, temp_west, ....
% sand_north_02,     ubar_south, vbar_south,zeta_south,...
% sand_north_03,     zeta_east, ubar_east,vbar_east,...
% sand_north_04,     zeta_west,ubar_west,vbar_west, ....
% sand_north_05,     temp_north, temp_south, temp_east, temp_west, ....
%    salt_north, salt_south, salt_east, salt_west, ....
%    u_north, u_south, u_east, u_west, .....
%    v_north, v_south, v_east, v_west,......
%   sand_south_01, sand_east_01, sand_west_01,....
%   sand_south_02, sand_east_02, sand_west_02,....
%   sand_south_03, sand_east_03, sand_west_03,....
%   sand_south_04, sand_east_04, sand_west_04,....
%   sand_south_05, sand_east_05, sand_west_05);

