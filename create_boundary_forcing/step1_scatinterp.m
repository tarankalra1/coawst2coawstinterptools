clear all ; close all ; clc ; clf; 
%  
%   This is the first code to execute before creating boundary
%   forcing file using scat interp
%   It saves 2 mat files for 3d and 4d COAWST variables 

%  - creates boundary for a COAWST grid that is refined
%    based on an existing COAWST coarse grid solution 
%    The refined grid is nested within coarse grid. 

% input -> 
% coarse grid -> grid_name_outer
% solution    -> ncml link
% refined grid -> grid_name_inner
% initial time 
% 
% Solution 
% output -> 3d and 4d var mat files  
% 
% This is currently set up to use opendap and nctoolbox.
%
% Written by: 05/20/2020 Tarandeep Kalra
% Adapted from Maitane Olabarrieta and 
% Christie Heggermiller edits
% 

% Define filenames.
% model_url : outer grid simulations
mdl_fname='https://geoport.whoi.edu/thredds/dodsC/vortexfs1/usgs/users/tkalra/barnegatsim/bbleh.ncml';
grid_name_outer = '/media/taran/DATADRIVE2/marsh_result/barnegat_bay/all_other_folders/runfiles_bbleh_zd_taran/grid/bbleh_grid_073d.nc';

setup_nctoolbox;
nc = ncgeodataset(mdl_fname);
time1 = nj_time(nc,'zeta');

%ocean_time=nc{'ocean_time'}(:); 
dstart=datenum(2015,04,15);

% These are grid parameters for the outer grid/bigger grid
% READ GRID PARAMETERS
coarse_masku = nc{'mask_u'}(:);
coarse_maskv = nc{'mask_v'}(:);
coarse_maskr = nc{'mask_rho'}(:);
angler = ncread(grid_name_outer,'angle');
angler = angler'; % Take the transpose of angler to be consistent for nc;
coarse_h = nc{'h'}(:);
coarse_Vtransform = nc{'Vtransform'}(:);
coarse_Vstretching = nc{'Vstretching'}(:);
coarse_hc = nc{'hc'}(:);
coarse_theta_s = nc{'theta_s'}(:);
coarse_theta_b = nc{'theta_b'}(:);
coarse_N = length(nc{'s_rho'}(:));
%
% This is lon and lat of the outer grid/biger grid saved 
%  
g = nc{'zeta'}(1,:,:).grid;
lon_rho = g.lon;
lat_rho = g.lat;
[nx ny]=size(lon_rho); 
grid_size_coarse=nx*ny;
%
% This is lon and lat of the outer grid/biger grid saved 
%  
lon_rho_coarse_3d=repmat(lon_rho,1,1,coarse_N);
lat_rho_coarse_3d=repmat(lat_rho,1,1,coarse_N);

close(nc);

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
ref_N   = 7;
Vtransform  = 1;
Vstretching = 1;

disp('getting roms grid dimensions ...');

Sinp.N           = ref_N;            % number of vertical levels
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
[nx_ref,ny_ref] = size(gn.lon_rho);
grid_size_ref=nx_ref*ny_ref       ;
%
% This is lon and lat of the outer grid/biger grid saved 
%  
lon_rho_ref_3d=repmat(gn.lon_rho,1,1,ref_N);
lat_rho_ref_3d=repmat(gn.lat_rho,1,1,ref_N);
zr_ref=gn.z_r; 

init_1=1; 
init_end=1; 

report=0;

%% Initialize the variables of interest
% these are temporary variables created for the interpolation to the
% smaller/refined grid ... 
tic
%% READ 3D VARIABLES FOR BRY
% 
% for scattered interolant convert arrays into 1 d vector
%Reshape all the data in column vector
lon_rho_coarse_col=reshape(lon_rho,[grid_size_coarse 1]);
lat_rho_coarse_col=reshape(lat_rho,[grid_size_coarse 1]); 
lon_rho_coarse_3d_col=reshape(lon_rho_coarse_3d,.........
                             [grid_size_coarse*coarse_N 1]);
lat_rho_coarse_3d_col=reshape(lat_rho_coarse_3d,.........
                             [grid_size_coarse*coarse_N 1]); 

lon_rho_ref_col   =reshape(gn.lon_rho,[grid_size_ref 1]);
lat_rho_ref_col   =reshape(gn.lat_rho,[grid_size_ref 1]); 
lon_rho_ref_3d_col=reshape(lon_rho_ref_3d,.........
                             [grid_size_ref*ref_N 1]);
lat_rho_ref_3d_col=reshape(lat_rho_ref_3d,.........
                             [grid_size_ref*ref_N 1]); 

for mm = init_1:init_end
    
    zeta_coarse(:,:) = double(squeeze(nc{'zeta'}(mm,:,:)));
    zeta_coarse=squeeze(zeta_coarse(:,:));
    
% Compute vertical elevations of the grid, this is time dependent
    zr_coarse = set_depth(coarse_Vtransform, coarse_Vstretching, ....
                   coarse_theta_s,    coarse_theta_b, .....
                   coarse_hc,         coarse_N,.........
                   5,                 coarse_h, zeta_coarse, report);
           
% interpolate using scatteredinterpolation 
    zeta_coarse(~coarse_maskr) = NaN;
    zeta_coarse = maplev(zeta_coarse);
    v=zeta_coarse;

%   grid size for 3D arrays
    grd_size=length(lon_rho_coarse_col);

%Reshape all the data in column vector 
    V_col=reshape(v,[grd_size 1]);

%Get the interpolant function 
    F_3d_coarse=scatteredInterpolant(lon_rho_coarse_col,lat_rho_coarse_col,V_col);

% SAVE the 3d var
    save('3dvar_bbleh_coarse.mat','F_3d_coarse')


% PROCEED TO 4D variables
% These are grid(nc temp) parameters for the outer grid/bigger grid
    temp_coarse=double(squeeze(nc{'temp'}(mm,:,:,:)));

% % 4d vars 
   for zz = 1:coarse_N
      aa = squeeze(temp_coarse(zz,:,:));
      aa(~coarse_maskr) = NaN;
      aa = maplev(aa);
      temp_tmp_coarse(:,:,zz) = aa;
      clear aa 
  end 
%     
  for zz = 1:coarse_N
      zr_coarse(:,:,zz) = maplev(squeeze(zr_coarse(:,:,zz)));
  end
 
 % Calculate the "u" velocity that is staggered on rho points.

    clear temp_coarse 

% convert the vertical coordinates in 1 d for scatinterp
   zr_coarse_3d_col=reshape( zr_coarse(:,:,1:coarse_N), ...
                             [grid_size_coarse*coarse_N 1] ) ; 

   zr_ref_3d_col   =reshape( zr_ref(:,:,1:ref_N), ....
                             [grid_size_ref*ref_N 1] ); 

   % Try the temperature variable		
   v=temp_tmp_coarse ; 
   % Saving the interpolant function for 4d vars
   grd_size=length(lon_rho_coarse_3d_col);
   %  Reshape all the data in column vector 
   V_col=reshape(v,[grd_size 1]);

  % Get the interpolant function 
   F_4d_coarse=scatteredInterpolant(lon_rho_coarse_3d_col,lat_rho_coarse_3d_col,......
                       zr_coarse_3d_col, V_col);

% SAVE the 4d var
   save('4dvar_bbleh_coarse.mat','F_4d_coarse')

end 
