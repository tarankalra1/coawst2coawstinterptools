% This routine :
%  - creates boundary and initial condition files for ROMS:
%    coawst_bdy.nc ; coawst_ini.nc
%    on a user-defined grid for a user-defined date.
%
% bry file includes zeta, ubar, vbar
% ini file includes zeta, u, v, ubar, vbar, temp, salt
%
% This is currently set up to use opendap and nctoolbox.
%
% written by Maitane Olabarrieta 01/14/2018
% Christie Heggermiller edits
% Tarandeep Kalra adding sediment bc's 

clear; close all;

% Define filenames.
% model_url : outer grid simulations
mdl_fname = 'http://geoport.whoi.edu/thredds/dodsC/sand/usgs/users/tkalra/bbleh/bbleh.ncml';
grid_name_outer = '/media/taran/DATADRIVE2/marsh_result/barnegat_bay/all_other_folders/runfiles_bbleh_zd_taran/grid/bbleh_grid_073d.nc';

setup_nctoolbox;
nc = ncgeodataset(mdl_fname);
time = nj_time(nc,'zeta');

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

% These are grid parameters for the refined grid/inner grid
% gridname : refined grid
grd_fname = '/media/taran/DATADRIVE2/marsh_result/barnegat_bay/easygrid/bbleh_reedy_grd.nc';
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
N       = 10;
Vtransform  = 1;
Vstretching = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('getting roms grid dimensions ...');

Sinp.N           = N;            % number of vertical levels
Sinp.Vtransform  = Vtransform;   % vertical transformation equation
Sinp.Vstretching = Vstretching;  % vertical stretching function
Sinp.theta_s     = theta_s;      % surface control parameter
Sinp.theta_b     = theta_b;      % bottom  control parameter
Sinp.Tcline      = Tcline;       % surface/bottom stretching width

if Vtransform == 1
    h = ncread(grd_fname,'h');
    hmin = min(h(:));
    hc = min(max(hmin,0),Tcline);
elseif Vtransform == 2
    h = ncread(grd_fname,'h');
    hmin = max(0.1,min(h(:)));
    hc = Tcline;
end

 % inputs pertaining to refined grids 
Sinp.hc = hc;                          % stretching width used in ROMS
gn = get_roms_grid(grd_fname,Sinp);
[nxr,nyr] = size(gn.lon_rho);
[nxu,nyu] = size(gn.lon_u);
[nxv,nyv] = size(gn.lon_v);

%% Read data
T =3;% length(time);

% Save 3D variables for interpolation from outer/bigger in temporary arrays
% READ VARIABLES
for mm = 1:T
    zeta_coawst(:,:,mm)   = double(squeeze(nc{'zeta'}(mm,:,:)));
    ubar_coawst(:,:,mm)   = double(squeeze(nc{'ubar'}(mm,:,:)));
    vbar_coawst(:,:,mm)   = double(squeeze(nc{'vbar'}(mm,:,:)));
end

% this is required for calcualtion of water level 
report=0;

%% Initialize the variables of interest
% these are temporary variables created for the interpolation to the
% smaller/refined grid ... 

u = zeros(nxu,nyu,N);
v = zeros(nxv,nyv,N);
temp = zeros(nxr,nyr,N);
salt = zeros(nxr,nyr,N);
ubar = zeros(nxu,nyu,T);
vbar = zeros(nxv,nyv,T);
zeta = zeros(nxr,nyr,T);

% Save the lon and lat of smaller/refined grid for interpolation 
lon_rho_romsz = repmat(gn.lon_rho,1,1,N);
lat_rho_romsz = repmat(gn.lat_rho,1,1,N);


%% READ 3D VARIABLES FOR BRY
% 
for mm = 2:T
    
    aa = squeeze(zeta_coawst(:,:,mm));
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
    
% These are grid(gn.lon_rho) parameters for the outer grid/bigger grid
    aa(~maskr) = NaN;
    aa = maplev(aa);
    zz = griddata(lon_rho(:),lat_rho(:),aa(:),gn.lon_rho,gn.lat_rho);
    zeta(:,:,mm) = zz;
    
    au = squeeze(ubar_coawst(:,:,mm));
    au(~masku) = NaN;
    au = maplev(au);
    %au = au.*coawst_hu;
    aur = u2rho_2d_mai(au);
    
    av = squeeze(vbar_coawst(:,:,mm));
    av(~maskv) = NaN;
    av = maplev(av);
    %av = av.*coawst_hv;
    avr = v2rho_2d_mai(av);
    
    % Compute northward and eastward velocities, important!
    vel = aur + avr.*sqrt(-1);
    vel = vel .* exp(sqrt(-1)*angler);
    velu = real(vel);
    velv = imag(vel);

% These are grid(gn.lon_rho) parameters for the outer grid/bigger grid
    velu1 = griddata(lon_rho(:),lat_rho(:),velu(:),gn.lon_rho,gn.lat_rho);
    velv1 = griddata(lon_rho(:),lat_rho(:),velv(:),gn.lon_rho,gn.lat_rho);
    
    % Rotate velocities to ROMS grid, important!
% These are grid(gn.angle) parameters for the outer grid/bigger grid
    ubar1 = velu1.*cos(gn.angle)+velv1.*sin(gn.angle);
    vbar1 = velv1.*cos(gn.angle)-velu1.*sin(gn.angle);
    ubar(:,:,mm)=rho2u_2d_mw(ubar1);  % defined at u points
    vbar(:,:,mm)=rho2v_2d_mw(vbar1);
    
    % Recalculate ubar and vbar based on refined grid depths.
%     ubar1 = rho2u_2d_mw(ubar1);
%     vbar1 = rho2v_2d_mw(vbar1);
%     ubar1 = ubar1./hu;
%     vbar1 = vbar1./hv;
   
    clear vel au av aur avr velu velv velu1 velv1 ubar1 vbar1
%
% For each 4d variable first save it in 3d from the big grid
% 
%    u_coawst = double(squeeze(nc{'u'}(mm,:,:,:)));
%    v_coawst = double(squeeze(nc{'v'}(mm,:,:,:)));

% These are grid(nc temp) parameters for the outer grid/bigger grid
    temp_coawst=double(squeeze(nc{'temp'}(mm,:,:,:)));
    sand01_coawst=double(squeeze(nc{'sand_01'}(mm,:,:,:)));

    for zz = 1:coawst_N
        aa = squeeze(temp_coawst(zz,:,:));
        aa(~maskr) = NaN;
        aa = maplev(aa);
        temp2(:,:,zz) = aa;
        clear aa
        
        aa = squeeze(sand01_coawst(zz,:,:));
        aa(~maskr) = NaN;
        aa = maplev(aa);
        sand01_tmp_coarse(:,:,zz) = aa;
        clear aa
        
%        aa = squeeze(sand02_coawst(zz,:,:));
%        aa(~maskr) = NaN;
%        aa = maplev(aa);
%        sand02_2(:,:,zz) = aa;
%        clear aa

%        aa = squeeze(sand03_coawst(zz,:,:));
%        aa(~maskr) = NaN;
%        aa = maplev(aa);
%        sand03_2(:,:,zz) = aa;
%        clear aa
%        aa = squeeze(salt_coawst(zz,:,:));
%        aa = squeeze(salt_coawst(zz,:,:));
%        aa(~maskr) = NaN;
%        aa = maplev(aa);
%        salt2(:,:,zz) = aa;
%        clear aa
        
%        aa = squeeze(u_coawst(zz,:,:));
%        aa(~masku) = NaN;
%        aa = maplev(aa);
%        u2(:,:,zz) = aa;
%        clear aa
%
%	aa = squeeze(v_coawst(zz,:,:));
%        aa(~maskv) = NaN;
%       aa = maplev(aa);
%        v2(:,:,zz) = aa;
%        clear aa
    end
%    size(zr)
    for zz = 1:coawst_N
        zr(:,:,zz) = maplev(squeeze(zr(:,:,zz)));
    end
    
% These are grid(nc temp) parameters for the outer grid/bigger grid
%  temp is the inner/refined grid 
    temp          = griddata(lon_rho_z,lat_rho_z,zr(:,:,1:7),temp2,.....
                             lon_rho_romsz,lat_rho_romsz,gn.z_r);
    sand01_tmp_ref= griddata(lon_rho_z,lat_rho_z,zr(:,:,1:7),.....
                             sand01_tmp_coarse,lon_rho_romsz,lat_rho_romsz,gn.z_r);
%   salt = griddata(lon_rho_z,lat_rho_z,zr(:,:,1:7),salt2,lon_rho_romsz,lat_rho_romsz,gn.z_r);
 
% take the "u" velocity that is staggered on rho points.    
%    for zz = 1:coawst_N
%        ur(:,:,zz) = u2rho_2d_mai(u2(:,:,zz));
%        vr(:,:,zz) = v2rho_2d_mai(v2(:,:,zz));
%    end
    
    % Compute Northward and Eastward velocities, angler is the angle for big grid
%    for zz = 1:coawst_N
%        vel = squeeze(ur(:,:,zz))+squeeze(vr(:,:,zz)).*sqrt(-1);
%        vel = vel.* exp(sqrt(-1) * angler);
%        ur(:,:,zz) = real(vel);
%        vr(:,:,zz) = imag(vel);
%    end
    % interpolation 
%    ur2 = griddata(lon_rho_z,lat_rho_z,zr(:,:,1:7),ur,lon_rho_romsz,lat_rho_romsz,gn.z_r);
%    vr2 = griddata(lon_rho_z,lat_rho_z,zr(:,:,1:7),vr,lon_rho_romsz,lat_rho_romsz,gn.z_r);
    
    % Rotate velocities to ROMS grid, important!, gn.angle is angle for refined grid
%    for zz = 1:N
%        ur2(:,:,zz) = squeeze(ur2(:,:,zz)).*cos(gn.angle)+squeeze(vr2(:,:,zz)).*sin(gn.angle);
%        vr2(:,:,zz) = squeeze(vr2(:,:,zz)).*cos(gn.angle)-squeeze(ur2(:,:,zz)).*sin(gn.angle);
%        u(:,:,zz) = rho2u_2d_mw(ur2(:,:,zz));  % defined at u points
%        v(:,:,zz) = rho2v_2d_mw(vr2(:,:,zz));  % defined at v points
%    end
    
    % Remove possible NaN values
    for zz = 1:N
        aa = squeeze(temp(:,:,zz));
        aa = maplev(aa);
        temp(:,:,zz) = aa;
        clear aa
        
        aa = squeeze(sand01_tmp_ref(:,:,zz));
        aa = maplev(aa);
        sand01_tmp_ref(:,:,zz) = aa;
        clear aa
        
%        aa = squeeze(u(:,:,zz));
%        aa = maplev(aa);
%        u(:,:,zz) = aa;
%        clear aa
        
%        aa = squeeze(v(:,:,zz));
%        aa = maplev(aa);
%        v(:,:,zz) = aa;
%        clear aa
    end
% temp_4d is the refined/inner grid
    temp_4d(:,:,:,mm)=temp(:,:,:); 
    sand01_4d(:,:,:,mm)=sand01_tmp_ref(:,:,:); 
end
% 
% %% CREATE BOUNDARY CONDITION
% 
time = time-datenum(1858,11,17);

zeta_north = squeeze(zeta(:,end,:));
ubar_north = squeeze(ubar(:,end,:));
vbar_north = squeeze(vbar(:,end,:));

zeta_south = squeeze(zeta(:,1,:));
ubar_south = squeeze(ubar(:,1,:));
vbar_south = squeeze(vbar(:,1,:));

zeta_east = squeeze(zeta(end,:,:));
ubar_east = squeeze(ubar(end,:,:));
vbar_east = squeeze(vbar(end,:,:));

zeta_west = squeeze(zeta(1,:,:));
ubar_west = squeeze(ubar(1,:,:));
vbar_west = squeeze(vbar(1,:,:));

temp_north = squeeze(temp_4d(:,end,:,:));
temp_south = squeeze(temp_4d(:,1,:,:));
temp_east  = squeeze(temp_4d(end,:,:,:));
temp_west  = squeeze(temp_4d(1,:,:,:));

sand_north_01=squeeze(sand01_4d(:,end,:,:));
sand_south_01=squeeze(sand01_4d(:,1,:,:));
sand_east_01=squeeze(sand01_4d(end,:,:,:));
sand_west_01=squeeze(sand01_4d(1,:,:,:));

create_roms_bry_from_coawst(grd_fname,bry_fname,time,...
    Sinp.theta_s,Sinp.theta_b,Sinp.Tcline,Sinp.Vtransform,Sinp.Vstretching,Sinp.N,...
    zeta_north,ubar_north,vbar_north,...
    ubar_south,vbar_south,zeta_south,...
    zeta_east,ubar_east,vbar_east,...
    zeta_west,ubar_west,vbar_west, ....
    temp_north, temp_south, temp_east, temp_west, ....
    sand_north_01, sand_south_01, sand_east_01, sand_west_01);
