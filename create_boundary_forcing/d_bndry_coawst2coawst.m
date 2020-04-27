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

clear; close all;

% Define filenames.
% model_url : outer grid simulations
%%%%%%%%%%%%%OUTER/COARSE GRID%%%%%%%%%%%%%%%%%%%%%%%%%%
mdl_fname = 'http://geoport.whoi.edu/thredds/dodsC/sand/usgs/users/tkalra/bbleh/bbleh.ncml';
grid_name_outer = '/media/taran/DATADRIVE2/marsh_result/barnegat_bay/all_other_folders/runfiles_bbleh_zd_taran/grid/bbleh_grid_073d.nc';

setup_nctoolbox;
nc = ncgeodataset(mdl_fname);
 time1 = nj_time(nc,'zeta');

time=nc{'ocean_time'}(:); 

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
N       = 10;
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
init_1=1;
init_end=2 ;

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
for mm = init_1:init_end
    zeta_coawst(:,:,mm)   = double(squeeze(nc{'zeta'}(mm,:,:)));
    ubar_coawst(:,:,mm)   = double(squeeze(nc{'ubar'}(mm,:,:)));
    vbar_coawst(:,:,mm)   = double(squeeze(nc{'vbar'}(mm,:,:)));
end

% this is required for calculation of water level 
report=0;

%% Initialize the variables of interest
% these are temporary variables created for the interpolation to the
% smaller/refined grid ... 

%% READ 3D VARIABLES FOR BRY
% 
for mm = init_1:init_end
    
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
    
% interpolate using griddata
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
% These are grid(nc temp) parameters for the outer grid/bigger grid
    temp_coawst=double(squeeze(nc{'temp'}(mm,:,:,:)));
    salt_coawst=double(squeeze(nc{'salt'}(mm,:,:,:)));
    sand01_coawst=double(squeeze(nc{'sand_01'}(mm,:,:,:)));
    sand02_coawst=double(squeeze(nc{'sand_02'}(mm,:,:,:)));
    sand03_coawst=double(squeeze(nc{'sand_03'}(mm,:,:,:)));
    sand04_coawst=double(squeeze(nc{'sand_04'}(mm,:,:,:)));
    sand05_coawst=double(squeeze(nc{'sand_05'}(mm,:,:,:)));
    u_coawst=double(squeeze(nc{'u'}(mm,:,:,:)));
    v_coawst=double(squeeze(nc{'v'}(mm,:,:,:))); 
    
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

        aa = squeeze(sand03_coawst(zz,:,:));
        aa(~maskr) = NaN;
        aa = maplev(aa);
        sand03_tmp_coarse(:,:,zz) = aa;
        clear aa

        aa = squeeze(sand04_coawst(zz,:,:));
        aa(~maskr) = NaN;
        aa = maplev(aa);
        sand04_tmp_coarse(:,:,zz) = aa;
        clear aa

        aa = squeeze(sand05_coawst(zz,:,:));
        aa(~maskr) = NaN;
        aa = maplev(aa);
        sand05_tmp_coarse(:,:,zz) = aa;
        clear aa

        aa = squeeze(temp_coawst(zz,:,:));
        aa(~maskr) = NaN;
        aa = maplev(aa);
        temp_coarse(:,:,zz) = aa;
        clear aa

        aa = squeeze(salt_coawst(zz,:,:));
        aa(~maskr) = NaN;
        aa = maplev(aa);
        salt_coarse(:,:,zz) = aa;
        clear aa ;

        aa = squeeze(u_coawst(zz,:,:));
        aa(~masku) = NaN;
        aa = maplev(aa);
        u2_tmp_coarse(:,:,zz) = aa;
        clear aa

        aa = squeeze(v_coawst(zz,:,:));
        aa(~maskv) = NaN;
        aa = maplev(aa);
        v2_tmp_coarse(:,:,zz) = aa;
        clear aa

    end 
    
    for zz = 1:coawst_N
        zr(:,:,zz) = maplev(squeeze(zr(:,:,zz)));
    end
    
 % Calculate the "u" velocity that is staggered on rho points.
    for zz = 1:coawst_N
        ur_coarse(:,:,zz) = u2rho_2d_mai(u2_tmp_coarse(:,:,zz));
        vr_coarse(:,:,zz) = v2rho_2d_mai(v2_tmp_coarse(:,:,zz));
    end
% Compute Northward and Eastward velocities, angler is the angle for big grid
    for zz = 1:coawst_N
        vel = squeeze(ur_coarse(:,:,zz))+squeeze(vr_coarse(:,:,zz)).*sqrt(-1);
        vel = vel.* exp(sqrt(-1) * angler);
        ur_tmp_coarse(:,:,zz) = real(vel);
        vr_tmp_coarse(:,:,zz) = imag(vel);
    end
    clear ur_coarse vr_coarse u2_tmp_coarse v2_tmp_coarse

% These are grid(nc temp) parameters for the outer grid/bigger grid
%  temp is the inner/refined grid 

%    sand01_tmp_ref= griddata(lon_rho_z,         lat_rho_z,      zr(:,:,1:7),.....
%                             sand01_tmp_coarse, lon_rho_romsz,  lat_rho_romsz, gn.z_r);
    
    % Rotate velocities to ROMS grid, important!, gn.angle is angle for refined grid
    % Remove possible NaN values
%     for zz = 1:N
%         
%         aa = squeeze(sand01_tmp_ref(:,:,zz));
%         aa = maplev(aa);
%         sand01_tmp_ref(:,:,zz) = aa;
%         clear aa
%         
%     end

    % These are grid(nc temp) parameters for the outer grid/bigger grid
%  temp is the inner/refined grid
    temp_ref       = griddata(lon_rho_z,    lat_rho_z,     zr(:,:,1:7), temp_coarse,.....
                             lon_rho_romsz,lat_rho_romsz, gn.z_r);

    salt_ref       = griddata(lon_rho_z,    lat_rho_z,     zr(:,:,1:7), salt_coarse,.....
                             lon_rho_romsz,lat_rho_romsz, gn.z_r);

    sand01_tmp_ref= griddata(lon_rho_z,         lat_rho_z,      zr(:,:,1:7),.....
                             sand01_tmp_coarse, lon_rho_romsz,  lat_rho_romsz, gn.z_r);

    sand02_tmp_ref= griddata(lon_rho_z,         lat_rho_z,      zr(:,:,1:7),.....
                            sand02_tmp_coarse, lon_rho_romsz,  lat_rho_romsz, gn.z_r);

    sand03_tmp_ref= griddata(lon_rho_z,         lat_rho_z,      zr(:,:,1:7),.....
                             sand03_tmp_coarse, lon_rho_romsz,  lat_rho_romsz, gn.z_r);

    sand04_tmp_ref= griddata(lon_rho_z,         lat_rho_z,      zr(:,:,1:7),.....
                             sand04_tmp_coarse, lon_rho_romsz,  lat_rho_romsz, gn.z_r);

    sand05_tmp_ref= griddata(lon_rho_z,         lat_rho_z,      zr(:,:,1:7),.....
                             sand05_tmp_coarse, lon_rho_romsz,  lat_rho_romsz, gn.z_r);

    % interpolation
    ur_ref = griddata(lon_rho_z,      lat_rho_z, zr(:,:,1:7),    ur_tmp_coarse,.....
                      lon_rho_romsz,  lat_rho_romsz,gn.z_r);
    vr_ref = griddata(lon_rho_z,      lat_rho_z,    zr(:,:,1:7), vr_tmp_coarse,.....
                       lon_rho_romsz, lat_rho_romsz,gn.z_r);

    % Rotate velocities to ROMS grid, important!, gn.angle is angle for refined grid
    for zz = 1:N
        ur2_ref(:,:,zz) = squeeze(ur_ref(:,:,zz)).*cos(gn.angle) + .....
                          squeeze(vr_ref(:,:,zz)).*sin(gn.angle);
        vr2_ref(:,:,zz) = squeeze(vr_ref(:,:,zz)).*cos(gn.angle) -.....
                          squeeze(ur_ref(:,:,zz)).*sin(gn.angle);
        u_tmp_ref(:,:,zz) = rho2u_2d_mw(ur2_ref(:,:,zz));  % defined at u points
        v_tmp_ref(:,:,zz) = rho2v_2d_mw(vr2_ref(:,:,zz));  % defined at v points
    end
    clear ur_ref vr_ref ur2_ref  vr2_ref

       % Remove possible NaN values
    for zz = 1:N
        aa = squeeze(temp_ref(:,:,zz));
        aa = maplev(aa);
        temp_ref(:,:,zz) = aa;
        clear aa

        aa = squeeze(salt_ref(:,:,zz));
        aa = maplev(aa);
        salt_ref(:,:,zz) = aa;
        clear aa

        aa = squeeze(sand01_tmp_ref(:,:,zz));
        aa = maplev(aa);
        sand01_tmp_ref(:,:,zz) = aa;
        clear aa

        aa = squeeze(sand02_tmp_ref(:,:,zz));
        aa = maplev(aa);
        sand02_tmp_ref(:,:,zz) = aa;
        clear aa

        aa = squeeze(sand03_tmp_ref(:,:,zz));
        aa = maplev(aa);
        sand03_tmp_ref(:,:,zz) = aa;
        clear aa

        aa = squeeze(sand04_tmp_ref(:,:,zz));
        aa = maplev(aa);
        sand04_tmp_ref(:,:,zz) = aa;
        clear aa

        aa = squeeze(sand05_tmp_ref(:,:,zz));
        aa = maplev(aa);
        sand05_tmp_ref(:,:,zz) = aa;
       clear aa

        aa = squeeze(u_tmp_ref(:,:,zz));
        aa = maplev(aa);
        u_tmp_ref(:,:,zz) = aa;
        clear aa

        aa = squeeze(v_tmp_ref(:,:,zz));
        aa = maplev(aa);
        v_tmp_ref(:,:,zz) = aa;
        clear aa
    end
% temp_4d is the refined/inner grid
    temp_4d(:,:,:,mm)=temp_ref(:,:,:);
    salt_4d(:,:,:,mm)=salt_ref(:,:,:);
    sand01_4d(:,:,:,mm)=sand01_tmp_ref(:,:,:);
    sand02_4d(:,:,:,mm)=sand02_tmp_ref(:,:,:);
    sand03_4d(:,:,:,mm)=sand03_tmp_ref(:,:,:);
    sand04_4d(:,:,:,mm)=sand04_tmp_ref(:,:,:);
    sand05_4d(:,:,:,mm)=sand05_tmp_ref(:,:,:);
    u_4d(:,:,:,mm)=u_tmp_ref(:,:,:) ;
    v_4d(:,:,:,mm)=v_tmp_ref(:,:,:) ;
end
clear temp_ref salt_ref sand01_tmp_ref sand02_tmp_ref
clear sand03_tmp_ref sand04_temp_ref sand05_tmp_ref
clear u_tmp_ref v_tmp_ref


% %% CREATE BOUNDARY CONDITION
% 
time=time(init_1:init_end); 

zeta_north = squeeze(zeta(:, end, init_1:init_end));
zeta_south = squeeze(zeta(:,  1,  init_1:init_end));
zeta_east = squeeze(zeta(end, :,  init_1:init_end));
zeta_west = squeeze(zeta(1,   :,  init_1:init_end));

ubar_north = squeeze(ubar(:, end, init_1:init_end));
ubar_south = squeeze(ubar(:,  1,  init_1:init_end));
ubar_east = squeeze(ubar(end, :,  init_1:init_end));
ubar_west = squeeze(ubar(1,   :,  init_1:init_end));

vbar_north = squeeze(vbar(:, end, init_1:init_end));
vbar_south = squeeze(vbar(:,  1,  init_1:init_end));
vbar_east = squeeze(vbar(end, :,  init_1:init_end));
vbar_west = squeeze(vbar(1,   :,  init_1:init_end));

temp_north = squeeze(temp_4d(:, end, :,init_1:init_end));
temp_south = squeeze(temp_4d(:,  1,  :,init_1:init_end))
temp_east  = squeeze(temp_4d(end,:,  :,init_1:init_end));
temp_west  = squeeze(temp_4d(1,  :,  :,init_1:init_end));

salt_north = squeeze(salt_4d(:, end, :,init_1:init_end));
salt_south = squeeze(salt_4d(:,  1,  :,init_1:init_end));
salt_east  = squeeze(salt_4d(end,:,  :,init_1:init_end));
salt_west  = squeeze(salt_4d(1,  :,  :,init_1:init_end));

u_north=squeeze(u_4d(:, end,:, init_1:init_end));
u_south=squeeze(u_4d(:,  1, :, init_1:init_end));
u_east=squeeze (u_4d(end, :, :, init_1:init_end));
u_west=squeeze (u_4d(1,   :, :, init_1:init_end));

v_north=squeeze(v_4d(:, end,:, init_1:init_end));
v_south=squeeze(v_4d(:,  1, :, init_1:init_end));
v_east=squeeze (v_4d(end, :, :, init_1:init_end));
v_west=squeeze (v_4d(1,   :, :, init_1:init_end));

sand_north_01=squeeze(sand01_4d(:,end,:,:));
sand_south_01=squeeze(sand01_4d(:,1,:,:));
sand_east_01=squeeze(sand01_4d(end,:,:,:));
sand_west_01=squeeze(sand01_4d(1,:,:,:));

sand_north_02=squeeze(sand02_4d(:,end,:,:));
sand_south_02=squeeze(sand02_4d(:,1,:,:));
sand_east_02=squeeze(sand02_4d(end,:,:,:));
sand_west_02=squeeze(sand02_4d(1,:,:,:));

sand_north_03=squeeze(sand03_4d(:,end,:,:));
sand_south_03=squeeze(sand03_4d(:,1,:,:));
sand_east_03=squeeze(sand03_4d(end,:,:,:));
sand_west_03=squeeze(sand03_4d(1,:,:,:));

sand_north_04=squeeze(sand04_4d(:,end,:,:));
sand_south_04=squeeze(sand04_4d(:,1,:,:));
sand_east_04=squeeze(sand04_4d(end,:,:,:));
sand_west_04=squeeze(sand04_4d(1,:,:,:));

sand_north_05=squeeze(sand05_4d(:,end,:,:));
sand_south_05=squeeze(sand05_4d(:,1,:,:));
sand_east_05=squeeze(sand05_4d(end,:,:,:));
sand_west_05=squeeze(sand05_4d(1,:,:,:));

create_roms_bry_from_coawst(grid_name_inner, bry_fname, time,...
    Sinp.theta_s,Sinp.theta_b,Sinp.Tcline,Sinp.Vtransform,Sinp.Vstretching,Sinp.N,...
    zeta_north,ubar_north,vbar_north,...
    ubar_south,vbar_south,zeta_south,...
    zeta_east,ubar_east,vbar_east,...
    zeta_west,ubar_west,vbar_west, ....
    temp_north, temp_south, temp_east, temp_west, ....
    salt_north, salt_south, salt_east, salt_west, ....
    u_north, u_south, u_east, u_west, .....
    v_north, v_south, v_east, v_west,......
    sand_north_01, sand_south_01, sand_east_01, sand_west_01,....
    sand_north_02, sand_south_02, sand_east_02, sand_west_02,....
    sand_north_03, sand_south_03, sand_east_03, sand_west_03,....
    sand_north_04, sand_south_04, sand_east_04, sand_west_04,....
    sand_north_05, sand_south_05, sand_east_05, sand_west_05);

