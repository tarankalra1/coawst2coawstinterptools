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
% CAH edits

clear; close all;

% Define filenames.
% model_url : outer grid simulations
%mdl_fname = 'http://geoport.whoi.edu/thredds/dodsC/sand/usgs/users/chegermiller/projects/matthew_coamps_SAB_TY_SSW_relwind_7218/Output_nest/roms_his.ncml';
mdl_fname = 'http://geoport.whoi.edu/thredds/dodsC/sand/usgs/users/tkalra/bbleh/bbleh.ncml';
grid_name_outer = '/media/taran/DATADRIVE2/marsh_result/barnegat_bay/all_other_folders/runfiles_bbleh_zd_taran/grid/bbleh_grid_073d.nc';

% gridname : refined grid
% grd_fname  = '/Users/Christie/Desktop/Research/USGS/iFMSIP/modeling/grids/GTM_COARSE_May12018_smooth_open_s.nc';
%grd_fname = '/media/taran/DATADRIVE2/marsh_result/barnegat_bay/all_other_folders/runfiles_bbleh_zd_taran/grid/bbleh_grid_073d.nc';
grd_fname = '/media/taran/DATADRIVE2/marsh_result/barnegat_bay/easygrid/bbleh_reedy_grd.nc';

% init_file : filename to create for initial conditions
%init_fname = '/Users/Christie/Desktop/Research/USGS/iFMSIP/modeling/prep_runs/data/GTM_COARSE_SAB_TY_SSW_relwind_7218_init.nc';
init_fname = '/media/taran/DATADRIVE2/marsh_result/barnegat_bay/all_other_folders/runfiles_bbleh_zd_taran/bbleh_init_0126.nc';
init_fname=  'reedy_init.nc'
% bry_file : filename to create for boundary conditions
%bry_fname  = '/Users/Christie/Desktop/Research/USGS/iFMSIP/modeling/prep_runs/data/GTM_COARSE_SAB_TY_SSW_relwind_7218_bry.nc';
bry_fname=  'reedy_bry.nc'
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

Sinp.hc = hc;                    % stretching width used in ROMS
gn = get_roms_grid(grd_fname,Sinp);
[nxr,nyr] = size(gn.lon_rho);
[nxu,nyu] = size(gn.lon_u);
[nxv,nyv] = size(gn.lon_v);

%% Read data

setup_nctoolbox;
nc = ncgeodataset(mdl_fname);
time = nj_time(nc,'zeta');

T =100;% length(time);

% GET HORIZONTAL GRID
g = nc{'zeta'}(1,:,:).grid;
lon_rho = g.lon;
lat_rho = g.lat;

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

% % TSK 
% coawst_N=coawst_N+1 ; 

lon_rho_z = repmat(lon_rho,1,1,coawst_N);
lat_rho_z = repmat(lat_rho,1,1,coawst_N);

% READ VARIABLES
for mm = 1:T
    zeta_coawst(:,:,mm)   = double(squeeze(nc{'zeta'}(mm,:,:)));
    ubar_coawst(:,:,mm)   = double(squeeze(nc{'ubar'}(mm,:,:)));
    vbar_coawst(:,:,mm)   = double(squeeze(nc{'vbar'}(mm,:,:)));
end

u_coawst    = double(squeeze(nc{'u'}(1,:,:,:)));
v_coawst    = double(squeeze(nc{'v'}(1,:,:,:)));
temp_coawst = double(squeeze(nc{'temp'}(1,:,:,:)));
salt_coawst = double(squeeze(nc{'salt'}(1,:,:,:)));

save coawst_vars.mat zeta_coawst ubar_coawst vbar_coawst u_coawst v_coawst temp_coawst salt_coawst -v7.3

close(nc);

%% Initialize the variables of interest

u = zeros(nxu,nyu,N);
v = zeros(nxv,nyv,N);
temp = zeros(nxr,nyr,N);
salt = zeros(nxr,nyr,N);
ubar = zeros(nxu,nyu,T);
vbar = zeros(nxv,nyv,T);
zeta = zeros(nxr,nyr,T);

lon_rho_romsz = repmat(gn.lon_rho,1,1,N);
lat_rho_romsz = repmat(gn.lat_rho,1,1,N);

%% Compute depth at u, v points for both grids and get rid of land or masked cells
% Depending on the setup and the outer/refined grids, this next section may
% or may not be necessary. In the next section, I have commented out the
% lines that refer to the variables created here, coawst_hu and coawst_hv.

% hu = rho2u_2d_mw(h);
% hu(~gn.mask_u) = NaN;
% hu(hu<=0) = NaN;
% hv = rho2v_2d_mw(h);
% hv(~gn.mask_v) = NaN;
% hv(hv<=0) = NaN;
%  
% coawst_hu = rho2u_2d_mai(coawst_h);
% coawst_hu(~masku) = NaN;
% coawst_hu(coawst_hu<=0) = NaN;
% coawst_hv = rho2v_2d_mai(coawst_h);
% coawst_hv(~maskv) = NaN;
% coawst_hv(coawst_hv<=0) = NaN;

%% READ 3D & 4D VARIABLES FOR INIT
%Vtransform=1; Vstretching =1; theta_s=0.0; theta_b=0.0;
%hc=0.0; igrid=5 ;
% 
 report=0;
% for t=tstart:tend
%  z(:,:,:,t)=set_depth(Vtransform, Vstretching, ...
%                       theta_s, theta_b, hc, N, ...
%                       igrid, h, squeeze(zeta(:,:,t)), report);
iinit_file=0;  %if init_file is not required 
if(iinit_file==1); 
 for mm = 1
    aa = sq(zeta_coawst(:,:,mm));
        
    % Compute vertical elevations of the grid, this is time dependent
    % These depths are negative
    if mod(mm-1,24) == 0
        zr = set_depth(coawst_Vtransform, coawst_Vstretching, ....
                       coawst_theta_s,    coawst_theta_b, .....
                       coawst_hc,               coawst_N,.........
                       5,                 coawst_h, aa, report);
    end
    size(zr)
     
    aa(~maskr) = NaN;
    aa = maplev(aa);
    zz = griddata(lon_rho(:),lat_rho(:),aa(:),gn.lon_rho,gn.lat_rho);
    zeta(:,:,mm) = zz;
    
    au = sq(ubar_coawst(:,:,mm));
    au(~masku) = NaN;
    au = maplev(au);
    % au = au.*coawst_hu;
    aur = u2rho_2d_mai(au);
    
    av = sq(vbar_coawst(:,:,mm));
    av(~maskv) = NaN;
    av = maplev(av);
    % av = av.*coawst_hv;
    avr = v2rho_2d_mai(av);
    
    % Compute northward and eastward velocities, important!
    vel = aur + avr.*sqrt(-1);
    vel = vel .* exp(sqrt(-1)*angler);
    velu = real(vel);
    velv = imag(vel);
    % Projecting 3-D vars
    velu1 = griddata(lon_rho(:),lat_rho(:),velu(:),gn.lon_rho,gn.lat_rho);
    velv1 = griddata(lon_rho(:),lat_rho(:),velv(:),gn.lon_rho,gn.lat_rho);
    
    % Rotate velocities to ROMS grid, important!
    ubar1 = velu1.*cos(gn.angle)+velv1.*sin(gn.angle);
    vbar1 = velv1.*cos(gn.angle)-velu1.*sin(gn.angle);
    ubar(:,:,mm)=rho2u_2d_mw(ubar1);  % defined at u points
    vbar(:,:,mm)=rho2v_2d_mw(vbar1); 
    
    % Recalculate ubar and vbar based on refined grid depths.
%     ubar1 = rho2u_2d_mw(ubar1);
%     vbar1 = rho2v_2d_mw(vbar1);
%     ubar1 = ubar1./hu;
%     vbar1 = vbar1./hv;
%     ubar(:,:,mm) = maplev(ubar1);
%     vbar(:,:,mm) = maplev(vbar1);
    
    clear vel au av aur avr velu velv velu1 velv1 ubar1 vbar1
    
    for zz = 1:coawst_N
        aa = squeeze(temp_coawst(zz,:,:));
        aa(~maskr) = NaN;
        aa = maplev(aa);
        temp2(:,:,zz) = aa;
        clear aa
        
        aa = squeeze(salt_coawst(zz,:,:));
        aa(~maskr) = NaN;
        aa = maplev(aa);
        salt2(:,:,zz) = aa;
        clear aa
        
        aa = squeeze(u_coawst(zz,:,:));
        aa(~masku) = NaN;
        aa = maplev(aa);
        u2(:,:,zz) = aa;
        clear aa
        
        aa = squeeze(v_coawst(zz,:,:));
        aa(~maskv) = NaN;
        aa = maplev(aa);
        v2(:,:,zz) = aa;
        clear aa
    end
    size(zr)
    for zz = 1:coawst_N
        zr(:,:,zz) = maplev(squeeze(zr(:,:,zz)));
    end
    
   temp = griddata(lon_rho_z,lat_rho_z,zr(:,:,1:7),temp2,lon_rho_romsz,lat_rho_romsz,gn.z_r);
   salt = griddata(lon_rho_z,lat_rho_z,zr(:,:,1:7),salt2,lon_rho_romsz,lat_rho_romsz,gn.z_r);
    
    for zz = 1:coawst_N
        ur(:,:,zz) = u2rho_2d_mai(u2(:,:,zz));
        vr(:,:,zz) = v2rho_2d_mai(v2(:,:,zz));
    end
    
    % Compute Northward and Eastward velocities
    for zz = 1:coawst_N
        vel = squeeze(ur(:,:,zz))+squeeze(vr(:,:,zz)).*sqrt(-1);
        vel = vel.* exp(sqrt(-1) * angler);
        ur(:,:,zz) = real(vel);
        vr(:,:,zz) = imag(vel);
    end
    
    ur2 = griddata(lon_rho_z,lat_rho_z,zr(:,:,1:7),ur,lon_rho_romsz,lat_rho_romsz,gn.z_r);
    vr2 = griddata(lon_rho_z,lat_rho_z,zr(:,:,1:7),vr,lon_rho_romsz,lat_rho_romsz,gn.z_r);
    
    % Rotate velocities to ROMS grid, important!
    for zz = 1:N
        ur2(:,:,zz) = squeeze(ur2(:,:,zz)).*cos(gn.angle)+squeeze(vr2(:,:,zz)).*sin(gn.angle);
        vr2(:,:,zz) = squeeze(vr2(:,:,zz)).*cos(gn.angle)-squeeze(ur2(:,:,zz)).*sin(gn.angle);
        u(:,:,zz) = rho2u_2d_mw(ur2(:,:,zz));  % defined at u points
        v(:,:,zz) = rho2v_2d_mw(vr2(:,:,zz));  % defined at v points
    end
    
    % Remove possible NaN values
    for zz = 1:N
        aa = squeeze(temp(:,:,zz));
        aa = maplev(aa);
        temp(:,:,zz) = aa;
        clear aa
        
        aa = squeeze(salt(:,:,zz));
        aa = maplev(aa);
        salt(:,:,zz) = aa;
        clear aa
        
        aa = squeeze(u(:,:,zz));
        aa = maplev(aa);
        u(:,:,zz) = aa;
        clear aa
        
        aa = squeeze(v(:,:,zz));
        aa = maplev(aa);
        v(:,:,zz) = aa;
        clear aa
    end
    
    %% CREATE INITIAL CONDITION
    init_time = time(1,1)-datenum(1858,11,17);
    create_roms_init_from_coawst(grd_fname,init_fname,init_time,...
        Sinp.theta_s,Sinp.theta_b,Sinp.Tcline,Sinp.Vtransform,Sinp.Vstretching,...
        Sinp.N,u,v,squeeze(ubar(:,:,1)),squeeze(vbar(:,:,1)),...
        temp,salt,squeeze(zeta(:,:,1)));
  end
end 
%% READ 3D VARIABLES FOR BRY
% 
for mm = 2:T
    
    aa = squeeze(zeta_coawst(:,:,mm));
    size(aa)
    % Compute vertical elevations of the grid, this is time dependent
    % These depths are negative
    if mod(mm-1,24) == 0;
         zr = set_depth(coawst_Vtransform, coawst_Vstretching, ....
                       coawst_theta_s,    coawst_theta_b, .....
                       coawst_hc,               coawst_N,.........
                       5,                 coawst_h, aa, report);
                   
%         zr = set_depth(coawst_Vtransform,coawst_Vstretching,coawst_theta_s,.......
%                        coawst_theta_b,coawst_hc,coawst_N, ...
%                        1,coawst_h,aa);
    end
    
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
    velu1 = griddata(lon_rho(:),lat_rho(:),velu(:),gn.lon_rho,gn.lat_rho);
    velv1 = griddata(lon_rho(:),lat_rho(:),velv(:),gn.lon_rho,gn.lat_rho);
    
    % Rotate velocities to ROMS grid, important!
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

create_roms_bry_from_coawst(grd_fname,bry_fname,time,...
    Sinp.theta_s,Sinp.theta_b,Sinp.Tcline,Sinp.Vtransform,Sinp.Vstretching,Sinp.N,...
    zeta_north,ubar_north,vbar_north,...
    ubar_south,vbar_south,zeta_south,...
    zeta_east,ubar_east,vbar_east,...
    zeta_west,ubar_west,vbar_west);
