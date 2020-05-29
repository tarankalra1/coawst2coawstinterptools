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

%[nxu,nyu] = size(gn.lon_u);
%[nxv,nyv] = size(gn.lon_v);

%%%%%%%%%%%%%%%%%Choose time steps for forcing file %%%%%%%%%%%%%%%%%%%%%%%5
%T =3;% length(time);
 
%u = zeros(nxu,nyu,N);
%v = zeros(nxv,nyv,N);
%temp = zeros(nxr,nyr,N);
%salt = zeros(nxr,nyr,N);
%ubar = zeros(nxu,nyu,init_end);
%vbar = zeros(nxv,nyv,init_end);
%zeta = zeros(nxr,nyr,init_end);

% Save the lon and lat of smaller/refined grid for interpolation 
%lon_rho_romsz = repmat(gn.lon_rho,1,1,N);
%lat_rho_romsz = repmat(gn.lat_rho,1,1,N);

% Save 3D variables for interpolation from outer/bigger in temporary arrays
% READ VARIABLES

% this is required for calculation of water level 
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
    [zeta_ref]=interp3d_tsk(lon_rho_coarse_col, lat_rho_coarse_col,...
                            zeta_coarse, lon_rho_ref_col, lat_rho_ref_col, ....
                            nx_ref, ny_ref);
  
    zeta_ref(:,:,mm-init_1+1)=zeta_ref; 
    
    ubar_coarse(:,:)=double(squeeze(nc{'ubar'}(mm,:,:)));
    au = ubar_coarse;
    au(~coarse_masku) = NaN;
    au = maplev(au);
    %au = au.*coawst_hu;
    aur = u2rho_2d_mai(au);

    vbar_coarse(:,:)=double(squeeze(nc{'vbar'}(mm,:,:)));
    av = vbar_coarse;
    av(~coarse_maskv) = NaN;
    av = maplev(av);
    %av = av.*coawst_hv;
    avr = v2rho_2d_mai(av);

    % Compute northward and eastward velocities, important!
    vel = aur + avr.*sqrt(-1);
    vel = vel .* exp(sqrt(-1)*angler);
    velu = real(vel);
    velv = imag(vel);

% These are grid(gn.lon_rho) parameters for the outer grid/bigger grid
%    velu1 = griddata(lon_rho(:),lat_rho(:),velu(:),gn.lon_rho,gn.lat_rho);
%    velv1 = griddata(lon_rho(:),lat_rho(:),velv(:),gn.lon_rho,gn.lat_rho);
    [velu1]=interp3d_tsk(lon_rho_coarse_col, lat_rho_coarse_col,...
                            velu, lon_rho_ref_col, lat_rho_ref_col, ....
                            nx_ref, ny_ref);
    [velv1]=interp3d_tsk(lon_rho_coarse_col, lat_rho_coarse_col,...
                            velv, lon_rho_ref_col, lat_rho_ref_col, ....
                            nx_ref, ny_ref);

% Rotate velocities to ROMS grid, important!
% These are grid(gn.angle) parameters for the outer grid/bigger grid
    ubar1 = velu1.*cos(gn.angle)+velv1.*sin(gn.angle);
    vbar1 = velv1.*cos(gn.angle)-velu1.*sin(gn.angle);

    mm
    mm-init_1+1
% This will be ubar, vbar for refined mesh
    ubar_ref(:,:,mm-init_1+1)=rho2u_2d_mw(ubar1);  % defined at u points
    vbar_ref(:,:,mm-init_1+1)=rho2v_2d_mw(vbar1);

   clear zeta zeta_coarse ubar_coarse vbar_coarse
   clear vel au av aur avr velu velv velu1 velv1 ubar1 vbar1
% 
% PROCEED TO 4D variables
% These are grid(nc temp) parameters for the outer grid/bigger grid
    temp_coarse=double(squeeze(nc{'temp'}(mm,:,:,:)));
    salt_coarse=double(squeeze(nc{'salt'}(mm,:,:,:)));
    sand01_coarse=double(squeeze(nc{'sand_01'}(mm,:,:,:)));
    sand02_coarse=double(squeeze(nc{'sand_02'}(mm,:,:,:)));
    sand03_coarse=double(squeeze(nc{'sand_03'}(mm,:,:,:)));
    sand04_coarse=double(squeeze(nc{'sand_04'}(mm,:,:,:)));
    sand05_coarse=double(squeeze(nc{'sand_05'}(mm,:,:,:)));
    u_coarse=double(squeeze(nc{'u'}(mm,:,:,:)));
    v_coarse=double(squeeze(nc{'v'}(mm,:,:,:)));

% % 4d vars 
   for zz = 1:coarse_N
      aa = squeeze(sand01_coarse(zz,:,:));
      aa(~coarse_maskr) = NaN;
      aa= maplev(aa);
      sand01_tmp_coarse(:,:,zz)=aa;
      clear aa
      
      aa = squeeze(sand02_coarse(zz,:,:));
      aa(~coarse_maskr) = NaN;
      aa= maplev(aa);
      sand02_tmp_coarse(:,:,zz)=aa;
      clear aa

      aa = squeeze(sand03_coarse(zz,:,:));
      aa(~coarse_maskr) = NaN;
      aa= maplev(aa);
      sand03_tmp_coarse(:,:,zz)=aa;
      clear aa

      aa = squeeze(sand04_coarse(zz,:,:));
      aa(~coarse_maskr) = NaN;
      aa= maplev(aa);
      sand04_tmp_coarse(:,:,zz)=aa;
      clear aa

      aa = squeeze(sand05_coarse(zz,:,:));
      aa(~coarse_maskr) = NaN;
      aa= maplev(aa);
      sand05_tmp_coarse(:,:,zz)=aa;
      clear aa

      aa = squeeze(temp_coarse(zz,:,:));
      aa(~coarse_maskr) = NaN;
      aa = maplev(aa);
      temp_tmp_coarse(:,:,zz) = aa;
      clear aa 
%       
      aa = squeeze(salt_coarse(zz,:,:));
      aa(~coarse_maskr) = NaN;
      aa = maplev(aa);
      salt_tmp_coarse(:,:,zz) = aa;
      clear aa ;

      aa = squeeze(u_coarse(zz,:,:));
      aa(~coarse_masku) = NaN;
      aa = maplev(aa);
      u2_tmp_coarse(:,:,zz) = aa;
      clear aa

      aa = squeeze(v_coarse(zz,:,:));
      aa(~coarse_maskv) = NaN;
      aa = maplev(aa);
      v2_tmp_coarse(:,:,zz) = aa;
      clear aa

  end 
%     
  for zz = 1:coarse_N
      zr_coarse(:,:,zz) = maplev(squeeze(zr_coarse(:,:,zz)));
  end
 
 % Calculate the "u" velocity that is staggered on rho points.
    for zz = 1:coarse_N
        ur_coarse(:,:,zz) = u2rho_2d_mai(u2_tmp_coarse(:,:,zz));
        vr_coarse(:,:,zz) = v2rho_2d_mai(v2_tmp_coarse(:,:,zz));
    end
% Compute Northward and Eastward velocities, angler is the angle for big grid
    for zz = 1:coarse_N
        vel = squeeze(ur_coarse(:,:,zz))+squeeze(vr_coarse(:,:,zz)).*sqrt(-1);
        vel = vel.* exp(sqrt(-1) * angler);
        ur_tmp_coarse(:,:,zz) = real(vel);
        vr_tmp_coarse(:,:,zz) = imag(vel);
    end

    clear ur_coarse vr_coarse u2_tmp_coarse v2_tmp_coarse
    clear temp_coarse salt_coarse
    clear sand01_coarse sand02_coarse sand03_coarse
    clear sand04_coarse sand05_coarse
    clear u_coarse v_coarse

%     
%  % Calculate the "u" velocity that is staggered on rho points.

   zr_coarse_3d_col=reshape( zr_coarse(:,:,1:coarse_N), ...
                             [grid_size_coarse*coarse_N 1] ) ; 

   zr_ref_3d_col   =reshape( zr_ref(:,:,1:ref_N), ....
                             [grid_size_ref*ref_N 1] ); 

   [sand01_tmp_ref]=interp4d_tsk(lon_rho_coarse_3d_col, lat_rho_coarse_3d_col,...
                                  zr_coarse_3d_col, ...
                                  sand01_tmp_coarse, .....
                                  lon_rho_ref_3d_col, lat_rho_ref_3d_col, ....
                                  zr_ref_3d_col,....
                                  nx_ref, ny_ref, ref_N);
	
   [sand02_tmp_ref]=interp4d_tsk(lon_rho_coarse_3d_col, lat_rho_coarse_3d_col,...
                                  zr_coarse_3d_col, ...
                                  sand02_tmp_coarse, .....
                                  lon_rho_ref_3d_col, lat_rho_ref_3d_col, ....
                                  zr_ref_3d_col,....
                                  nx_ref, ny_ref, ref_N);		    
     
   [sand03_tmp_ref]=interp4d_tsk(lon_rho_coarse_3d_col, lat_rho_coarse_3d_col,...
                                  zr_coarse_3d_col, ...
                                  sand03_tmp_coarse, .....
                                  lon_rho_ref_3d_col, lat_rho_ref_3d_col, ....
                                  zr_ref_3d_col,....
                                  nx_ref, ny_ref, ref_N);

   [sand04_tmp_ref]=interp4d_tsk(lon_rho_coarse_3d_col, lat_rho_coarse_3d_col,...
                                  zr_coarse_3d_col, ...
                                  sand04_tmp_coarse, .....
                                  lon_rho_ref_3d_col, lat_rho_ref_3d_col, ....
                                  zr_ref_3d_col,....
                                  nx_ref, ny_ref, ref_N);

   [sand05_tmp_ref]=interp4d_tsk(lon_rho_coarse_3d_col, lat_rho_coarse_3d_col,...
                                  zr_coarse_3d_col, ...
                                  sand05_tmp_coarse, .....
                                  lon_rho_ref_3d_col, lat_rho_ref_3d_col, ....
                                  zr_ref_3d_col,....
                                  nx_ref, ny_ref, ref_N);

   [temp_tmp_ref]=interp4d_tsk(lon_rho_coarse_3d_col, lat_rho_coarse_3d_col,...
                                  zr_coarse_3d_col, ...
                                  temp_tmp_coarse, .....
                                  lon_rho_ref_3d_col, lat_rho_ref_3d_col, ....
                                  zr_ref_3d_col,....
                                  nx_ref, ny_ref, ref_N);

   [salt_tmp_ref]=interp4d_tsk(lon_rho_coarse_3d_col, lat_rho_coarse_3d_col,...
                                  zr_coarse_3d_col, ...
                                  salt_tmp_coarse, .....
                                  lon_rho_ref_3d_col, lat_rho_ref_3d_col, ....
                                  zr_ref_3d_col,....
                                  nx_ref, ny_ref, ref_N);		 
%
   [ur_ref]=interp4d_tsk(lon_rho_coarse_3d_col, lat_rho_coarse_3d_col,...
                                  zr_coarse_3d_col, ...
                                  ur_tmp_coarse, .....
                                  lon_rho_ref_3d_col, lat_rho_ref_3d_col, ....
                                  zr_ref_3d_col,....
                                  nx_ref, ny_ref, ref_N);		 

   [vr_ref]=interp4d_tsk(lon_rho_coarse_3d_col, lat_rho_coarse_3d_col,...
                                  zr_coarse_3d_col, ...
                                  vr_tmp_coarse, .....
                                  lon_rho_ref_3d_col, lat_rho_ref_3d_col, ....
                                  zr_ref_3d_col,....
                                  nx_ref, ny_ref, ref_N);		 


    % Rotate velocities to ROMS grid, important!, gn.angle is angle for refined grid
    for zz = 1:ref_N
        ur2_ref(:,:,zz) = squeeze(ur_ref(:,:,zz)).*cos(gn.angle) + .....
                          squeeze(vr_ref(:,:,zz)).*sin(gn.angle);
        vr2_ref(:,:,zz) = squeeze(vr_ref(:,:,zz)).*cos(gn.angle) -.....
                          squeeze(ur_ref(:,:,zz)).*sin(gn.angle);
        u_tmp_ref(:,:,zz) = rho2u_2d_mw(ur2_ref(:,:,zz));  % defined at u points
        v_tmp_ref(:,:,zz) = rho2v_2d_mw(vr2_ref(:,:,zz));  % defined at v points
    end
    
    clear ur_ref vr_ref ur2_ref vr2_ref
    clear temp_tmp_coarse salt_tmp_coarse
    clear sand01_tmp_coarse sand02_tmp_coarse sand03_tmp_coarse 
    clear sand04_tmp_coarse sand05_tmp_coarse
    clear ur_tmp_coarse vr_tmp_coarse

       % Remove possible NaN values
%     for zz = 1:ref_N
%         aa = squeeze(sand01_tmp_ref(:,:,zz));
%         aa = maplev(aa);
%         sand01_tmp_ref(:,:,zz) = aa;
%         clear aa
%     end

% temp_4d is the refined/inner grid
% temp_4d is the refined/inner grid
    temp_4d(:,:,:,mm-init_1+1)=temp_tmp_ref(:,:,:);
    salt_4d(:,:,:,mm-init_1+1)=salt_tmp_ref(:,:,:);
    sand01_4d(:,:,:,mm-init_1+1)=sand01_tmp_ref(:,:,:);
    sand02_4d(:,:,:,mm-init_1+1)=sand02_tmp_ref(:,:,:);
    sand03_4d(:,:,:,mm-init_1+1)=sand03_tmp_ref(:,:,:);
    sand04_4d(:,:,:,mm-init_1+1)=sand04_tmp_ref(:,:,:);
    sand05_4d(:,:,:,mm-init_1+1)=sand05_tmp_ref(:,:,:);
    u_4d(:,:,:,mm-init_1+1)=u_tmp_ref(:,:,:) ;
    v_4d(:,:,:,mm-init_1+1)=v_tmp_ref(:,:,:) ;
end
clear temp_ref salt_ref sand01_tmp_ref sand02_tmp_ref
clear sand03_tmp_ref sand04_temp_ref sand05_tmp_ref
clear ur_tmp_ref vr_tmp_ref

toc    
%clear temp_ref salt_ref sand01_tmp_ref sand02_tmp_ref
%clear sand03_tmp_ref sand04_temp_ref sand05_tmp_ref
%clear u_tmp_ref v_tmp_ref

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
  zeta_north = squeeze(zeta_ref(:, end,      count_1:count_end));
  ubar_north = squeeze(ubar_ref(:, end,      count_1:count_end));
  vbar_north = squeeze(vbar_ref(:, end,      count_1:count_end));
  
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
  zeta_south = squeeze(zeta_ref(:,  1,        count_1:count_end));
  ubar_south = squeeze(ubar_ref(:,  1,        count_1:count_end));
  vbar_south = squeeze(vbar_ref(:,  1,        count_1:count_end));
  
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
  zeta_east = squeeze(zeta_ref(end, :,       count_1:count_end));
  ubar_east = squeeze(ubar_ref(end, :,       count_1:count_end));
  vbar_east = squeeze(vbar_ref(end, :,       count_1:count_end));
  
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
  zeta_west = squeeze(zeta_ref(1,   :,       count_1:count_end));
  ubar_west = squeeze(ubar_ref(1,   :,       count_1:count_end)); 
  vbar_west = squeeze(vbar_ref(1,   :,       count_1:count_end)); 
  
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

toc
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

