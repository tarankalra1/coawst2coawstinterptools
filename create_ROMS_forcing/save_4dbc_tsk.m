function [zeta_north, ubar_north, vbar_north, temp_north, ....
          salt_north, u_north, v_north, ......
	  sand_north_01, sand_north_02, sand_north_03, .....
	  sand_north_04, sand_north_05,........
	  zeta_south, ubar_south, vbar_south, .....
	  temp_south, salt_south, u_south, v_south, ......
	  sand_south_01, sand_south_02, sand_south_03, .....
	  sand_south_04, sand_south_05,........
      zeta_east, ubar_east, vbar_east, temp_east, ....
      salt_east, u_east, v_east, ......
	  sand_east_01, sand_east_02, sand_east_03, .....
	  sand_east_04, sand_east_05,........
      zeta_west, ubar_west, vbar_west, temp_west, ....
      salt_west, u_west, v_west, ......
	  sand_west_01, sand_west_02, sand_west_03, .....
	  sand_west_04, sand_west_05]=save_4dbc_tsk(zeta_ref, .....
                            ubar_ref, vbar_ref, .....
                            temp_4d, salt_4d, u_4d, v_4d,......
                            sand01_4d, sand02_4d, sand03_4d,....
                            sand04_4d, sand05_4d, count_1, count_end, obc);
                        
if(obc.north==1)	
   disp('in north')
%   zeta_north = squeeze(zeta_ref(:, end,      count_1:count_end));
%   ubar_north = squeeze(ubar_ref(:, end,      count_1:count_end));
%   vbar_north = squeeze(vbar_ref(:, end,      count_1:count_end));
% 
%   temp_north = squeeze(temp_4d(:, end, :,count_1:count_end));
%   salt_north = squeeze(salt_4d(:, end, :,count_1:count_end));
%   u_north=squeeze(u_4d(:, end,:,         count_1:count_end));
%   v_north=squeeze(v_4d(:, end,:,         count_1:count_end));
% 
%   sand_north_01=squeeze(sand01_4d(:,end,:,count_1:count_end));
%   sand_north_02=squeeze(sand02_4d(:,end,:,count_1:count_end));
%   sand_north_03=squeeze(sand03_4d(:,end,:,count_1:count_end));
%   sand_north_04=squeeze(sand04_4d(:,end,:,count_1:count_end));
%   sand_north_05=squeeze(sand05_4d(:,end,:,count_1:count_end));

  zeta_north = squeeze(zeta_ref(:,       count_1:count_end));
  ubar_north = squeeze(ubar_ref(:,       count_1:count_end));
  vbar_north = squeeze(vbar_ref(:,       count_1:count_end));

  temp_north = squeeze(temp_4d(:, :,count_1:count_end));
  salt_north = squeeze(salt_4d(:, :,count_1:count_end));
  u_north=squeeze(u_4d(:, :,         count_1:count_end));
  v_north=squeeze(v_4d(:, :,         count_1:count_end));

  sand_north_01=squeeze(sand01_4d(:,:,count_1:count_end));
  sand_north_02=squeeze(sand02_4d(:,:,count_1:count_end));
  sand_north_03=squeeze(sand03_4d(:,:,count_1:count_end));
  sand_north_04=squeeze(sand04_4d(:,:,count_1:count_end));
  sand_north_05=squeeze(sand05_4d(:,:,count_1:count_end));
else(obc.north==0)
  zeta_north=0.0; 
  ubar_north=0.0; 
  vbar_north=0.0; 
  temp_north=0.0 ; 
  salt_north=0.0 ;
  u_north=0.0; 
  v_north=0.0; 
  sand_north_01=0.0; 
  sand_north_02=0.0; 
  sand_north_03=0.0; 
  sand_north_04=0.0; 
  sand_north_05=0.0; 
end

if(obc.south==1)
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

else(obc.south==0)
    disp('here')
  zeta_south=0.0; 
  ubar_south=0.0; 
  vbar_south=0.0; 
  temp_south=0.0 ; 
  salt_south=0.0 ;
  u_south=0.0; 
  v_south=0.0; 
  sand_south_01=0.0; 
  sand_south_02=0.0; 
  sand_south_03=0.0; 
  sand_south_04=0.0; 
  sand_south_05=0.0; 
end

if(obc.east==1)
  zeta_east = squeeze(zeta_ref(:,       count_1:count_end));
  ubar_east = squeeze(ubar_ref(:,       count_1:count_end));
  vbar_east = squeeze(vbar_ref(:,       count_1:count_end));

  temp_east  = squeeze(temp_4d(:,  :,count_1:count_end));
  salt_east  = squeeze(salt_4d(:,  :,count_1:count_end));
  u_east=squeeze (u_4d(:, :,        count_1:count_end));
  v_east=squeeze (v_4d(:, :,        count_1:count_end));

  sand_east_01=squeeze(sand01_4d(:,:,count_1:count_end));
  sand_east_02=squeeze(sand02_4d(:,:,count_1:count_end));
  sand_east_03=squeeze(sand03_4d(:,:,count_1:count_end));
  sand_east_04=squeeze(sand04_4d(:,:,count_1:count_end));
  sand_east_05=squeeze(sand05_4d(:,:,count_1:count_end));

else(obc.east==0)
  zeta_east=0.0; 
  ubar_east=0.0; 
  vbar_east=0.0; 
  temp_east=0.0 ; 
  salt_east=0.0 ;
  u_east=0.0; 
  v_east=0.0; 
  sand_east_01=0.0; 
  sand_east_02=0.0; 
  sand_east_03=0.0; 
  sand_east_04=0.0; 
  sand_east_05=0.0; 
end

if(obc.west==1)
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

else(obc.west==0)
  zeta_west=0.0; 
  ubar_west=0.0; 
  vbar_west=0.0; 
  temp_west=0.0 ; 
  salt_west=0.0 ;
  u_west=0.0; 
  v_west=0.0; 
  sand_west_01=0.0; 
  sand_west_02=0.0; 
  sand_west_03=0.0; 
  sand_west_04=0.0; 
  sand_west_05=0.0; 
end
