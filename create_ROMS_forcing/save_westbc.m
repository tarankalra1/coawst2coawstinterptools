function [zeta_west,ubar_west,vbar_west,temp_west,salt_west,......
          u_west, v_west, sand_west_01, sand_west_02, .....
          sand_west_03, sand_west_04, sand_west_05]=save_westbc(obc,...
            zeta_ref, ubar_ref, vbar_ref, u_4d, v_4d, .....
            temp_4d, salt_4d, .....
            sand01_4d, sand02_4d, sand03_4d, sand04_4d, sand05_4d, count_1, count_end)

if(obc.west==1)
   zeta_west = squeeze(zeta_ref(:,       count_1:count_end));
   ubar_west = squeeze(ubar_ref(:,       count_1:count_end));
   vbar_west = squeeze(vbar_ref(:,       count_1:count_end));
 
   temp_west  = squeeze(temp_4d(:,  :,count_1:count_end));
   salt_west  = squeeze(salt_4d(:,  :,count_1:count_end));
   u_west=squeeze (u_4d(:, :,        count_1:count_end));
   v_west=squeeze (v_4d(:, :,        count_1:count_end));
 
   sand_west_01=squeeze(sand01_4d(:,:,count_1:count_end));
   sand_west_02=squeeze(sand02_4d(:,:,count_1:count_end));
   sand_west_03=squeeze(sand03_4d(:,:,count_1:count_end));
   sand_west_04=squeeze(sand04_4d(:,:,count_1:count_end));
   sand_west_05=squeeze(sand05_4d(:,:,count_1:count_end));
 
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

