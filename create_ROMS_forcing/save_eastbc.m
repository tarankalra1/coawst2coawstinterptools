function [zeta_east,ubar_east,vbar_east,temp_east,salt_east,......
          u_east, v_east, sand_east_01, sand_east_02, .....
          sand_east_03, sand_east_04, sand_east_05]=save_eastbc(obc,...
            zeta_ref, ubar_ref, vbar_ref, u_4d, v_4d, .....
            temp_4d, salt_4d, .....
            sand01_4d, sand02_4d, sand03_4d, sand04_4d, sand05_4d, count_1, count_end)

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

