function [sand02_tmp_coarse, sand03_tmp_coarse, .....
      sand04_tmp_coarse, sand05_tmp_coarse, .....
	  u2_tmp_coarse, v2_tmp_coarse, .....
      zr_coarse]=maplev_4dvar_tsk(sand02_coarse, sand03_coarse, ......
                                  sand04_coarse, sand05_coarse, .......
                                  u_coarse, v_coarse, zr_coarse, .....
                                  coarse_N,.....
                                  coarse_maskr, coarse_masku, coarse_maskv)
	     
for zz = 1:coarse_N
%       aa = squeeze(sand01_coarse(zz,:,:));
%       aa(~coarse_maskr) = NaN;
%       aa= maplev(aa);
%       sand01_tmp_coarse(:,:,zz)=aa;
%       clear aa
      
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
 
%      aa = squeeze(temp_coarse(zz,:,:));
%      aa(~coarse_maskr) = NaN;
%      aa = maplev(aa);
%      temp_tmp_coarse(:,:,zz) = aa;
%      clear aa
%       
%      aa = squeeze(salt_coarse(zz,:,:));
%      aa(~coarse_maskr) = NaN;
%      aa = maplev(aa);
%      salt_tmp_coarse(:,:,zz) = aa;
%      clear aa ;
%       

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
%
