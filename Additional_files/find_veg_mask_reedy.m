clear all ; close all ; clc ;
format long 
% written by Tarandeep Kalra July 2, 2020
% Zafer shared veg contour boundaries 
% using polyout function 

% Loading Zafer's file
load('/media/taran/DATADRIVE2/marsh_result/barnegat_bay/veg_fill/dc_rc_veg.mat');
r1=rc_veg(:,1); r2=rc_veg(:,2);
% 
 % POLYOUT function to get areas with veg
xx=rc_veg(:,1);
yy=rc_veg(:,2);
pgon=polyshape(xx,yy);

%polyout = holes(pgon) % get the holes and then you can get their verticees

figure(1)
plot(pgon)
veg_vert=pgon(1:end).Vertices;
%save('veg_vert_poly_dinner.mat','polyout')
 
 url=('bbleh_reedy_grd.nc'); 
 lon_rho=ncread(url,'lon_rho')  ;
 lat_rho=ncread(url,'lat_rho')  ; 
 h      =ncread(url,'h')       ;
 veg_mask=lon_rho*0.0; 
 
% 
 [mx, my] = size(lon_rho);
 grid_size_coarse=mx*my        ;
% 
 veg_mask_1d=reshape(veg_mask,[grid_size_coarse 1]); 
 lon_rho_1d=reshape(lon_rho,[grid_size_coarse 1]);
 lat_rho_1d=reshape(lat_rho,[grid_size_coarse 1]); 
% 
 xv=lon_rho_1d;
 yv=lat_rho_1d; 
% 
% 
 for i=1:length(pgon)
%     
     xx=pgon(i).Vertices(1:end,1);
     yy=pgon(i).Vertices(1:end,2);
%     
     in = inpolygon(xv,yv,xx,yy);
     % have a temporrary maask to idenitfy regions with veg 
     veg_mask_1d(in)=1 ;
 end 
 veg_mask_2d=reshape(veg_mask_1d,mx,my);
% %  
% % %  
% figure(2)
% plot(xv,yv,'g+')
% hold on 
% plot(xx,yy,'r+');
% hold on
% plot(xv(in), yv(in),'.r')
% 
figure(3)
 pcolorjw(lon_rho, lat_rho, veg_mask_2d)
  print('-dpng','-r200','veg_mask_reedy_Zafer_bdry.png')

% 
% figure(3)
% pcolorjw(lon_rho, lat_rho, h)
% 
% figure(4)
% plot(r1, r2)
% %plot(xv(~in),yv(~in),'.b')
% 
% %figure(1)
% %plot(lon_rho_1d, lat_rho_1d, 'bp')                                          % Plot All Points
% %hold on
% %plot(xx, yy, '-r')                                          % Plot Polygon
% %plot(xcoord, ycoord, 'gp')                                  % Overplot ‘inon’ Points
% %hold off
% %plot(polyout)
% %hold on
% %plot(xx,yy,'r+')
