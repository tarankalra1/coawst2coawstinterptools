ncclear
clear all ; close all ; clc; 
setup_nctoolbox
echo off
 % Original code by Zafer Defne
 % to compute the TPAR files from a given COAWST solution that is avaliable
 % on thredds to drive a smaller/refined COAWST grid
 % Modified by Tarandeep S Kalra and added comments 
 
% To get TPAR Files the big grid should have save Hwave, PWave_top,Dwave 
% Geoport link to big grid 
url= 'http://geoport.whoi.edu/thredds/dodsC/sand/usgs/users/tkalra/bbleh/bbleh.ncml';
nc=ncgeodataset(url);
m=nc{'mask_rho'}(:);
%lat=nc{'lat_rho'}(:);
%
% This is lon and lat of the outer grid/biger grid saved 
%  
g = nc{'zeta'}(1,:,:).grid;
lon_rho = g.lon;
lat_rho = g.lat;

lon_rho(m==0)=nan;
lat_rho(m==0)=nan;
dstart=datenum(2015,04,15);

%time=nc{'ocean_time'}(:); 
t_insec(:)=nc{'ocean_time'}(:);
t_indays(:)=t_insec(:)/(3600*24);
% Converting time in days to julian time units;
t(:)=t_indays(:) +dstart;

% The time period for which the TPAR files need to be generated
t1=datenum(2015,07,01);
t2=datenum(2015,07,03);
ti1=near(t,t1);
ti2=near(t,t2);

time=t(ti1:ti2);

% These are grid parameters for the refined grid/inner grid
% gridname : refined grid
grid_name_inner = 'C:\Users\tkalra\Desktop\barnegatbay\easygrid\bbleh_reedy_grd.nc';

lat1=ncread(grid_name_inner,'lat_rho');
lon1=ncread(grid_name_inner,'lon_rho'); 
m1=ncread(grid_name_inner,'mask_rho'); 

% % 
% use the wetting and drying mask to nan out points
%lat1(m1==0)=nan;
%lon1(m1==0)=nan;

 [ii, jj]=size(lon1);
% 
% User input: Grid specific  
% 
% % go through the eastern boundary i.e. j=1 and all i points along the
% eastern boundary edge beacuse we want wave forcing from eastern edge. 

jj=1; 
 for i=1:ii
     [ixcn(i), iycn(i)]=find_nearest_point(lon1(i,jj),lat1(i,jj), lon_rho,....
                                            lat_rho, m, 'auto');
 end
%discretize along the eastern edge and that would be the total number of tpar files
% the small grid 134/2 = 67 

% for 134 points we chose to create 67 segments=67 TPAR files ..for each
% TPAR file need to find one nearest point from which data can be obtained
% from big grid to drive wave forcing in the small grid 
% therefore use the "mode" function that gets the most frequent nearest
% point data 

% also 67 segments would imply that there are 2 points each segment =134
% points along the eastern edge so information in TPAR file would be stored
% as 1:2 and 3:4 and 5:6... and so on 133:134


 mult=2; 
for i=1:67
     etac(i)=mode(iycn(i*mult-1:i*mult));
     xic(i)=mode(ixcn (i*mult-1:i*mult));
end

 for i=1:length(etac)
   TPAR(1:length(time),2)=squeeze(nc{'Hwave'}(ti1:ti2,etac(i),xic(i)));

%   TPAR(1:length(time),3)=squeeze(nc{'Pwave_top'}(ti1:ti2,etac(i),xic(i)));
   
   TPAR(1:length(time),4)=squeeze(nc{'Dwave'}(ti1:ti2,etac(i),xic(i)));

    TPAR(1:length(time),5)=20;
    TPAR_time(1:length(time),1)=str2num(datestr(time,'yyyymmdd.HHMM'));
    l=sprintf('%d',i);
    ofile=['TPAR',l,'.txt'];
    fid=fopen(ofile,'w');
    fprintf(fid,'TPAR \n');
    for wavet=1:length(time)
        fprintf(fid,'%8.4f',TPAR_time(wavet));
        fprintf(fid,'         %3.2f        %3.2f     %3.1f       %2.f\n',TPAR(wavet,2:5)');
    end
    fclose(fid);
end
