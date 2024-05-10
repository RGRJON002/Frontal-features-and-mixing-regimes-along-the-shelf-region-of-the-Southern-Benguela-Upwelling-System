%%%% This script is inteded only to plot the composite 2x2 plot showing
%%%% avergae summer and winter circulation dynamics integrated over the
%%%% surface to 300 m and the normilzed forward FTLE computed over 2014 
%%%% to 2017 with an integration time of 30 days with a 10 day window

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER LIBRARIES and PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /home/jono/CROCO/croco_tools/UTILITIES/m_map1.4h
addpath /media/data/DATASETS_CROCOTOOLS/m_map1.4f
addpath /home/jono/Documents/MATLAB/CMOCEAN
addpath '/home/jono/Documents/MATLAB/CDT-master/cdt'
addpath '/home/jono/Documents/MATLAB/CDT-master/cdt/cdt_data'
addpath '/home/jono/Documents/MATLAB/tight_subplot'
addpath '/home/jono/Documents/MATLAB/LCS-Tool-master'
addpath '/usr/local/MATLAB/R2020a/toolbox/matlab/imagesci/'

% Will process FTLE and the STREAM seperately due to memory issues

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FTLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Laod my data

load('FTLE.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Average for summer and winter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Y,MO,D] = datevec(mydate');  
% Average all variable data for our seasons

summer = [1, 2, 12];  winter = [6, 7, 8]; 

ind_sum = [];
for i = 1:length(summer)
    ind_sum = [ind_sum; find(summer(i) == MO)];
end

ind_win = [];
for i = 1:length(winter)
    ind_win = [ind_win; find(winter(i) == MO)];
end

%%%% Compute the averages

FTLE_sum = nanmean(myFTLE(:,:,ind_sum),3);
FTLE_win = nanmean(myFTLE(:,:,ind_win),3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add topography
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CROCO_path = '/media/jono/SBUS/SBUS_3km/CHPC_OUTPUT';
Ymin = 2005;

%%%%%%%%%
% CROCO
%%%%%%%%%

% Once-off gridding data and dimensions
file = strcat(CROCO_path,'/','croco_avg_Y',string(Ymin),'M',string(1),'.nc');
mask=ncread(file,'mask_rho');
mask(mask==0)=nan;
CROCO_top = ncread(file,'h').*mask;
CROCO_lat=ncread(strcat(CROCO_path,'/','croco_avg_Y',string(Ymin),'M',string(1),'.nc'),'lat_rho');
CROCO_lon=ncread(strcat(CROCO_path,'/','croco_avg_Y',string(Ymin),'M',string(1),'.nc'),'lon_rho');

% Create domain equal in resolution to size(FTLE) noting the zoom

domain = [15,20;-36,-28];
resolutionX = 300;
[resolutionY,deltaX] = equal_resolution(domain,resolutionX);
resolution = [resolutionX,resolutionY];

lat_min = -36;
lat_max = -28;
lon_min = 15;
lon_max = 20;

[~,idx_lat_min]=min(abs(CROCO_lat(1,:)-lat_min));
[~,idx_lat_max]=min(abs(CROCO_lat(1,:)-lat_max));
[~,idx_lon_min]=min(abs(CROCO_lon(:,1)-lon_min));
[~,idx_lon_max]=min(abs(CROCO_lon(:,1)-lon_max));

% Index the variables Spatially

CROCO_lat = CROCO_lat(idx_lon_min:idx_lon_max,idx_lat_min:idx_lat_max);
CROCO_lon = CROCO_lon(idx_lon_min:idx_lon_max,idx_lat_min:idx_lat_max);
CROCO_top = CROCO_top(idx_lon_min:idx_lon_max,idx_lat_min:idx_lat_max);
mask      = mask(idx_lon_min:idx_lon_max,idx_lat_min:idx_lat_max);

% Increase resolution

lat_tmp = min(CROCO_lat(1,:)):(max(CROCO_lat(1,:))-min(CROCO_lat(1,:)))/resolutionY:max(CROCO_lat(1,:));
lon_tmp = min(CROCO_lon(:,1)):(max(CROCO_lon(:,1))-min(CROCO_lon(:,1)))/resolutionX:max(CROCO_lon(:,1));

% Clean

lat_tmp(end) = [];
lon_tmp(end) = [];

[XX,YY] = meshgrid(lon_tmp',lat_tmp');

% Interpolate FTLE to CROCO grid

FTLE_sum = interp2(XX,YY,FTLE_sum,CROCO_lon',CROCO_lat');
FTLE_sum(FTLE_sum==0) = nan;

FTLE_win = interp2(XX,YY,FTLE_win,CROCO_lon',CROCO_lat');
FTLE_win(FTLE_win==0) = nan;

clear myFTLE    % Save space

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STREAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Summer

load('summer_transport_croco.mat'); 

% Average all variable data for our seasons

% Index U, V and psi

U_sum = double(nanmean(U_sum,3));
V_sum = double(nanmean(V_sum,3));
PSI_sum = double(nanmean(PSI_sum,3));

% Winter

load('winter_transport_croco.mat'); 

% Average all variable data for our seasons

% Index U, V and psi

U_win = double(nanmean(U_win,3));
V_win = double(nanmean(V_win,3));
PSI_win = double(nanmean(PSI_win,3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

start   % Have to tun to get libraries for STREAM
%crocotools_param
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FTLE

cmin = 0;
cmax = 0.8;
levels = [cmin:0.1:cmax];
mydepths = [200,300,500];

figure
m_proj('miller','long',[domain(1,:)]...
    ,'lat',[domain(2,:)]);
subplot(2,2,3);
m_contourf(CROCO_lon,CROCO_lat,FTLE_sum',levels,'ShowText','on');
m_gshhs_l('patch',[.7 .7 .7],'edgecolor','none');
m_grid('box','fancy','linest','none','tickdir','out','fontsize',13);
cmocean('curl',length(levels))
caxis([cmin cmax])
cRange=caxis;
hold on
m_contour(CROCO_lon,CROCO_lat,CROCO_top,mydepths,'Color','g','LineWidth',1,'LineStyle','--');
caxis(cRange)
title("(c) FTLE [DJF]",'FontSize', 14)

subplot(2,2,4)
m_contourf(CROCO_lon,CROCO_lat,FTLE_win',levels,'ShowText','on');
m_gshhs_l('patch',[.7 .7 .7],'edgecolor','none');
m_grid('box','fancy','linest','none','tickdir','out','fontsize',13);
cmocean('curl',length(levels))
caxis([cmin cmax])
cRange=caxis;
hold on
m_contour(CROCO_lon,CROCO_lat,CROCO_top,mydepths,'Color','g','LineWidth',1,'LineStyle','--');
caxis(cRange)
title("(d) FTLE [JJA]",'FontSize', 14)

% Create colorbar
hold on
ca = colorbar('eastOutside');
ca.Position = ca.Position + 1e-10;
ca.Label.String = 'FTLE [normalised]';
ca.FontSize = 12;
caxis([cmin cmax]);

% STREAM
npts=[2 2 2 2];

lonmin=15;
lonmax=20;
latmin=-36;
latmax=-28;

topo=1;

subplot(2,2,1)  % Summer

m_proj('mercator',...
       'lon',[lonmin lonmax],...
       'lat',[latmin latmax]);
[x,y]=m_ll2xy(lon,lat,'clip','off');
hold on

psi_r=1e-6*psi2rho(PSI_sum).*mask;
u = U_sum;
v = V_sum;

colmax = 22;
colmin = 6;
dcol = 1;

[C1,h1]=m_contour(rempoints(lon,npts),rempoints(lat,npts),rempoints(psi_r,npts)...
           ,[dcol:dcol:colmax],'k');
if ~isempty(C1)
 % clabel(C1,h1,'LabelSpacing',1000,'Rotation',0)
  hf1=add_streamarrows(C1,x,y,u2rho_2d(u),v2rho_2d(v));
  hold on
end

[C2,h2]=m_contour(rempoints(lon,npts),rempoints(lat,npts),rempoints(-psi_r,npts)...
           ,[dcol:dcol:colmin],'k');
if ~isempty(C2)
 % clabel(C2,h2,'LabelSpacing',1000,'Rotation',0)
  hf2=add_streamarrows(C2,x,y,u2rho_2d(u),v2rho_2d(v));
  hold on
end

[C3,h3]=m_contour(rempoints(lon,npts),rempoints(lat,npts),rempoints(psi_r,npts)...
           ,[0 0],'k');
if ~isempty(C3)
 % clabel(C3,h3,'LabelSpacing',1000,'Rotation',0)
  hf3=add_streamarrows(C3,x,y,u2rho_2d(u),v2rho_2d(v));
  set(h3,'Color','k','Linewidth',1.2);
  set(hf3,'Linewidth',1.2);
  hold on
end

[C4,h4]=m_contour(rempoints(lon,npts),rempoints(lat,npts),rempoints(psi_r,npts)...
           ,[20:0.1:21],'k');
if ~isempty(C4)
 % clabel(C4,h4,'LabelSpacing',1000,'Rotation',0)
  hf4=add_streamarrows(C4,x,y,u2rho_2d(u),v2rho_2d(v));
  set(h4,'Color','k','Linewidth',1.2);
  set(hf4,'Linewidth',1.2);
  hold on
end

if topo==1
  [C5,h5]=m_contour(rempoints(lon,npts),rempoints(lat,npts),rempoints(h,npts)...
           ,[200 300 500],'--g');
  %clabel(C5,h5,'LabelSpacing',1000,'Rotation',0)
end

m_usercoast('coastline_l.mat','patch',[.9 .9 .9]);
hold off
m_gshhs_h('patch',[.8 .8 .8],'edgecolor','none');
m_grid('box','fancy','xtick',5,'ytick',5,'tickdir','out');
set(findobj('tag','m_grid_color'),'facecolor','white')

title(['(a) ','Streamlines [DJF]'],'FontSize',14)
hold off

subplot(2,2,2)

m_proj('mercator',...
       'lon',[lonmin lonmax],...
       'lat',[latmin latmax]);
[x,y]=m_ll2xy(lon,lat,'clip','off');
hold on

psi_r=1e-6*psi2rho(PSI_win).*mask;
u = U_win;
v = V_win;

colmax = 22;  % 11 winter
colmin = 6;
dcol=1;

[C1,h1]=m_contour(rempoints(lon,npts),rempoints(lat,npts),rempoints(psi_r,npts)...
           ,[dcol:dcol:colmax],'k');
if ~isempty(C1)
 % clabel(C1,h1,'LabelSpacing',1000,'Rotation',0)
  hf1=add_streamarrows(C1,x,y,u2rho_2d(u),v2rho_2d(v));
  hold on
end

[C2,h2]=m_contour(rempoints(lon,npts),rempoints(lat,npts),rempoints(-psi_r,npts)...
           ,[dcol:dcol:colmin],'k');
if ~isempty(C2)
 % clabel(C2,h2,'LabelSpacing',1000,'Rotation',0)
  hf2=add_streamarrows(C2,x,y,u2rho_2d(u),v2rho_2d(v));
  hold on
end

[C3,h3]=m_contour(rempoints(lon,npts),rempoints(lat,npts),rempoints(psi_r,npts)...
           ,[0 0],'k');
if ~isempty(C3)
 % clabel(C3,h3,'LabelSpacing',1000,'Rotation',0)
  hf3=add_streamarrows(C3,x,y,u2rho_2d(u),v2rho_2d(v));
  set(h3,'Color','k','Linewidth',1.2);
  set(hf3,'Linewidth',1.2);
  hold on
end

[C4,h4]=m_contour(rempoints(lon,npts),rempoints(lat,npts),rempoints(psi_r,npts)...
           ,[20:0.1:21],'k');
if ~isempty(C4)
 % clabel(C4,h4,'LabelSpacing',1000,'Rotation',0)
  hf4=add_streamarrows(C4,x,y,u2rho_2d(u),v2rho_2d(v));
  set(h4,'Color','k','Linewidth',1.2);
  set(hf4,'Linewidth',1.2);
  hold on
end

if topo==1
  [C5,h5]=m_contour(rempoints(lon,npts),rempoints(lat,npts),rempoints(h,npts)...
           ,[200 300 500],'--g');
  %clabel(C5,h5,'LabelSpacing',1000,'Rotation',0)
end

m_usercoast('coastline_l.mat','patch',[.9 .9 .9]);
hold off
m_gshhs_h('patch',[.8 .8 .8],'edgecolor','none');
m_grid('box','fancy','xtick',5,'ytick',5,'tickdir','out');
set(findobj('tag','m_grid_color'),'facecolor','white')

title(['(b) ','Streamlines [JJA]'],'FontSize',14)
hold off

% Save the figure

set(gcf, 'InvertHardcopy', 'off')
print('-f1','FTLE_STREAM','-dpng','-r600');




