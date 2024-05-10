%%%% This script will plot the entire model domain as well as show the
%%%% zoomed region of interest, position of the SAMBA mooring, the SHBML as
%%%% well as show Cape Town and Cape Columbine

%%%%%%%%%%%%%%%%%%
% PROCESS CROCO
%%%%%%%%%%%%%%%%%%

addpath /home/jono/Documents/MATLAB/m_map
addpath /home/jono/Documents/MATLAB/CMOCEAN
addpath '/home/jono/Documents/MATLAB/CDT-master/cdt'
addpath '/home/jono/Documents/MATLAB/CDT-master/cdt/cdt_data'
addpath '/home/jono/Documents/MATLAB/tight_subplot'
addpath '/media/data/OBS_DATA/SHBML DATA'
addpath('/usr/local/MATLAB/R2020a/toolbox/matlab/imagesci/')

%%%%%%%%%%%%%%
% INPUTS
%%%%%%%%%%%%%%

CROCO_path = '/media/jono/SBUS/SBUS_3km/CHPC_OUTPUT';
Ymin = 2004;
Ymax = 2018;
Yorig = 1990;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SAMBA mooring data is located at

lon_M4 = 17.8803;
lat_M4 = -34.2924;

% SHBML

load 'SHBML_new.mat'

% Stations 8 and 11 are the Columbine stations
% From offshore to inshore: 1..2..3..4..5..6..7..9..10..12..13..14 (land)

for k = 1:14
    OBS = struct2cell(load(strcat('STATION_',string(k),'.mat')));
    OBS = OBS{1,1};
    % Get the individual variables from the file
    lon = OBS.lon;
    lat = OBS.lat;
    % Process the data and index date first    
    Latitude(:,k) = lat;
    Longitude(:,k) = lon;  
end

% SBUS mooring (70 m)

lat_SH_mooring = -32.329;
lon_SH_mooring = 18.318;

% Elands Bay Mooring (20 m)

lat_EL_mooring = -32.292;
lon_EL_mooring = 18.318;

%
mask=ncread(strcat(CROCO_path,'/','croco_avg_Y',string(Ymin),'M',string(1),'.nc'),'mask_rho');
mask(mask==0)=nan;
CROCO_lat=ncread(strcat(CROCO_path,'/','croco_avg_Y',string(Ymin),'M',string(1),'.nc'),'lat_rho');
CROCO_lon=ncread(strcat(CROCO_path,'/','croco_avg_Y',string(Ymin),'M',string(1),'.nc'),'lon_rho');
croco_top = ncread(strcat(CROCO_path,'/','croco_avg_Y',string(Ymin),'M',string(1),'.nc'),'h');
croco_top = croco_top.*mask;

lat_min = -36;
lat_max = -28;
lon_min = 15;
lon_max = 20;

% Figure

figure
m_proj('lambert','long',[lon_min lon_max],'lat',[lat_min lat_max]);
[CS,CH]= m_contourf(CROCO_lon,CROCO_lat,croco_top,[0:200:1000 1000:1000:5000],'edgecolor','none');
m_gshhs_f('patch',[.7 .7 .7],'edgecolor','none');
caxis([0 5000])
cRange=caxis;
hold on
m_contour(CROCO_lon,CROCO_lat,croco_top,[200,300,500],'Color','g','LineWidth',1,'LineStyle','--');
caxis(cRange)
hold on
h1=m_line(Longitude,Latitude,'marker','s','color',[0 .5 0],...
          'linest','none','markerfacecolor','w','clip','point');
hold on
h2 = m_line(lon_M4,lat_M4,'marker','o','color','r','linewi',2,...
          'linest','none','markersize',8,'markerfacecolor','w');
hold on
h3 = m_line(lon_SH_mooring, lat_SH_mooring,'marker','v','color','k','linewi',2,...
          'linest','none','markersize',8,'markerfacecolor','w');
hold on                   
h4 = m_line(lon_EL_mooring, lat_EL_mooring,'marker','o','color','b','linewi',2,...
          'linest','none','markersize',6,'markerfacecolor','w');
hold on                   
m_grid('linest','none','tickdir','out','box','fancy',...
    'fontsize',16,'BackgroundColor', [0.7 0.7 0.7]);
legend([h1(1),h2(1),h3(1),h4(1)],'SHBML','M4 SAMBA mooring','70 m Mooring','20 m Mooring','location','southwest','Fontsize',14);


colormap(flipud(m_colmap('blues')));  
caxis([000 5000]);

[ax,h]=m_contfbar([.55 .75],.8,CS,CH,'endpiece','no','axfrac',.05);
title(ax,'meters')

hold off
% Small subplot showing whole model domain

axes('position',[[0.182028265851786 0.441836155519825 0.200000000000001 0.200000000000002]])
m_proj('miller','long',[double(min(min(CROCO_lon))) double(max(max(CROCO_lon)))]...
    ,'lat',[double(min(min(CROCO_lat))) double(max(max(CROCO_lat)))]);
m_pcolor(CROCO_lon, CROCO_lat, croco_top);
caxis([000 5000]);
shading interp;
m_grid('linest','none','tickdir','out','box','fancy',...
    'fontsize',11,'BackgroundColor', [0.7 0.7 0.7]);
ax=gca;
set(gca, 'Color', [0.5 0.5 0.5])


set(gcf, 'InvertHardcopy', 'off')
print('-f1','Context','-dpng','-r600');
