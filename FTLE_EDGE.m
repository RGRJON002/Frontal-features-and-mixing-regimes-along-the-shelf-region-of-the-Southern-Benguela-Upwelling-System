% Figure of mean edge probabilities for summer and winter as well the mean
% FTLE field for 2014-2017

addpath '/home/jono/Documents/MATLAB/LCS-Tool-master'
addpath '/home/jono/Documents/PhD/CHAPTER_2'
addpath /home/jono/Documents/MATLAB/CMOCEAN
addpath /home/jono/CROCO/croco_tools/UTILITIES/m_map1.4h
addpath /media/data/DATASETS_CROCOTOOLS/m_map1.4f
addpath '/usr/local/MATLAB/R2020a/toolbox/matlab/imagesci/'
addpath '/media/data/DIAGNOSTICS/FRONTS'
addpath '/home/jono/CROCO/croco_tools/Preprocessing_tools'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load FTLE_zoom_nonorm.mat

[Y,MO,D] = datevec(mydate(1,:)');  
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

FTLE_sum = myFTLE(:,:,ind_sum);
sum_date = mydate(1,ind_sum);
FTLE_win = myFTLE(:,:,ind_win);
win_date = mydate(1,ind_win);

clear myFTLE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MY EDGE DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CROCO_path = '/media/jono/SBUS/SBUS_3km/CHPC_OUTPUT';
Ymin = 2004;
Ymax = 2018;
Yorig = 1990;

TIME = [];

for i = Ymin:Ymax
    for j = 1:12
    file = strcat(CROCO_path,'/','croco_avg_Y',string(i),'M',string(j),'.nc');
    if exist(file, 'file')
        disp(file)
        disp('Reading data')
        % Read in the time array
        time = ncread(file,'time');
        % Store arrays
        disp('time')
        TIME = cat(1,TIME,time);
    else
        disp(strcat('No data for',file))
    end
    end
end
    
clear time

[~, ia, ~] = unique(TIME);
TIME = TIME(ia);

CROCO_time = datetime(Yorig,1,1) + seconds(TIME);

load EDGE_MAT.mat
EDGE = pickle_data;
EDGE = permute(EDGE,[2 1 3]);
clear pickle_data

summer = [1, 2, 12];  winter = [6, 7, 8]; 

[Y,MO,D] = datevec(CROCO_time); 

ind_sum = [];
for i = 1:length(summer)
    ind_sum = [ind_sum; find(summer(i) == MO)];
end

ind_win = [];
for i = 1:length(winter)
    ind_win = [ind_win; find(winter(i) == MO)];
end

% Need to Patch to make up the size of the domain

EDGE(end+1,:,:) = EDGE(end,:,:);  % Basic patch to get array same size
EDGE(:,end+1,:) = EDGE(:,end,:);

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

% Now index all data to smaller array

lat_min = -33.5;
lat_max = -28;
lon_min = 15;
lon_max = 20;

[~,idx_lat_min]=min(abs(CROCO_lat(1,:)-lat_min));
[~,idx_lat_max]=min(abs(CROCO_lat(1,:)-lat_max));
[~,idx_lon_min]=min(abs(CROCO_lon(:,1)-lon_min));
[~,idx_lon_max]=min(abs(CROCO_lon(:,1)-lon_max));

CROCO_lat = CROCO_lat(idx_lon_min:idx_lon_max,idx_lat_min:idx_lat_max);
CROCO_lon = CROCO_lon(idx_lon_min:idx_lon_max,idx_lat_min:idx_lat_max);
CROCO_top = CROCO_top(idx_lon_min:idx_lon_max,idx_lat_min:idx_lat_max);
mask      = mask(idx_lon_min:idx_lon_max,idx_lat_min:idx_lat_max);

% Velocity data and EDGE
EDGE = EDGE(idx_lon_min:idx_lon_max,idx_lat_min:idx_lat_max,:);

% Create domain equal in resolution to size(FTLE) noting the zoom

domain = [15,20;-33.5,-28];
resolutionX = size(CROCO_lat,1)*5;
[resolutionY,deltaX] = equal_resolution(domain,resolutionX);
resolution = [resolutionX,resolutionY];

% Increase resolution

lat_tmp = min(CROCO_lat(1,:)):(max(CROCO_lat(1,:))-min(CROCO_lat(1,:)))/resolutionY:max(CROCO_lat(1,:));
lon_tmp = min(CROCO_lon(:,1)):(max(CROCO_lon(:,1))-min(CROCO_lon(:,1)))/resolutionX:max(CROCO_lon(:,1));

% Clean

lat_tmp(end) = [];
lon_tmp(end) = [];

[XX,YY] = meshgrid(lon_tmp',lat_tmp');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolate FTLE to CROCO grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:size(FTLE_sum,3)
    FTLE_SUM(:,:,i) = interp2(XX,YY,FTLE_sum(:,:,i),CROCO_lon',CROCO_lat');
end

FTLE_SUM(FTLE_SUM==0) = nan;

for i = 1:size(FTLE_win,3)
    FTLE_WIN(:,:,i) = interp2(XX,YY,FTLE_win(:,:,i),CROCO_lon',CROCO_lat');
end

FTLE_WIN(FTLE_WIN==0) = nan;

% EDGE average probability field

PROB_sum = nansum(EDGE(:,:,ind_sum),3)./size(EDGE(:,:,ind_sum),3);
PROB_win = nansum(EDGE(:,:,ind_win),3)./size(EDGE(:,:,ind_win),3);

%clear EDGE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

domain = [15,19;-33.5,-28];

figure
% Probability
cmin = 0;
cmax = 0.3;
levels = [cmin:0.05:cmax];
mydepths = [200,300,500];

m_proj('miller','long',[domain(1,:)]...
    ,'lat',[domain(2,:)]);
subplot(2,2,1);
m_pcolor(CROCO_lon,CROCO_lat,PROB_sum);
shading interp
m_gshhs_h('patch',[.7 .7 .7],'edgecolor','none');
m_grid('box','fancy','linest','none','tickdir','out','fontsize',13);
cmocean('dense',length(levels))
caxis([cmin cmax])
cRange=caxis;
hold on
m_contour(CROCO_lon,CROCO_lat,CROCO_top,mydepths,'Color','g','LineWidth',1,'LineStyle','--');
caxis(cRange)
title("(a) Summer [DJF]",'FontSize', 14)

subplot(2,2,2)
m_pcolor(CROCO_lon,CROCO_lat,PROB_win);
shading interp
m_gshhs_h('patch',[.7 .7 .7],'edgecolor','none');
m_grid('box','fancy','linest','none','tickdir','out','fontsize',13);
cmocean('dense',length(levels))
caxis([cmin cmax])
cRange=caxis;
hold on
m_contour(CROCO_lon,CROCO_lat,CROCO_top,mydepths,'Color','g','LineWidth',1,'LineStyle','--');
caxis(cRange)
title("(b) Winter [JJA]",'FontSize', 14)

% Create colorbar
hold on
ca = colorbar('eastOutside');
ca.Position = ca.Position + 1e-10;
ca.Label.String = 'Probability';
ca.FontSize = 12;
caxis([cmin cmax]);

hold off

% FTLE
cmin = 0;
cmax = 0.2;
levels = [cmin:0.05:cmax];
mydepths = [200,300,500];
m_proj('miller','long',[domain(1,:)]...
    ,'lat',[domain(2,:)]);
subplot(2,2,3);
m_contourf(CROCO_lon,CROCO_lat,nanmean(FTLE_SUM,3)',levels,'ShowText','on');
m_gshhs_h('patch',[.7 .7 .7],'edgecolor','none');
m_grid('box','fancy','linest','none','tickdir','out','fontsize',13);
cmocean('curl',length(levels))
caxis([cmin cmax])
cRange=caxis;
hold on
m_contour(CROCO_lon,CROCO_lat,CROCO_top,mydepths,'Color','g','LineWidth',1,'LineStyle','--');
caxis(cRange)
title("(c) Summer [DJF]",'FontSize', 14)

subplot(2,2,4)
m_contourf(CROCO_lon,CROCO_lat,nanmean(FTLE_WIN,3)',levels,'ShowText','on');
m_gshhs_h('patch',[.7 .7 .7],'edgecolor','none');
m_grid('box','fancy','linest','none','tickdir','out','fontsize',13);
cmocean('curl',length(levels))
caxis([cmin cmax])
cRange=caxis;
hold on
m_contour(CROCO_lon,CROCO_lat,CROCO_top,mydepths,'Color','g','LineWidth',1,'LineStyle','--');
caxis(cRange)
title("(d) Winter [JJA]",'FontSize', 14)

% Create colorbar
hold on
ca = colorbar('eastOutside');
ca.Position = ca.Position + 1e-10;
ca.Label.String = 'FTLE [days^{-1}]';
ca.FontSize = 12;
caxis([cmin cmax]);

% Save the figure

set(gcf, 'InvertHardcopy', 'off')
print('-f1','FTLE_EDGE_nonorm','-dpng','-r600');

