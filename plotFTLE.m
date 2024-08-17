%%%% This script will plot the seasonal mean FTLE field from the daily 
%%%% .mat file. Also, a spatial average of the FTLE field over the shelf
%%%% to show a time-series plot will be useful to highlight the difference
%%%% in mixing regimes along the shelf. Furtherr, the edges will also be
%%%% overlain to show a probability map. This is similar to the FTLE_EDGE
%%%% script in the CHAPTER_1 directory, only this script adds some minor
%%%% imporvments to the plot

addpath /home/jono/CROCO/croco_tools/UTILITIES/m_map1.4h
addpath /media/data/DATASETS_CROCOTOOLS/m_map1.4f
addpath /home/jono/Documents/MATLAB/CMOCEAN
addpath /usr/local/MATLAB/R2020a/toolbox/matlab/imagesci/
addpath '/home/jono/Documents/MATLAB/CDT-master/cdt'
addpath '/home/jono/Documents/MATLAB/CDT-master/cdt/cdt_data'
addpath '/media/data/DIAGNOSTICS/FRONTS'  % EDGE FILE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS CROCO grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Need one CROCO file for the topography data in order to do the mask

CROCO_path = '/media/jono/SBUS/SBUS_3km/CHPC_OUTPUT';
Ymin = 2004;
Ymax = 2018;
Yorig = 1990;

% POCESS TIME DATA

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

%%%%%%%%%
% CROCO
%%%%%%%%%

% Once-off gridding data and dimensions
file = strcat(CROCO_path,'/','croco_avg_Y',string(Ymin),'M',string(1),'.nc');
% Read in the mask
mask=double(ncread(file,'mask_rho'));
mask(mask==0)=nan;   % land = 0, make nan
% Define my region, will focus on the coastal domain
CROCO_lat=ncread(strcat(CROCO_path,'/','croco_avg_Y',string(Ymin),'M',string(1),'.nc'),'lat_rho');
CROCO_lon=ncread(strcat(CROCO_path,'/','croco_avg_Y',string(Ymin),'M',string(1),'.nc'),'lon_rho');
CROCO_top = ncread(strcat(CROCO_path,'/','croco_avg_Y',string(Ymin),'M',string(1),'.nc'),'h');
%
CROCO_top=CROCO_top.*mask;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subset the domain of interest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The whole shelf region: needed in order to process EDGE first as FTLE is
% already on the smaller grid

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

%% Load in the EDGE DATA

load EDGE_ABS_MAT.mat
EDGE = pickle_data;
EDGE = permute(EDGE,[2 1 3]);
clear pickle_data

% Need to Patch to make up the size of the domain

EDGE(end+1,:,:) = EDGE(end,:,:);  % Basic patch to get array same size
EDGE(:,end+1,:) = EDGE(:,end,:);

% We know now the edge map will have the same dimensions. Just need to cut
% down the lon and lat variable

% SEASONAL MEAN

PROB_sum = sum(EDGE(:,:,ind_sum),3,'omitnan')./size(EDGE(:,:,ind_sum),3);
PROB_win = sum(EDGE(:,:,ind_win),3,'omitnan')./size(EDGE(:,:,ind_win),3);

% Cut-dwon to FTLE region

lat_min = -33.5;
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
PROB_sum = PROB_sum(idx_lon_min:idx_lon_max,idx_lat_min:idx_lat_max);
PROB_win = PROB_win(idx_lon_min:idx_lon_max,idx_lat_min:idx_lat_max);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FTLE

% Chose a file for now
load FTLE_zoom_nonorm_2014.mat  % Contains the time-series of FTLE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%
% DATA
%%%%%%%%%

mydate = mydate(1,:)';

summer = [1, 2, 12];  winter = [6, 7, 8]; 

[Y,MO,D] = datevec(mydate); 

ind_sum = [];
for i = 1:length(summer)
    ind_sum = [ind_sum; find(summer(i) == MO)];
end

ind_win = [];
for i = 1:length(winter)
    ind_win = [ind_win; find(winter(i) == MO)];
end

FTLE_SUM = mean(myFTLE(:,:,ind_sum),3,'omitnan');
FTLE_WIN = mean(myFTLE(:,:,ind_win),3,'omitnan');
%% 

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
cmax = 0.25;
levels = [cmin:0.05:cmax];
mydepths = [200,300,500];

m_proj('miller','long',[domain(1,:)]...
    ,'lat',[domain(2,:)]);
subplot(2,2,3);
m_contourf(CROCO_lon,CROCO_lat,FTLE_SUM',levels,'ShowText','on');
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
m_contourf(CROCO_lon,CROCO_lat,FTLE_WIN',levels,'ShowText','on');
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
%print('-f1','FTLE_EDGE_nonorm','-dpng','-r600');

%% Compute diagnostics using the 500 isobath to delineate the shelf region

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get mask for shelf region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Quick
croco_top = CROCO_top;

% Whole shelf
shelf_depth = 500;
h_mask = croco_top;
h_mask(h_mask > shelf_depth) = NaN;
mask_canny = ~isnan(h_mask);
mask_canny = double(mask_canny); 
mask_canny(mask_canny ==0) = nan;

% Nearshore
shal_depth = 200;
h_mask = croco_top;
h_mask(h_mask > shal_depth) = NaN;
shal_mask = ~isnan(h_mask);
shal_mask = double(shal_mask); 
shal_mask(shal_mask ==0) = nan;

% Shelf-break
h_mask = croco_top;
h_mask(h_mask > shelf_depth) = NaN;
h_mask(h_mask < shal_depth) = NaN;
shelf_mask = ~isnan(h_mask);
shelf_mask = double(shelf_mask); 
shelf_mask(shelf_mask ==0) = nan;

%%%%%%%%%%%%%%%%%%%%%%%% 
% Compute time-series
%%%%%%%%%%%%%%%%%%%%%%%%

% Start with the whole shelf
FTLE_TS_all = squeeze(mean(myFTLE.*mask_canny',[1,2],'omitnan'));
% Now do the upwelling region
FTLE_TS_up = squeeze(mean(myFTLE.*shal_mask',[1,2],'omitnan'));
% Now do the shelf-break
FTLE_TS_shf = squeeze(mean(myFTLE.*shelf_mask',[1,2],'omitnan'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save these arays
data = struct('mydate',mydate,'FTLE_TS_all',FTLE_TS_all,...
    'FTLE_TS_up',FTLE_TS_up,'FTLE_TS_shf',FTLE_TS_shf);  % Create data structure
save('FTLE_TS.mat','-struct','data'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Clim
FTLE_TS_all =  climatology(FTLE_TS_all,mydate(1,:),'monthly');
FTLE_TS_up =  climatology(FTLE_TS_up,mydate(1,:),'monthly');
FTLE_TS_shf =  climatology(FTLE_TS_shf,mydate(1,:),'monthly');

% Now do the plot
figure
xrange = [1:12];
plot(xrange,FTLE_TS_all,'--x','LineWidth',2,'Color','k')
hold on
plot(xrange,FTLE_TS_up,'-o','LineWidth',2,'Color','r')
hold on
plot(xrange,FTLE_TS_shf,'-s','LineWidth',2,'Color','b');
xticks([xrange])
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
xlabel('Months','FontSize',13)
ylabel('FTLE [day^-1]','FontSize',13)
legend('Whole shelf','Upwelling region','Shelf-break region')

% Save the figure

set(gcf, 'InvertHardcopy', 'off')
print('-f1','FTLE_TS','-dpng','-r600');

