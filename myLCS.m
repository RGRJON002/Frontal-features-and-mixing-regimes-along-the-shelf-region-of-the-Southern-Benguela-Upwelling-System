%%%% This script will plot LCS and overrlay frontal edges for summer and
%%%% winter. Main focus os for events in summer

addpath '/home/jono/Documents/MATLAB/LCS-Tool-master'
addpath '/home/jono/Documents/PhD/CHAPTER_2'
addpath /home/jono/Documents/MATLAB/CMOCEAN
addpath /home/jono/CROCO/croco_tools/UTILITIES/m_map1.4h
addpath /media/data/DATASETS_CROCOTOOLS/m_map1.4f
addpath '/usr/local/MATLAB/R2020a/toolbox/matlab/imagesci/'
addpath '/media/data/DIAGNOSTICS/FRONTS'
addpath '/home/jono/CROCO/croco_tools/Preprocessing_tools'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add velcoity field and edge positions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% VELOCITY DATA

CROCO_path = '/media/jono/SBUS/SBUS_3km/CHPC_OUTPUT';
Ymin = 2014;
Ymax = 2014;
Yorig = 1990;

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
CROCO_top=CROCO_top.*mask;
% Only need to index lon as lat bounds and eastern extent stay the same

u = [];
v = [];
TIME = [];

for i = Ymin:Ymax
    for j = 1:12
    file = strcat(CROCO_path,'/','croco_avg_Y',string(i),'M',string(j),'.nc');
    if exist(file, 'file')
        disp(file)
        disp('Reading data')
        % Read in the surface velocity
        ut = ncread(file,'u',[1,1,60,1],[inf,inf,1,inf]);
        ut = squeeze(ut);
        ut = double(ut);
        vt = ncread(file,'v',[1,1,60,1],[inf,inf,1,inf]);
        vt = squeeze(vt);
        vt = double(vt);
        % Read in the time array
        time = ncread(file,'time');
        % Store arrays
        disp('Storing U, V and time')
        TIME = cat(1,TIME,time);
        u = cat(3,u,ut);
        v = cat(3,v,vt);
    else
        disp(strcat('No data for',file))
    end
    end
end
    
clear ut vt time

% Clean-up

[~, ia, ~] = unique(TIME);
u = u(:,:,ia);
v = v(:,:,ia);
TIME = TIME(ia);

% Get U and V to be the same grid

for i = 1:size(u,3)
    V(:,:,i) = v2rho_2d(squeeze(v(:,:,i))')'.*mask;
    U(:,:,i) = u2rho_2d(squeeze(u(:,:,i))')'.*mask;
end
 
u = U;
v = V;

clear U V

% Up to here, we have the velocity data over the same time domain of the
% FTLE calculated field

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do the edge data now
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% First need to index the edge map

start = datetime('01-Jan-2014 12:00:00');
back = datetime('31-Dec-2014 12:00:00');

ind_start = find(start==CROCO_time);
ind_back = find(back == CROCO_time);

EDGE = pickle_data(:,:,ind_start:ind_back); % Now we have the same time frame
EDGE = permute(EDGE,[2 1 3]);
clear pickle_data

% Now the EDGEs are on the same time domain as the FTLE calculation
%%  Needs to be seperate, now need to be spatially indexed to same domain as FTLE

% Index the spatial region to the FTLE grid

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

% Velocity data

u = u(idx_lon_min:idx_lon_max,idx_lat_min:idx_lat_max,:).*mask;
v = v(idx_lon_min:idx_lon_max,idx_lat_min:idx_lat_max,:).*mask;

% Edge data

% Need to Patch to make up the size of the domain

EDGE(end+1,:,:) = EDGE(end,:,:);  % Basic patch to get array same size
EDGE(:,end+1,:) = EDGE(:,end,:);

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

u = u(idx_lon_min:idx_lon_max,idx_lat_min:idx_lat_max,:).*mask;
v = v(idx_lon_min:idx_lon_max,idx_lat_min:idx_lat_max,:).*mask;
EDGE = EDGE(idx_lon_min:idx_lon_max,idx_lat_min:idx_lat_max,:);

%% Now get the same time-frames

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Average for summer and winter FTLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load FTLE_zoom_nonorm_2014.mat

% No normalisation has been done yet on the FTLE field to identify LCS
% Need to loop over 

for x = 1:size(myFTLE,3)
    tmp = myFTLE(:,:,x);
    FTLE_norm(:,:,x) = tmp./max(tmp(:));
end

myFTLE = FTLE_norm;
clear FTLE_norm

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

FTLE_SUM = myFTLE(:,:,ind_sum);
sum_date = mydate(1,ind_sum);
FTLE_WIN = myFTLE(:,:,ind_win);
win_date = mydate(1,ind_win);

clear myFTLE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Velcoity + EDGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TIME = TIME(ind_start:ind_back); 
CROCO_time = datetime(Yorig,1,1) + seconds(TIME);
date_CROCO = datenum(CROCO_time);
date_CROCO = datetime(date_CROCO, 'ConvertFrom', 'datenum', 'Format', 'dd-MM-yyyy');

% Summer
for j = 1:length(sum_date)
    dt = abs((datenum(date_CROCO)-datenum(sum_date(j)))); 
    [val,idx] = min(dt);
    sum_IDX(:,j) = idx;
end

% Winter
for j = 1:length(win_date)
    dt = abs((datenum(date_CROCO)-datenum(win_date(j)))); 
    [val,idx] = min(dt);
    win_IDX(:,j) = idx;
end

% Get the corresponding first day of the FTLE computation

% Velocity and EDGE
% summer
u_sum = u(:,:,sum_IDX);
v_sum = v(:,:,sum_IDX);
EDGE_sum = EDGE(:,:,sum_IDX);
date_sum = CROCO_time(sum_IDX);
%winter
u_win = u(:,:,win_IDX);
v_win = v(:,:,win_IDX);
EDGE_win = EDGE(:,:,win_IDX);
date_win = CROCO_time(win_IDX);

%%
% Create some composite plots of some summer and winter data points:
% Note, must offset the u, v and edge by one relative to ftle.

index = [1,23,33,2,23,35];
domain = [15,19;-33.5,-28];

figure
for i = 1:length(index)
    subplot(2,3,i)
    cmin = 0.5;
    cmax = 0.9;
    mydepths = [200,300,500];
    m_proj('miller','long',[domain(1,:)]...
        ,'lat',[domain(2,:)]);
    if i <= 3
        m_pcolor(CROCO_lon,CROCO_lat,FTLE_SUM(:,:,index(i))');
    else
        m_pcolor(CROCO_lon,CROCO_lat,FTLE_WIN(:,:,index(i))');
    end
    shading interp
    m_gshhs_h('patch',[.7 .7 .7],'edgecolor','none');
    m_grid('box','fancy','linest','none','tickdir','out','fontsize',13);
    colormap(flipud(hot))
    caxis([cmin cmax])
    cRange=caxis;
    hold on
    m_contour(CROCO_lon,CROCO_lat,CROCO_top,mydepths,'Color','g','LineWidth',2.5,'LineStyle','--');
    caxis(cRange)
    hold on
    scale_factor = 1.5;
    skp = 50;
    if i <= 3
        ugos_plot = u_sum(:,:,index(i));
        vgos_plot = v_sum(:,:,index(i));
    else
        ugos_plot = u_win(:,:,index(i));
        vgos_plot = v_win(:,:,index(i));
    end
    qp2 = m_quiver(CROCO_lon(1:skp:end),CROCO_lat(1:skp:end), ugos_plot(1:skp:end)...
         ,vgos_plot(1:skp:end),'k');
    set(qp2, 'AutoScale','on','AutoScaleFactor',scale_factor);
    hold on
    if i <= 3
        m_contourf(CROCO_lon,CROCO_lat,EDGE_sum(:,:,index(i)),'Color','b','LineWidth',0.001,'LineStyle',':');
    else
        m_contourf(CROCO_lon,CROCO_lat,EDGE_win(:,:,index(i)),'Color', 'b','LineWidth',0.001,'LineStyle',':');
    end
    caxis(cRange)
    hold on
    if i <=3
        tmp = convertStringsToChars(string(date_sum(index(i))));
        tmp = convertCharsToStrings(tmp(1:11));
        title(strcat(tmp),'FontSize', 14)
    else
        tmp = convertStringsToChars(string(date_win(index(i))));
        tmp = convertCharsToStrings(tmp(1:11));
        title(strcat(tmp),'FontSize', 14)
    end
end

% Create colorbar
hold on
ca = colorbar('southoutside');
ca.Position = ca.Position + 1e-10;
ca.Label.String = 'FTLE';
ca.FontSize = 12;
caxis([cmin cmax]);

% Save figure
% Save the figure

set(gcf, 'InvertHardcopy', 'off')
print('-f1','FTLE_snap','-dpng','-r600');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Specific figure with edges overlaid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We want to show a summer and winter example with the FTLE field and then overliad 
% the frontal positions overlaid to highlight the point

% Summer

index = [1,55,1,55];
label = ["(a)","(b)","(c)","(d)"];

figure
for i = 1:length(index)
    subplot(2,2,i)
    cmin = 0.5;
    cmax = 0.8;
    mydepths = [200,300,500];
    m_proj('miller','long',[domain(1,:)]...
        ,'lat',[domain(2,:)]);
    m_pcolor(CROCO_lon,CROCO_lat,FTLE_SUM(:,:,index(i))');
    shading interp
    m_gshhs_h('patch',[.7 .7 .7],'edgecolor','none');
    m_grid('box','fancy','linest','none','tickdir','out','fontsize',13);
    colormap(flipud(hot))
    caxis([cmin cmax])
    cRange=caxis;
    hold on
    m_contour(CROCO_lon,CROCO_lat,CROCO_top,mydepths,'Color','g','LineWidth',2.5,'LineStyle','--');
    caxis(cRange)
    hold on
    scale_factor = 1.5;
    skp = 50;
    ugos_plot = u_sum(:,:,index(i));
    vgos_plot = v_sum(:,:,index(i));
    qp2 = m_quiver(CROCO_lon(1:skp:end),CROCO_lat(1:skp:end), ugos_plot(1:skp:end)...
         ,vgos_plot(1:skp:end),'k');
    set(qp2, 'AutoScale','on','AutoScaleFactor',scale_factor);
    hold on
    if i > 2
        m_contour(CROCO_lon,CROCO_lat,EDGE_sum(:,:,index(i)),'Color','b','LineWidth',1,'LineStyle','-');
    end
    caxis(cRange)
    hold on
    if i <= 2
        tmp = convertStringsToChars(string(date_sum(index(i))));
        tmp = convertCharsToStrings(tmp(1:11));
        title(strcat(label(i)," ", tmp),'FontSize', 14)
    else
        title(strcat(label(i)," ","Frontal edges"),'FontSize', 14)
    end
end

% Create colorbar
hold on
ca = colorbar('southoutside');
ca.Position = ca.Position + 1e-10;
ca.Label.String = 'FTLE';
ca.FontSize = 12;
caxis([cmin cmax]);

set(gcf, 'InvertHardcopy', 'off')
print('-f1','FTLE_EDGE_summer','-dpng','-r600');

fig1 = gcf;   % This is the position of the figure after adjusting it

% Winter

index = [2,35,2,35];
label = ["(a)","(b)","(c)","(d)"];

figure
for i = 1:length(index)
    subplot(2,2,i)
    cmin = 0.5;
    cmax = 0.8;
    mydepths = [200,300,500];
    m_proj('miller','long',[domain(1,:)]...
        ,'lat',[domain(2,:)]);
    m_pcolor(CROCO_lon,CROCO_lat,FTLE_WIN(:,:,index(i))');
    shading interp
    m_gshhs_h('patch',[.7 .7 .7],'edgecolor','none');
    m_grid('box','fancy','linest','none','tickdir','out','fontsize',13);
    colormap(flipud(hot))
    caxis([cmin cmax])
    cRange=caxis;
    hold on
    m_contour(CROCO_lon,CROCO_lat,CROCO_top,mydepths,'Color','g','LineWidth',2.5,'LineStyle','--');
    caxis(cRange)
    hold on
    scale_factor = 1.5;
    skp = 50;
    ugos_plot = u_win(:,:,index(i));
    vgos_plot = v_win(:,:,index(i));
    qp2 = m_quiver(CROCO_lon(1:skp:end),CROCO_lat(1:skp:end), ugos_plot(1:skp:end)...
         ,vgos_plot(1:skp:end),'k');
    set(qp2, 'AutoScale','on','AutoScaleFactor',scale_factor);
    hold on
    if i > 2
        m_contour(CROCO_lon,CROCO_lat,EDGE_win(:,:,index(i)),'Color','b','LineWidth',1,'LineStyle','-');
    end
    caxis(cRange)
    hold on
    if i <= 2
        tmp = convertStringsToChars(string(date_win(index(i))));
        tmp = convertCharsToStrings(tmp(1:11));
        title(strcat(label(i)," ", tmp),'FontSize', 14)
    else
        title(strcat(label(i)," ","Frontal edges"),'FontSize', 14)
    end
end

% Create colorbar
hold on
ca = colorbar('southoutside');
ca.Position = ca.Position + 1e-10;
ca.Label.String = 'FTLE';
ca.FontSize = 12;
caxis([cmin cmax]);

set(gcf,'Position',fig1.Position)

set(gcf, 'InvertHardcopy', 'off')
print('-f1','FTLE_EDGE_win','-dpng','-r600');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cmin = 0;
cmax = 0.8;
levels = [cmin:0.1:cmax];
mydepths = [200,300,500];

figure
m_proj('miller','long',[domain(1,:)]...
    ,'lat',[domain(2,:)]);
subplot(1,2,1);
m_contourf(CROCO_lon,CROCO_lat,nanmean(FTLE_SUM,3)',levels,'ShowText','off');
m_gshhs_h('patch',[.7 .7 .7],'edgecolor','none');
m_grid('box','fancy','linest','none','tickdir','out','fontsize',13);
%colormap(flipud('gray'))
cmocean('curl',length(levels))
caxis([cmin cmax])
cRange=caxis;
hold on
m_contour(CROCO_lon,CROCO_lat,CROCO_top,mydepths,'Color','g','LineWidth',1,'LineStyle','--');
caxis(cRange)
title("CROCO [DJF]",'FontSize', 14)
subplot(1,2,2)
m_contourf(CROCO_lon,CROCO_lat,nanmean(FTLE_WIN,3)',levels,'ShowText','off');
m_gshhs_h('patch',[.7 .7 .7],'edgecolor','none');
m_grid('box','fancy','linest','none','tickdir','out','fontsize',13);
cmocean('curl',length(levels))
caxis([cmin cmax])
cRange=caxis;
hold on
m_contour(CROCO_lon,CROCO_lat,CROCO_top,mydepths,'Color','g','LineWidth',1,'LineStyle','--');
caxis(cRange)
title("CROCO [JJA]",'FontSize', 14)

% Create colorbar
hold on
ca = colorbar('eastOutside');
ca.Position = ca.Position + 1e-10;
%ca.Label.String = 'FTLE (days)';
ca.FontSize = 12;
caxis([cmin cmax]);

% Save the figure

set(gcf, 'InvertHardcopy', 'off')
print('-f1','TT','-dpng','-r600');