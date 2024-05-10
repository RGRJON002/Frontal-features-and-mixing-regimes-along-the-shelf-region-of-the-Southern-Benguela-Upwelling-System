%%%% Composite plot overlaying the streamfunction transport, vertical
%%%% velocities and frontal edges for summer and winter UF and SBF,
%%%% Restrict domain to only that north of Cape Columbine

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
addpath /home/jono/CROCO/croco_tools/Diagnostic_tools/EKE
addpath '/media/data/DIAGNOSTICS/FRONTS'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

mask_plot = mask;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get mask for shelf region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Nearshore
shelf_depth = 500;
h_mask = CROCO_top;
h_mask(h_mask > shelf_depth) = NaN;
mask_canny = ~isnan(h_mask);
mask_canny = double(mask_canny); 
mask_canny(mask_canny ==0) = nan;

% Shelf-break front

sb_mask = mask_canny;
sb_mask(h_mask < 200)  = NaN;
sb_mask(h_mask > 500) = NaN;

% Upwelling front

up_mask = mask_canny; 
up_mask(h_mask > 200) = NaN;

% Correct for slight mismatch in topography

sb_mask(26:35,221:233) = 1;
up_mask(26:35,221:233) = NaN;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDGE DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% Loop over the lat dimension and find the max value

EDGE_sum = nansum(EDGE(:,:,ind_sum),3);
EDGE_win = nansum(EDGE(:,:,ind_win),3);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Positions of the fronts
%%%%%%%%%%%%%%%%%%%%%%%%%%

% UF

% Upwelling region
UF_sum = zeros(180,341);
UF_win = zeros(size(UF_sum));
for i = 1:size(up_mask,2)
    % Summer
    tmp = EDGE_sum.*up_mask;
    ind = find(max(tmp(:,i)) == tmp(:,i));
    if length(ind) > 1
        % Will take the mean of the max's 
        ind = max(ind);
        UF_sum(ind,i) = 1;
    else
        UF_sum(ind,i) = 1;
    end
    % Winter
    tmp = EDGE_win.*up_mask;
    ind = find(max(tmp(:,i)) == tmp(:,i));
    if length(ind) > 1
        % Will take the mean of the max's 
        ind = max(ind);
        UF_win(ind,i) = 1;
    else
        UF_win(ind,i) = 1;
    end
end

% SBF

% Upwelling region
SBF_sum = zeros(180,341);
SBF_win = zeros(size(SBF_sum));
for i = 1:size(up_mask,2)
    % Summer
    tmp = EDGE_sum.*sb_mask;
    ind = find(max(tmp(:,i)) == tmp(:,i));
    if length(ind) > 1
        % Will take the mean of the max's 
        ind = max(ind);
        SBF_sum(ind,i) = 1;
    else
        SBF_sum(ind,i) = 1;
    end
    % Winter
    tmp = EDGE_win.*sb_mask;
    ind = find(max(tmp(:,i)) == tmp(:,i));
    if length(ind) > 1
        % Will take the mean of the max's 
        ind = max(ind);
        SBF_win(ind,i) = 1;
    else
        SBF_win(ind,i) = 1;
    end
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STREAM DATA                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

start   % Have to tun to get libraries for STREAM

% Summer

load('summer_transport_croco.mat'); 

% Average all variable data for our seasons

% Index U, V and psi

U_sum = double(nanmean(U_sum,3));
V_sum = double(nanmean(V_sum,3));
PSI_sum = double(nanmean(PSI_sum,3));

U_sum = u2rho_2d(U_sum)';
V_sum = v2rho_2d(V_sum)';
PSI_sum = 1e-6*psi2rho(PSI_sum)';

U_sum = U_sum(idx_lon_min:idx_lon_max,idx_lat_min:idx_lat_max);
V_sum = V_sum(idx_lon_min:idx_lon_max,idx_lat_min:idx_lat_max);
PSI_sum = PSI_sum(idx_lon_min:idx_lon_max,idx_lat_min:idx_lat_max);

% Winter

load('winter_transport_croco.mat'); 

% Average all variable data for our seasons

% Index U, V and psi

U_win = double(nanmean(U_win,3));
V_win = double(nanmean(V_win,3));
PSI_win = double(nanmean(PSI_win,3));

U_win = u2rho_2d(U_win)';
V_win = v2rho_2d(V_win)';
PSI_win = 1e-6*psi2rho(PSI_win)';

U_win = U_win(idx_lon_min:idx_lon_max,idx_lat_min:idx_lat_max);
V_win = V_win(idx_lon_min:idx_lon_max,idx_lat_min:idx_lat_max);
PSI_win = PSI_win(idx_lon_min:idx_lon_max,idx_lat_min:idx_lat_max);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vertical velocity data          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load W_1.mat
% load W_2.mat
% load W_3.mat
% load W_4.mat
% 
% W = cat(3,W_1,W_2,W_3,W_4);
% 
% clear W_1 W_2 W_3 W_4
% 
% W_sum = nanmean(W(:,:,ind_sum),3)';
% W_win = nanmean(W(:,:,ind_win),3)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lat_min = -33.5;
lat_max = -28;
lon_min = 15;
lon_max = 20;

lon = CROCO_lon;
lat = CROCO_lat;

figure
npts=[2 2 2 2];
topo=1;

subplot(1,2,1)  % Summer

m_proj('mercator',...
       'lon',[lon_min lon_max],...
       'lat',[lat_min lat_max]);
[x,y]=m_ll2xy(double(lon),double(lat),'clip','off');

hold on

psi_r= PSI_sum;
u = U_sum;
v = V_sum;

colmax = 22;
colmin = 6;
dcol = 1;

[C1,h1]=m_contour(rempoints(lon,npts),rempoints(lat,npts),rempoints(psi_r,npts)...
           ,[dcol:dcol:colmax],'k');
if ~isempty(C1)
 % clabel(C1,h1,'LabelSpacing',1000,'Rotation',0)
  hf1=add_streamarrows(C1,x',y',u',v');
  hold on
end

[C2,h2]=m_contour(rempoints(lon,npts),rempoints(lat,npts),rempoints(-psi_r,npts)...
           ,[dcol:dcol:colmin],'k');
if ~isempty(C2)
 % clabel(C2,h2,'LabelSpacing',1000,'Rotation',0)
  hf2=add_streamarrows(C2,x',y',u',v');
  hold on
end

[C3,h3]=m_contour(rempoints(lon,npts),rempoints(lat,npts),rempoints(psi_r,npts)...
           ,[0 0],'k');
if ~isempty(C3)
 % clabel(C3,h3,'LabelSpacing',1000,'Rotation',0)
  hf3=add_streamarrows(C3,x',y',u',v');
  set(h3,'Color','k','Linewidth',1.2);
  set(hf3,'Linewidth',1.2);
  hold on
end

[C4,h4]=m_contour(rempoints(lon,npts),rempoints(lat,npts),rempoints(psi_r,npts)...
           ,[20:0.1:21],'k');
if ~isempty(C4)
 % clabel(C4,h4,'LabelSpacing',1000,'Rotation',0)
  hf4=add_streamarrows(C4,x',y',u',v');
  set(h4,'Color','k','Linewidth',1.2);
  set(hf4,'Linewidth',1.2);
  hold on
end

if topo==1
  [C5,h5]=m_contour(rempoints(lon,npts),rempoints(lat,npts),rempoints(CROCO_top,npts)...
           ,[200 300 500],'--g');
  %clabel(C5,h5,'LabelSpacing',1000,'Rotation',0)
end

m_usercoast('coastline_l.mat','patch',[.9 .9 .9]);
hold off
m_gshhs_h('patch',[.8 .8 .8],'edgecolor','none');
m_grid('box','fancy','xtick',5,'ytick',5,'tickdir','out');
set(findobj('tag','m_grid_color'),'facecolor','white')
hold on
m_contour(rempoints(lon,npts),rempoints(lat,npts),rempoints(UF_sum,npts)...
           ,'-r');
hold on
m_contour(rempoints(lon,npts),rempoints(lat,npts),rempoints(SBF_sum,npts)...
           ,'-b');
title(['(a) ','Streamlines [DJF]'],'FontSize',14)
hold on

% Add cross-section locations

str = {'(1)','(2)','(3)'};

lon_tmp = CROCO_lon(:,1);
lat_tmp = CROCO_lat(1,:);

locs = [-32,-31.1,-29];

for i = 1:length(locs)
    lat_m = locs(i);
    LTmin_index = min(find(lat_tmp>=lat_m));
    m_line([lon_tmp(min(find(CROCO_top(:,LTmin_index)<=1000))) lon_tmp(min(find(isnan(mask_plot(:,LTmin_index)))))]...
        ,[lat_tmp(LTmin_index) lat_tmp(LTmin_index)],'linewi',2,'color','magenta');
    m_text(double(lon_tmp(min(find(isnan(mask_plot(:,LTmin_index)))))+0.01),double(lat_m),str{i},'Fontsize',15)
    hold on
end

% Custom legend

qw{1} = m_plot(nan,nan, 'r--','LineWidth',2);
qw{2} = m_plot(nan,nan ,'b--','LineWidth',2);
qw{3} = m_plot(nan,nan ,'m-','LineWidth',2);
l = legend([qw{:}], {'UF','SBF','Section'}, 'location', 'best');
l.FontSize = 12;

hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,2,2)  % Winter

m_proj('mercator',...
       'lon',[lon_min lon_max],...
       'lat',[lat_min lat_max]);
[x,y]=m_ll2xy(double(lon),double(lat),'clip','off');

hold on

psi_r= PSI_win;
u = U_win;
v = V_win;

colmax = 22;
colmin = 6;
dcol = 1;

[C1,h1]=m_contour(rempoints(lon,npts),rempoints(lat,npts),rempoints(psi_r,npts)...
           ,[dcol:dcol:colmax],'k');
if ~isempty(C1)
 % clabel(C1,h1,'LabelSpacing',1000,'Rotation',0)
  hf1=add_streamarrows(C1,x',y',u',v');
  hold on
end

[C2,h2]=m_contour(rempoints(lon,npts),rempoints(lat,npts),rempoints(-psi_r,npts)...
           ,[dcol:dcol:colmin],'k');
if ~isempty(C2)
 % clabel(C2,h2,'LabelSpacing',1000,'Rotation',0)
  hf2=add_streamarrows(C2,x',y',u',v');
  hold on
end

[C3,h3]=m_contour(rempoints(lon,npts),rempoints(lat,npts),rempoints(psi_r,npts)...
           ,[0 0],'k');
if ~isempty(C3)
 % clabel(C3,h3,'LabelSpacing',1000,'Rotation',0)
  hf3=add_streamarrows(C3,x',y',u',v');
  set(h3,'Color','k','Linewidth',1.2);
  set(hf3,'Linewidth',1.2);
  hold on
end

[C4,h4]=m_contour(rempoints(lon,npts),rempoints(lat,npts),rempoints(psi_r,npts)...
           ,[20:0.1:21],'k');
if ~isempty(C4)
 % clabel(C4,h4,'LabelSpacing',1000,'Rotation',0)
  hf4=add_streamarrows(C4,x',y',u',v');
  set(h4,'Color','k','Linewidth',1.2);
  set(hf4,'Linewidth',1.2);
  hold on
end

if topo==1
  [C5,h5]=m_contour(rempoints(lon,npts),rempoints(lat,npts),rempoints(CROCO_top,npts)...
           ,[200 300 500],'--g');
  %clabel(C5,h5,'LabelSpacing',1000,'Rotation',0)
end

m_usercoast('coastline_l.mat','patch',[.9 .9 .9]);
hold off
m_gshhs_h('patch',[.8 .8 .8],'edgecolor','none');
m_grid('box','fancy','xtick',5,'ytick',5,'tickdir','out');
set(findobj('tag','m_grid_color'),'facecolor','white')

hold on
m_contour(rempoints(lon,npts),rempoints(lat,npts),rempoints(UF_win,npts)...
           ,'-r');
hold on
m_contour(rempoints(lon,npts),rempoints(lat,npts),rempoints(SBF_win,npts)...
           ,'-b');
title(['(b) ','Streamlines [JJA]'],'FontSize',14)
hold on

% Add labels and section

str = {'(1)','(2)','(3)'};

lon_tmp = CROCO_lon(:,1);
lat_tmp = CROCO_lat(1,:);

locs = [-32,-31.1,-29];

for i = 1:length(locs)
    lat_m = locs(i);
    LTmin_index = min(find(lat_tmp>=lat_m));
    m_line([lon_tmp(min(find(CROCO_top(:,LTmin_index)<=1000))) lon_tmp(min(find(isnan(mask_plot(:,LTmin_index)))))]...
        ,[lat_tmp(LTmin_index) lat_tmp(LTmin_index)],'linewi',2,'color','magenta');
    m_text(double(lon_tmp(min(find(isnan(mask_plot(:,LTmin_index)))))+0.01),double(lat_m),str{i},'Fontsize',15)
    hold on
end
hold off

% Save the figure

set(gcf, 'InvertHardcopy', 'off')
print('-f1','STREAM_EDGE','-dpng','-r600');
