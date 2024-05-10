%%%%% Look at wind-stress and curl as well as vertical velocities

addpath /home/jono/CROCO/croco_tools/UTILITIES/m_map1.4h
addpath /media/data/DATASETS_CROCOTOOLS/m_map1.4f
addpath /home/jono/Documents/MATLAB/CMOCEAN
addpath '/home/jono/Documents/MATLAB/CDT-master/cdt'
addpath '/home/jono/Documents/MATLAB/CDT-master/cdt/cdt_data'
addpath '/home/jono/Documents/MATLAB/tight_subplot'
addpath('/usr/local/MATLAB/R2020a/toolbox/matlab/imagesci/')

CROCO_path = '/media/data/CHPC_SBUS_3km/';
Ymin = 2005;
Ymax = 2010;
Yorig = 1990;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%
% CROCO
%%%%%%%%%

% Once-off gridding data and dimensions
file = strcat(CROCO_path,'/','croco_avg_Y',string(Ymin),'M',string(1),'.nc');
% Read in the mask
mask=ncread(file,'mask_rho');
mask(mask==0)=nan;
% Read in topo
h = ncread(file,'h');
h = h.*mask;
% Define my region, will focus on the coastal domain
CROCO_lat=ncread(strcat(CROCO_path,'/','croco_avg_Y',string(Ymin),'M',string(1),'.nc'),'lat_rho');
CROCO_lon=ncread(strcat(CROCO_path,'/','croco_avg_Y',string(Ymin),'M',string(1),'.nc'),'lon_rho');
% Only need to index lon as lat bounds and eastern extent stay the same

%lon_min = 14;
%[~,idx_lon_min]=min(abs(CROCO_lon(:,1)-lon_min));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the U and V stress

U_stress = [];
V_stress = [];
TIME = [];
for i = Ymin:Ymax
    for j = 1:12
    file = strcat(CROCO_path,'/','croco_avg_Y',string(i),'M',string(j),'.nc');
    if exist(file, 'file')
        disp(file)
        disp('Reading data')
        % Read in the time array
        time = ncread(file,'time');
        % Read in the wind-stress
        u_stress = ncread(file,'sustr',[idx_lon_min 1 1],[inf inf inf]);
        v_stress = ncread(file,'svstr',[idx_lon_min 1 1],[inf inf inf]);
        % Fix the format
        u_stress = (u_stress(:,2:end,:)+u_stress(:,1:end-1,:))/2; 
        v_stress = (v_stress(2:end,:,:)+v_stress(1:end-1,:,:))/2;
        %  Store arrays
        disp('Storing U_stress, V_stress and time')
        TIME = cat(1,TIME,time);
        U_stress = cat(3,U_stress,u_stress);
        V_stress = cat(3,V_stress,v_stress);
    else
        disp(strcat('No data for',file))
    end
    clear u_stress v_stress time 
    end
end

W = [];
for i = Ymin:Ymax
    for j = 1:12
    file = strcat(CROCO_path,'/','croco_avg_Y',string(i),'M',string(j),'.nc');
    if exist(file, 'file')
        disp(file)
        disp('Reading data')
        % Get the vertical grid for that month
        grd = getmydepth(convertStringsToChars(file),200);
        disp('Processing grid')
        %grd = grd(idx_lon_min:end,:,:,:);
        % Read in the vertical velocity data
        addpath('/usr/local/MATLAB/R2020a/toolbox/matlab/imagesci/')
        %w = ncread(file,'w',[idx_lon_min 1 1 1],[inf inf inf inf]);
        w = ncread(file,'w');
        w = bsxfun(@rdivide, squeeze(nansum(w.*grd,3)), h);
        % Fix the format
        %w = (w(:,2:end,:)+w(:,1:end-1,:))/2; 
        %w = (w(2:end,:,:)+w(1:end-1,:,:))/2;
        %  Store arrays
        disp('Storing w and time')
        W = cat(3,W,w);
    else
        disp(strcat('No data for',file))
    end
    clear w 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POST-PROCESSING  GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define my region, will focus on the coastal domain
croco_top = ncread(strcat(CROCO_path,'/','croco_avg_Y',string(Ymin),'M',string(1),'.nc'),'h',[idx_lon_min 1 ],[inf inf]);
mask = mask(idx_lon_min:end,:);
croco_top = croco_top.*mask;

% Edit to new grid

CROCO_lon = CROCO_lon(idx_lon_min:end,:);
CROCO_lat = CROCO_lat(idx_lon_min:end,:);

new_lat=(CROCO_lat(1,2:end)+CROCO_lat(1,1:end-1))/2;
new_lon=(CROCO_lon(2:end,1)+CROCO_lon(1:end-1,1))/2;
[CROCO_lon,CROCO_lat] = meshgrid(new_lon,new_lat);
CROCO_lon=CROCO_lon';
CROCO_lat=CROCO_lat';
clear new_lon new_lat

croco_top = (croco_top(2:end,:,:)+croco_top(1:end-1,:,:))/2;
croco_top = (croco_top(:,2:end,:)+croco_top(:,1:end-1,:))/2; 

mask = (mask(2:end,:,:)+mask(1:end-1,:,:))/2;
mask = (mask(:,2:end,:)+mask(:,1:end-1,:))/2;

lat_min = -36;
lat_max = -28;
lon_min = 15;
lon_max = 20;

[~,idx_lat_min]=min(abs(CROCO_lat(1,:)-lat_min));
[~,idx_lat_max]=min(abs(CROCO_lat(1,:)-lat_max));
[~,idx_lon_min]=min(abs(CROCO_lon(:,1)-lon_min));
[~,idx_lon_max]=min(abs(CROCO_lon(:,1)-lon_max));

%%%%%%%%%%%%%%%%%%%%%%
% POST-PROCESS DATA
%%%%%%%%%%%%%%%%%%%%%%

% Convert the time data

CROCO_time = datetime(Yorig,1,1) + seconds(TIME);
[Y,MO,D] = datevec(CROCO_time);  
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

% Subset the data and find seasonal averages

w_summer = double(nanmean(W(:,:,ind_sum),3)).*mask;
w_winter = double(nanmean(W(:,:,ind_win),3)).*mask;

Taux_sum = double(nanmean(U_stress(:,:,ind_sum),3)).*mask;
Tauy_sum = double(nanmean(V_stress(:,:,ind_sum),3)).*mask;

Taux_win = double(nanmean(U_stress(:,:,ind_win),3)).*mask;
Tauy_win = double(nanmean(V_stress(:,:,ind_win),3)).*mask;

Curl_sum = cdtcurl(CROCO_lon, CROCO_lat, Taux_sum, Tauy_sum); 
Curl_win = cdtcurl(CROCO_lon, CROCO_lat, Taux_win, Tauy_win); 

%Curl_sum = Curl_sum.*sign(CROCO_lat);
%Curl_win = Curl_win.*sign(CROCO_lat);

%%%% First we plot the season

cmin = -10e-7;
cmax = 10e-7;
scale_factor = 1.5;
mydepths = [200,300,500];
skp = 10;

figure

m_proj('miller','long',[double(CROCO_lon(idx_lon_min,1)) double(CROCO_lon(idx_lon_max,1))]...
    ,'lat',[double(CROCO_lat(1,idx_lat_min)) double(CROCO_lat(1,idx_lat_max))]);
subplot(1,2,1)
m_pcolor(CROCO_lon,CROCO_lat,Curl_sum);
shading flat
cmocean('balance')
m_gshhs_h('patch',[.8 .8 .8],'edgecolor','none');
m_grid('box','fancy','linest','none','tickdir','out','fontsize',13);
caxis([cmin cmax])
cRange=caxis;
hold on
qp1 = m_quiver(double(CROCO_lon(1:skp:end,1:skp:end)),double(CROCO_lat(1:skp:end,1:skp:end)),...
    double(Taux_sum(1:skp:end,1:skp:end)),...
    double(Tauy_sum(1:skp:end,1:skp:end)),'Color',[0,0,0]);
set(qp1, 'AutoScale','on','AutoScaleFactor',scale_factor);
hold on
m_contour(CROCO_lon,CROCO_lat,croco_top,mydepths,'Color','y','LineWidth',1,'LineStyle','--');
caxis(cRange)
subplot(1,2,2)
m_pcolor(CROCO_lon,CROCO_lat,Curl_win);
shading flat
cmocean('balance')
m_gshhs_f('patch',[.8 .8 .8],'edgecolor','none');
m_grid('box','fancy','linest','none','tickdir','out','fontsize',13);
caxis([cmin cmax])
cRange=caxis;
hold on
qp1 = m_quiver(double(CROCO_lon(1:skp:end,1:skp:end)),double(CROCO_lat(1:skp:end,1:skp:end)),...
    double(Taux_win(1:skp:end,1:skp:end)),...
    double(Tauy_win(1:skp:end,1:skp:end)),'Color',[0,0,0]);
set(qp1, 'AutoScale','on','AutoScaleFactor',scale_factor);
hold on
m_contour(CROCO_lon,CROCO_lat,croco_top,mydepths,'Color','y','LineWidth',1,'LineStyle','--');
caxis(cRange)

% Create colorbar
hold on
ca = colorbar('eastOutside');
ca.Position = ca.Position + 1e-10;
ca.Label.String = 'Wind-stress curl (N/m^3)';
ca.FontSize = 12;
caxis([cmin cmax]);

% Save the figure

set(gcf, 'InvertHardcopy', 'off')
print('-f1','Wind_stress','-dpng','-r600');

% Plot verticl velocities overlayed with wind stress curl

cmin = -3e-4;
cmax = 3e-4;
scale_factor = 1.5;
mydepths = [200,300,500];
skp = 10;

figure

m_proj('miller','long',[double(CROCO_lon(idx_lon_min,1)) double(CROCO_lon(idx_lon_max,1))]...
    ,'lat',[double(CROCO_lat(1,idx_lat_min)) double(CROCO_lat(1,idx_lat_max))]);
subplot(1,2,1)
m_pcolor(CROCO_lon,CROCO_lat,w_summer);
shading flat
cmocean('balance')
m_gshhs_h('patch',[.8 .8 .8],'edgecolor','none');
m_grid('box','fancy','linest','none','tickdir','out','fontsize',13);
caxis([cmin cmax])
cRange=caxis;
hold on
qp1 = m_quiver(double(CROCO_lon(1:skp:end,1:skp:end)),double(CROCO_lat(1:skp:end,1:skp:end)),...
    double(Taux_sum(1:skp:end,1:skp:end)),...
    double(Tauy_sum(1:skp:end,1:skp:end)),'Color',[0,0,0]);
set(qp1, 'AutoScale','on','AutoScaleFactor',scale_factor);
hold on
m_contour(CROCO_lon,CROCO_lat,croco_top,mydepths,'Color','k','LineWidth',1,'LineStyle','--');
caxis(cRange)
subplot(1,2,2)
m_pcolor(CROCO_lon,CROCO_lat,w_winter);
shading flat
cmocean('balance')
m_gshhs_f('patch',[.8 .8 .8],'edgecolor','none');
m_grid('box','fancy','linest','none','tickdir','out','fontsize',13);
caxis([cmin cmax])
cRange=caxis;
hold on
qp1 = m_quiver(double(CROCO_lon(1:skp:end,1:skp:end)),double(CROCO_lat(1:skp:end,1:skp:end)),...
    double(Taux_win(1:skp:end,1:skp:end)),...
    double(Tauy_win(1:skp:end,1:skp:end)),'Color',[0,0,0]);
set(qp1, 'AutoScale','on','AutoScaleFactor',scale_factor);
hold on
m_contour(CROCO_lon,CROCO_lat,croco_top,mydepths,'Color','k','LineWidth',1,'LineStyle','--');
caxis(cRange)

% Create colorbar
hold on
ca = colorbar('eastOutside');
ca.Position = ca.Position + 1e-10;
ca.Label.String = 'Vertical velocities (m/s)';
ca.FontSize = 12;
caxis([cmin cmax]);

set(gcf, 'InvertHardcopy', 'off')
print('-f1','W_plot','-dpng','-r600');