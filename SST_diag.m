%%%% This script will plot the SST, U, V data all as composites as well as
%%%% create animations

addpath /home/jono/CROCO/croco_tools/UTILITIES/m_map1.4h
addpath /media/data/DATASETS_CROCOTOOLS/m_map1.4f
addpath /home/jono/Documents/MATLAB/CMOCEAN
addpath '/home/jono/Documents/MATLAB/CDT-master/cdt'
addpath '/home/jono/Documents/MATLAB/CDT-master/cdt/cdt_data'
addpath '/home/jono/Documents/MATLAB/tight_subplot'

CROCO_path = '/media/jono/SBUS/SBUS_3km/CHPC_OUTPUT';
Ymin = 2004;
Ymax = 2018;
Yorig = 1990;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%
% CROCO
%%%%%%%%%
g = 9.81;
SST = [];
U = [];
V = [];
TIME = [];

for i = Ymin:Ymax
    for j = 1:12
    file = strcat(CROCO_path,'/','croco_avg_Y',string(i),'M',string(j),'.nc');
    if exist(file, 'file')
        disp(file)
        disp('Reading data')
        % Read in the mask
        mask=ncread(file,'mask_rho');
        mask(mask==0)=nan;
        % extract coriolis
        f = ncread(file,'f');
        % extract dx and dy
        dx = ncread(file,'pm');
        dy = ncread(file,'pn');
        ssh = ncread(file,'zeta');
        ssh=ssh.*mask;  % mask out the land
        % Calculate vg
        dn=ssh(2:end,:,:)-ssh(1:end-1,:,:);
        ff = (f(2:end,:)+f(1:end-1,:))/2;
        dxx = (dx(2:end,:)+ dx(1:end-1,:))/2;
        vg = g./ff .* dn.*dxx;
        % Calculate ug
        dn=ssh(:,2:end,:)-ssh(:,1:end-1,:);
        ff = (f(:,2:end)+f(:,1:end-1))/2;
        dyy = (dy(:,2:end)+ dy(:,1:end-1))/2;
        ug = - g./ff .* dn.*dyy;
        %%% Do interpolation to get the same grid
        ug = (ug(2:end,:,:)+ug(1:end-1,:,:))/2;
        vg = (vg(:,2:end,:)+vg(:,1:end-1,:))/2;
        % Read in the surface sst data
        sst = squeeze(ncread(file,'temp',[1 1 60 1],[inf inf 1 inf]));
        sst = sst.*mask;
        sst = (sst(2:end,:,:)+sst(1:end-1,:,:))/2;
        sst = (sst(:,2:end,:)+sst(:,1:end-1,:))/2;        
        % Read in the time array
        time = ncread(file,'time');
        % Store arrays
        disp('Storing SST, U, V and time')
        SST = cat(3,SST,sst);
        TIME = cat(1,TIME,time);
        U = cat(3,U,ug);
        V = cat(3,V,vg);
    else
        disp(strcat('No data for',file))
    end
    end
end
    
% Clean-up 
clear sst ssh dx dy ug vg ff dn dxx dyy time

% Read in the lon, lat and topography data

CROCO_lat=ncread(strcat(CROCO_path,'/','croco_avg_Y',string(Ymin),'M',string(1),'.nc'),'lat_rho');
CROCO_lon=ncread(strcat(CROCO_path,'/','croco_avg_Y',string(Ymin),'M',string(1),'.nc'),'lon_rho');
croco_top = ncread(strcat(CROCO_path,'/','croco_avg_Y',string(Ymin),'M',string(1),'.nc'),'h');

croco_top = croco_top.*mask;

% Edit to new grid

new_lat=(CROCO_lat(1,2:end)+CROCO_lat(1,1:end-1))/2;
new_lon=(CROCO_lon(2:end,1)+CROCO_lon(1:end-1,1))/2;

[CROCO_lon,CROCO_lat] = meshgrid(new_lon,new_lat);
CROCO_lon=CROCO_lon';
CROCO_lat=CROCO_lat';

clear new_lon new_lat

croco_top = (croco_top(2:end,:,:)+croco_top(1:end-1,:,:))/2;
croco_top = (croco_top(:,2:end,:)+croco_top(:,1:end-1,:))/2; 

%%%%%%%%%%%%%%%%%%%%%%
% PROCESS DATA
%%%%%%%%%%%%%%%%%%%%%%

% Convert the time data

CROCO_time = datetime(Yorig,1,1) + seconds(TIME);
[Y,MO,D] = datevec(CROCO_time);    

% Average all sst data for our seasons

summer = [1, 2, 12];  winter = [6, 7, 8]; 

ind_sum = [];
for i = 1:length(summer)
    ind_sum = [ind_sum; find(summer(i) == MO)];
end

ind_win = [];
for i = 1:length(winter)
    ind_win = [ind_win; find(winter(i) == MO)];
end

sst_summer = nanmean(double(SST(:,:,ind_sum)),3);
sst_winter = nanmean(double(SST(:,:,ind_win)),3);

% Do the same for the velocity data

ugos_summer = nanmean(double(U(:,:,ind_sum)),3);
vgos_summer = nanmean(double(V(:,:,ind_sum)),3);

ugos_winter = nanmean(double(U(:,:,ind_win)),3);
vgos_winter = nanmean(double(V(:,:,ind_win)),3);

%%%% First we plot the season

% Define my region, will focus on the coastal domain

lat_min = -36;
lat_max = -28;
lon_min = 15;
lon_max = 20;

[~,idx_lat_min]=min(abs(CROCO_lat(1,:)-lat_min));
[~,idx_lat_max]=min(abs(CROCO_lat(1,:)-lat_max));
[~,idx_lon_min]=min(abs(CROCO_lon(:,1)-lon_min));
[~,idx_lon_max]=min(abs(CROCO_lon(:,1)-lon_max));

cmin = 12;
cmax = 24;
levels = [cmin:0.5:cmax];
scale_factor = 1.5;
mydepths = [200,300,500];
skp = 6;

figure

m_proj('miller','long',[double(CROCO_lon(idx_lon_min,1)) double(CROCO_lon(idx_lon_max,1))]...
    ,'lat',[double(CROCO_lat(1,idx_lat_min)) double(CROCO_lat(1,idx_lat_max))]);
subplot(1,2,1)
m_pcolor(CROCO_lon,CROCO_lat,sst_summer);
shading flat
cmocean('thermal',length(levels))
m_gshhs_f('patch',[.8 .8 .8],'edgecolor','none');
m_grid('box','fancy','linest','none','tickdir','out','fontsize',13);
caxis([cmin cmax])
cRange=caxis;
hold on
qp1 = m_quiver(double(CROCO_lon(1:skp:end,1:skp:end)),double(CROCO_lat(1:skp:end,1:skp:end)),...
    double(ugos_summer(1:skp:end,1:skp:end)),...
    double(vgos_summer(1:skp:end,1:skp:end)),'Color',[0,0,0]);
set(qp1, 'AutoScale','on','AutoScaleFactor',scale_factor);
hold on
m_contour(CROCO_lon,CROCO_lat,croco_top,mydepths,'Color',[0.6824,0.9686,0.6118],'LineWidth',1,'LineStyle','--');
caxis(cRange)
subplot(1,2,2)
m_pcolor(CROCO_lon,CROCO_lat,sst_winter);
shading flat
cmocean('thermal',length(levels))
m_gshhs_f('patch',[.8 .8 .8],'edgecolor','none');
m_grid('box','fancy','linest','none','tickdir','out','fontsize',13);
caxis([cmin cmax])
cRange=caxis;
hold on
qp1 = m_quiver(double(CROCO_lon(1:skp:end,1:skp:end)),double(CROCO_lat(1:skp:end,1:skp:end)),...
    double(ugos_winter(1:skp:end,1:skp:end)),...
    double(vgos_winter(1:skp:end,1:skp:end)),'Color',[0,0,0]);
set(qp1, 'AutoScale','on','AutoScaleFactor',scale_factor);
hold on
m_contour(CROCO_lon,CROCO_lat,croco_top,mydepths,'Color',[0.6824,0.9686,0.6118],'LineWidth',1,'LineStyle','--');
caxis(cRange)

% Create colorbar
hold on
ca = colorbar('eastOutside');
ca.Position = ca.Position + 1e-10;
ca.Label.String = 'SST (°C)';
ca.FontSize = 12;
caxis([cmin cmax]);

% Save the figure

set(gcf, 'InvertHardcopy', 'off')
print('-f1','SST_Currents_sumvswin','-dpng','-r600');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE ANIMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Construct a date array

date = datestr(CROCO_time,'dd-mm-yyyy');

% Pick my dates

start_anim = 2010; 
end_anim = 2010;
fget = 3;     % Days to increment by

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ind_start = min(find(start_anim == Y));
ind_end = max(find(end_anim == Y));
scrsz = get(0,'ScreenSize'); %getting the screensize of the 1 screen
figure('Position',[1 0 1920 1200],'MenuBar','none','ToolBar','none','resize','off')
vidfile = VideoWriter('testmovie.jpeg');
vidfile.FrameRate = 2;
open(vidfile);
for ind = ind_start:fget:ind_end
    clf
    m_proj('miller','long',[double(min(CROCO_lon(:,1))) double(max(CROCO_lon(:,1)))]...
    ,'lat',[double(min(CROCO_lat(1,:))) double(max(CROCO_lat(1,:)))]);
    m_pcolor(CROCO_lon,CROCO_lat,SST(:,:,ind));
    shading flat
    cmocean('thermal',length(levels))
    m_gshhs_f('patch',[.8 .8 .8],'edgecolor','none');
    m_grid('box','fancy','linest','none','tickdir','out','fontsize',13);
    caxis([cmin cmax])
    cRange=caxis;
    hold on
    qp1 = m_quiver(double(CROCO_lon(1:skp:end,1:skp:end)),double(CROCO_lat(1:skp:end,1:skp:end)),...
        double(U(1:skp:end,1:skp:end,ind)),...
        double(V(1:skp:end,1:skp:end,ind)),'Color',[0,0,0]);
    set(qp1, 'AutoScale','on','AutoScaleFactor',scale_factor);
    hold on
    m_contour(CROCO_lon,CROCO_lat,croco_top,mydepths,'Color',[0.6824,0.9686,0.6118],'LineWidth',1,'LineStyle','--');
    caxis(cRange)
    hold on
    title(string(date(ind,:)),'FontSize',15)
    hold on
    ca = colorbar('eastOutside');
    ca.Label.String = 'SST (°C)';
    ca.FontSize = 12;
    caxis([cmin cmax]);
    drawnow
    pause(0.2) 
    F = getframe(gcf); 
    writeVideo(vidfile,F);
end
close(vidfile)


% To put int
normu = U./norm(U);
normv = V./norm(V);


