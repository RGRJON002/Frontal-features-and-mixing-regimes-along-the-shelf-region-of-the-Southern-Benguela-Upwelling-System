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

% Once-off gridding data and dimensions
file = strcat(CROCO_path,'/','croco_avg_Y',string(Ymin),'M',string(1),'.nc');
% Read in the mask
mask=ncread(file,'mask_rho');
mask(mask==0)=nan;
% Define my region, will focus on the coastal domain
CROCO_lat=ncread(strcat(CROCO_path,'/','croco_avg_Y',string(Ymin),'M',string(1),'.nc'),'lat_rho');
CROCO_lon=ncread(strcat(CROCO_path,'/','croco_avg_Y',string(Ymin),'M',string(1),'.nc'),'lon_rho');
% Only need to index lon as lat bounds and eastern extent stay the same

lon_min = 14;
[~,idx_lon_min]=min(abs(CROCO_lon(:,1)-lon_min));
mask = mask(idx_lon_min:end,:);

%%%%%%%%%
% CROCO
%%%%%%%%%
g = 9.81;
U = [];
V = [];
TIME = [];
for i = Ymin:Ymax
    for j = 1:12
    file = strcat(CROCO_path,'/','croco_avg_Y',string(i),'M',string(j),'.nc');
    if exist(file, 'file')
        disp(file)
        disp('Reading data')
        % extract coriolis
        f = ncread(file,'f',[idx_lon_min 1],[inf inf]);
        % extract dx and dy
        dx = ncread(file,'pm',[idx_lon_min 1],[inf inf]);
        dy = ncread(file,'pn',[idx_lon_min 1 ],[inf inf]);
        ssh = ncread(file,'zeta',[idx_lon_min 1 1],[inf inf inf]);
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
        % Read in the time array
        time = ncread(file,'time');
        % Store arrays
        disp('Storing U, V and time')
        TIME = cat(1,TIME,time);
        U = cat(3,U,ug);
        V = cat(3,V,vg);
    else
        disp(strcat('No data for',file))
    end
         % Clean-up 
        %clear ssh dx dy ug vg ff dn dxx dyy time
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POST-PROCESSING  GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define my region, will focus on the coastal domain
croco_top = ncread(strcat(CROCO_path,'/','croco_avg_Y',string(Ymin),'M',string(1),'.nc'),'h',[idx_lon_min 1 ],[inf inf]);
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

% Compute long-term means from U and V 

Ubar = nanmean(U,3);
Vbar = nanmean(V,3);

% equation: EKE = [(u'.^2 + v'.^2)]/2
% Where u' = u - ubar

% Calculate seasonal EKE

% For every u and v point compute EKE
for xx = 1:length(TIME)
    tmpu = U(:,:,xx) - Ubar;
    tmpv = V(:,:,xx) - Vbar;
    EKE(:,:,xx) = 0.5 .* (tmpu.^2 + tmpv.^2);
end

clear tmpu tmpv

% Convert EKE from m^2 s^-2 to cm^2 s^-2

EKE = EKE.*100^2;

% Calculate the curl of the geostrophic currents

for yy = 1:length(TIME)
    tmpu = U(:,:,yy);
    tmpv = V(:,:,yy);
    cav(:,:,yy) = curl(CROCO_lon',CROCO_lat', tmpu', tmpv');
end

clear tmpu tmpv

cav = permute(cav,[2,1,3]);

% Focus on the coastal domain

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

% Index EKE, U, V and curl

EKE_summer = double(nanmean(EKE(:,:,ind_sum),3)).*mask;
EKE_winter = double(nanmean(EKE(:,:,ind_win),3)).*mask;

U_sum = double(nanmean(U(:,:,ind_sum),3)).*mask;
V_sum = double(nanmean(V(:,:,ind_sum),3)).*mask;

U_win = double(nanmean(U(:,:,ind_win),3)).*mask;
V_win = double(nanmean(V(:,:,ind_win),3)).*mask;

curl_sum = double(nanmean(cav(:,:,ind_sum),3)).*mask;
curl_win = double(nanmean(cav(:,:,ind_win),3)).*mask;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT THE FIGURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% First we plot the season for EKE

cmin = 0;
cmax = ceil(round(max(max(nanmean(EKE,3))))/100)*100;
levels = [cmin:200:cmax];
mydepths = [200,300,500];

figure

m_proj('miller','long',[double(CROCO_lon(idx_lon_min,1)) double(CROCO_lon(idx_lon_max,1))]...
    ,'lat',[double(CROCO_lat(1,idx_lat_min)) double(CROCO_lat(1,idx_lat_max))]);
subplot(1,2,1)
m_contourf(CROCO_lon,CROCO_lat,EKE_summer,levels,'ShowText','on');
shading flat
cmocean('curl')
m_gshhs_h('patch',[.8 .8 .8],'edgecolor','none');
m_grid('box','fancy','linest','none','tickdir','out','fontsize',13);
caxis([cmin cmax])
cRange=caxis;
hold on
m_contour(CROCO_lon,CROCO_lat,croco_top,mydepths,'Color','k','LineWidth',1,'LineStyle','--');
caxis(cRange)
subplot(1,2,2)
m_contourf(CROCO_lon,CROCO_lat,EKE_winter,levels,'ShowText','on');
shading flat
cmocean('curl')
m_gshhs_h('patch',[.8 .8 .8],'edgecolor','none');
m_grid('box','fancy','linest','none','tickdir','out','fontsize',13);
caxis([cmin cmax])
cRange=caxis;
hold on
m_contour(CROCO_lon,CROCO_lat,croco_top,mydepths,'Color','k','LineWidth',1,'LineStyle','--');
caxis(cRange)

% Create colorbar
hold on
ca = colorbar('eastOutside');
ca.Position = ca.Position + 1e-10;
ca.Label.String = 'EKE (cm/m^2)';
ca.FontSize = 12;
caxis([cmin cmax]);

% Save the figure

set(gcf, 'InvertHardcopy', 'off')
print('-f1','EKE','-dpng','-r600');

%%%% Plot the current curl with the currents

cmin = -1;
cmax = 1;
scale_factor = 1.5;
mydepths = [200,300,500];
skp = 10;

figure

m_proj('miller','long',[double(CROCO_lon(idx_lon_min,1)) double(CROCO_lon(idx_lon_max,1))]...
    ,'lat',[double(CROCO_lat(1,idx_lat_min)) double(CROCO_lat(1,idx_lat_max))]);
subplot(1,2,1)
m_pcolor(CROCO_lon,CROCO_lat,curl_sum);
shading flat
cmocean('balance')
m_gshhs_h('patch',[.8 .8 .8],'edgecolor','none');
m_grid('box','fancy','linest','none','tickdir','out','fontsize',13);
caxis([cmin cmax])
cRange=caxis;
hold on
qp1 = m_quiver(double(CROCO_lon(1:skp:end,1:skp:end)),double(CROCO_lat(1:skp:end,1:skp:end)),...
    double(U_sum(1:skp:end,1:skp:end)),...
    double(V_sum(1:skp:end,1:skp:end)),'Color',[0,0,0]);
set(qp1, 'AutoScale','on','AutoScaleFactor',scale_factor);
hold on
m_contour(CROCO_lon,CROCO_lat,croco_top,mydepths,'Color','k','LineWidth',1,'LineStyle','--');
caxis(cRange)
subplot(1,2,2)
m_pcolor(CROCO_lon,CROCO_lat,curl_win);
shading flat
cmocean('balance')
m_gshhs_f('patch',[.8 .8 .8],'edgecolor','none');
m_grid('box','fancy','linest','none','tickdir','out','fontsize',13);
caxis([cmin cmax])
cRange=caxis;
hold on
qp1 = m_quiver(double(CROCO_lon(1:skp:end,1:skp:end)),double(CROCO_lat(1:skp:end,1:skp:end)),...
    double(U_win(1:skp:end,1:skp:end)),...
    double(V_win(1:skp:end,1:skp:end)),'Color',[0,0,0]);
set(qp1, 'AutoScale','on','AutoScaleFactor',scale_factor);
hold on
m_contour(CROCO_lon,CROCO_lat,croco_top,mydepths,'Color','k','LineWidth',1,'LineStyle','--');
caxis(cRange)

% Create colorbar
hold on
ca = colorbar('eastOutside');
ca.Position = ca.Position + 1e-10;
ca.Label.String = 'Curl (rad/s)';
ca.FontSize = 12;
caxis([cmin cmax]);

set(gcf, 'InvertHardcopy', 'off')
print('-f2','Current_curl','-dpng','-r600');


    
