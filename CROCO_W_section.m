%%%%%%%%%%%%% This script will construct highly detailed horizontal
%%%%%%%%%%%%% transects across a given latitude up to a chosen isobath
%%%%%%%%%%%%% Only do the vertical dimension
%%%%%%%%%%%%% Have summer and winter: so 6 subplots for two each region

%%% make sure that you have croco_tools in your matlab path (I use a
%%% start.m script

%% roms data

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

% Once-off gridding data and dimensions
file = strcat(CROCO_path,'/','croco_avg_Y',string(Ymin),'M',string(1),'.nc');
% Read in the mask
mask=ncread(file,'mask_rho');
mask(mask==0)=nan;
% Define my region, will focus on the coastal domain
CROCO_lat=ncread(strcat(CROCO_path,'/','croco_avg_Y',string(Ymin),'M',string(1),'.nc'),'lat_rho');
CROCO_lon=ncread(strcat(CROCO_path,'/','croco_avg_Y',string(Ymin),'M',string(1),'.nc'),'lon_rho');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESSING  GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define my region, will focus on the coastal domain
croco_top = ncread(strcat(CROCO_path,'/','croco_avg_Y',string(Ymin),'M',string(1),'.nc'),'h');
croco_top = croco_top.*mask;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find the long and lat coordinate pairs for an isobath

mydepth = 1000;
lon_tmp = CROCO_lon(:,1);
lat_tmp = CROCO_lat(1,:);

for ii = 1:size(croco_top,2)
 ind = min(find(croco_top(:,ii)<=mydepth));
 if isempty(ind)
    data(:,ii) = NaN;
 else
    data(:,ii) = lon_tmp(ind);
 end
end

loc_300 = [data' lat_tmp'];

% Pick a latitude to go across

myLat = -34.2924;
LT_index = min(find(lat_tmp>=myLat));
myLon_min= loc_300(find(lat_tmp(LT_index) == loc_300(:,2))); 
myLon_max = lon_tmp(min(find(isnan(mask(:,LT_index)))));

% Plot the transect location

cmin = min(min(croco_top)); % Minimum and maximum temperature bounds
cmax = max(max(croco_top));
levels = [cmin:250:cmax];

figure
m_proj('miller','long',[double(min(min((CROCO_lon)))) double(max(max(CROCO_lon)))]...
    ,'lat',[double(min(min(CROCO_lat))) double(max(max(CROCO_lat)))]);
m_contourf(double(CROCO_lon),double(CROCO_lat),double(croco_top),double(levels),'ShowText','on');
m_gshhs_f('patch',[.7 .7 .7],'edgecolor','none');
m_grid('box','fancy','linest','none','tickdir','out','fontsize',13);
colormap([1 1 1])
hold on
m_line([lon_tmp(min(find(croco_top(:,LT_index)<=mydepth))) lon_tmp(min(find(isnan(mask(:,LT_index)))))]...
    ,[lat_tmp(LT_index) lat_tmp(LT_index)],'linewi',1.5,'color','k');

%% a vertical section

start

lonsec=double([myLon_min myLon_max]);
latsec=double([myLat myLat]);

W = [];
V = [];
U = [];
SALT = [];
TEMP = [];
X = [];
Z =[];
TIME = [];

for i = Ymin:Ymax
    for j = 1:12
    file = strcat(CROCO_path,'/','croco_avg_Y',string(i),'M',string(j),'.nc');
    file = convertStringsToChars(file);
    if exist(file, 'file')
        disp(file)
        disp('Reading data')
        % Get the time domain 
        time = ncread(file,'time');
        % Get temperature
        for k = 1:length(time)
            [x(:,:,k),z(:,:,k),temp(:,:,k)] = get_section(file,file,lonsec,latsec,'temp',k);
            [~,~,v(:,:,k)] = get_section(file,file,lonsec,latsec,'v',k);
            [~,~,w(:,:,k)] = get_section(file,file,lonsec,latsec,'w',k);
            [~,~,u(:,:,k)] = get_section(file,file,lonsec,latsec,'u',k);
            [~,~,salt(:,:,k)] = get_section(file,file,lonsec,latsec,'salt',k);
        end
        %  Store arrays
        disp('Storing for the vertical sections')
        U = cat(3,U,u);
        V = cat(3,V,v);
        W = cat(3,W,w);
        TEMP = cat(3,TEMP,temp);
        SALT = cat(3,SALT,salt);
        X = cat(3,X,x);
        Z = cat(3,Z,z);
        TIME = cat(1,TIME,time);
    else
        disp(strcat('No data for',file))
    end
    clear w u v salt temp time
    end
end

% Clean-up

[~, ia, ~] = unique(TIME);
TEMP = TEMP(:,:,ia);
SALT = SALT(:,:,ia);
U = U(:,:,ia);
V = V(:,:,ia);
W = W(:,:,ia);
X = X(:,:,ia);
Z = Z(:,:,ia);
TIME = TIME(ia);

% Create a structure to back-up run
name = 'SAMBA';
data = struct('X',X,'Z',Z,'U',U,'V',V,'W',W,'TEMP',TEMP,'SALT',SALT,'TIME',TIME);  % Create data structure
save(strcat(name,'',string(-34)),'-struct','data');    % Save structure

%%
% Plot a gradient of the surface SST

%sst_sum = nanmean(TEMP_sum,3);
%dn = sst_sum(end,2:end)-sst_sum(end,1:end-1);
%dis_sum = nanmean(X,3);
%z_sum = nanmean(Z,3);
%dx = dis_sum(end,2:end) - dis_sum(end,1:end-1);
%
%grad_sum = dn./dx;
%
%grad_sum = [grad_sum(1) + (grad_sum(2) - grad_sum(1)), grad_sum];
%plot(fliplr(-dis_sum),grad_sum)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESIGN FIGURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The figure will have six frames: top will be
% The summer: showing Vertical velocity, meridional and zonal data 
% and the bottom will be the winter data. To give greater resolution, the 
% upper 50 m for all plots will be enhanced

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADD the Fronts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numbers appear in the order they are entered:

% St Helena Bay: 111, 112; 46; 62
% Mid Shelf: 104, 118; 31, 41
% Namaqua: 133, 155; 58, 80

%%%%%%%%%%%%%%%%%%%%%%%%

figure

% Vertical velocities
cmin = -2e-5;
cmax = 2e-5;
numlevels = 10;
levels = [cmin:(cmax-cmin)/numlevels:cmax];
scale_factor = 1;
mytemps = [10,12,14,16,18];
myrho = [1024:1:1028];
mydepths = [200,300,500];
surface = 50; % Depth we want to cut the plots at 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ST HELENA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load St_Helena-32.mat

SB_sum = 111;
SB_win = 112;

UP_sum = 46;
UP_win = 62;

% Get the density and MLD data

for i = 1:length(TIME)
    for j = 1:size(TEMP,2)
        z=flipud(squeeze(Z(:,j,i)));
        t=flipud(squeeze(TEMP(:,j,i)));
        [MLD_T(j,i),qe,imf]=get_mld(z,t);
    end
end

for k = 1:length(TIME)
    [rho(:,:,k),bvf(:,:,k)] = rho_eos(squeeze(TEMP(:,:,k)),squeeze(SALT(:,:,k)),squeeze(Z(:,:,k)));
end

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

U_sum = U(:,:,ind_sum);
U_win = U(:,:,ind_win);

V_sum = V(:,:,ind_sum);
V_win = V(:,:,ind_win);

W_sum = W(:,:,ind_sum);
W_win = W(:,:,ind_win);

TEMP_sum = TEMP(:,:,ind_sum);
TEMP_win = TEMP(:,:,ind_win);

RHO_sum = rho(:,:,ind_sum);
RHO_win = rho(:,:,ind_win);

MLD_sum = MLD_T(:,ind_sum);
MLD_win = MLD_T(:,ind_win);

X_sum = X(:,:,ind_sum);
X_win = X(:,:,ind_win);

Z_sum = Z(:,:,ind_sum);
Z_win = Z(:,:, ind_win);

%%%%%%%%%%%%%%%%%%%%%%

subplot(6,3,1)
X_plot = nanmean(X_sum,3);
Z_plot = nanmean(Z_sum,3);
pcolor(fliplr(-X_plot),Z_plot,nanmean(W_sum,3))%,levels,'Showtext','off')
xl = xline(-UP_sum,'-');
xl.LineWidth = 4;
xl.Color = [0.75, 0, 0.75];
xl = xline(-SB_sum,'-');
xl.LineWidth = 4;
xl.Color = [0.0039, 0.5294, 0.0745];
hold on
cmocean('balance',length(levels))
shading interp
caxis([cmin cmax])
cRange=caxis;
hold on
[C_W,h_W] = contour(fliplr(-X_plot),Z_plot,nanmean(W_sum,3),levels,'Color','k','LineWidth',1,'LineStyle','-');
clabel(C_W,h_W,'FontSize',13,'Color','k')
caxis(cRange)
hold on
[C_T,h_T] = contour(fliplr(-X_plot),Z_plot,nanmean(TEMP_sum,3),mytemps,'Color','r','LineWidth',1.5,'LineStyle','--');
clabel(C_T,h_T,'FontSize',13,'Color','red')
caxis(cRange)
hold on
[C_R,h_R] = contour(fliplr(-X_plot),Z_plot,nanmean(RHO_sum,3),myrho,'Color','g','LineWidth',1.5,'LineStyle','--');
clabel(C_R,h_R,'FontSize',13,'Color','green')
caxis(cRange)
hold on
plot(fliplr(-X_plot(1,:)),nanmean(MLD_sum,2),'Color','y','LineWidth',1.5,'LineStyle','-');
hold on 
plot(fliplr(-X_plot(1,:)),Z_plot(1,:),'-k','LineWidth',2);
caxis(cRange)
hold on
ax = gca;
set(gca,'color',[0.8 0.8 0.8])
ax.FontSize = 13;
set(gca,'XTick',[])
title("(a) St Helena Bay: 32^{\circ} S [DJF]",'FontSize',12)
%xlabel('Distance along section [km]')
%ylabel('Depth [m]')
ylim([-surface 0])

clear C_T h_T C_R h_R C_W h_W

h1 = subplot(6,3,[4 7]);
X_plot = nanmean(X_sum,3);
Z_plot = nanmean(Z_sum,3);
pcolor(fliplr(-X_plot),Z_plot,nanmean(W_sum,3))%,levels,'Showtext','off')
cmocean('balance',length(levels))
shading interp
caxis([cmin cmax])
cRange=caxis;
hold on
[C_W,h_W] = contour(fliplr(-X_plot),Z_plot,nanmean(W_sum,3),levels,'Color','k','LineWidth',1,'LineStyle','-');
clabel(C_W,h_W,'FontSize',13,'Color','k')
caxis(cRange)
hold on
[C_T,h_T] = contour(fliplr(-X_plot),Z_plot,nanmean(TEMP_sum,3),mytemps,'Color','r','LineWidth',1.5,'LineStyle','--');
clabel(C_T,h_T,'FontSize',13,'Color','red')
caxis(cRange)
hold on
[C_R,h_R] = contour(fliplr(-X_plot),Z_plot,nanmean(RHO_sum,3),myrho,'Color','g','LineWidth',1.5,'LineStyle','--');
clabel(C_R,h_R,'FontSize',13,'Color','green')
caxis(cRange)
hold on
plot(fliplr(-X_plot(1,:)),nanmean(MLD_sum,2),'Color','y','LineWidth',1.5,'LineStyle','-');
hold on 
plot(fliplr(-X_plot(1,:)),Z_plot(1,:),'-k','LineWidth',2);
caxis(cRange)
hold on
ax = gca;
set(gca,'color',[0.8 0.8 0.8])
ax.FontSize = 13;
%xlabel('Distance along section [km]')
ylab = ylabel('Depth [m]');
ylab.Position(2) = 0; % change vertical position of ylabel
ylab.Position(1) = ylab.Position(1) + 60;
ylim([-1000 -surface])
pos = get(h1,'Position');
pos(2) = pos(2) + 0.03 ;                         % shift up
set( h1, 'Position', pos ) ;

% Add the legend

qw{1} = plot(nan,nan, 'Color',[0.75, 0, 0.75],'LineWidth',2);
qw{2} = plot(nan,nan ,'Color',[0.0039, 0.5294, 0.0745],'LineWidth',2);
l = legend(h1,[qw{:}], {'UF','SBF',}, 'location', 'best');
l.FontSize = 12;

h1 = subplot(6,3,[13 16]);
X_plot = nanmean(X_win,3);
Z_plot = nanmean(Z_win,3);
pcolor(fliplr(-X_plot),Z_plot,nanmean(W_win,3))%,levels,'Showtext','off')
cmocean('balance',length(levels))
shading interp
caxis([cmin cmax])
cRange=caxis;
hold on
[C_W,h_W] = contour(fliplr(-X_plot),Z_plot,nanmean(W_win,3),levels,'Color','k','LineWidth',1,'LineStyle','-');
clabel(C_W,h_W,'FontSize',13,'Color','k')
caxis(cRange)
hold on
[C_T,h_T] = contour(fliplr(-X_plot),Z_plot,nanmean(TEMP_win,3),mytemps,'Color','r','LineWidth',1.5,'LineStyle','--');
clabel(C_T,h_T,'FontSize',13,'Color','red')
caxis(cRange)
hold on
[C_R,h_R] = contour(fliplr(-X_plot),Z_plot,nanmean(RHO_win,3),myrho,'Color','g','LineWidth',1.5,'LineStyle','--');
clabel(C_R,h_R,'FontSize',13,'Color','green')
caxis(cRange)
hold on
plot(fliplr(-X_plot(1,:)),nanmean(MLD_win,2),'Color','y','LineWidth',1.5,'LineStyle','-');
hold on 
plot(fliplr(-X_plot(1,:)),Z_plot(1,:),'-k','LineWidth',2);
caxis(cRange)
hold on
ax = gca;
set(gca,'color',[0.8 0.8 0.8])
ax.FontSize = 13;
xlabel('Distance along section [km]')
ylab = ylabel('Depth [m]');
ylab.Position(2) = 0; % change vertical position of ylabel
ylab.Position(1) = ylab.Position(1) + 60;
ylim([-1000 -surface])
pos = get(h1,'Position');
pos(2) = pos(2) - 0.02 ;                         % shift down
set( h1, 'Position', pos ) ;

h1 = subplot(6,3,10);
X_plot = nanmean(X_win,3);
Z_plot = nanmean(Z_win,3);
pcolor(fliplr(-X_plot),Z_plot,nanmean(W_win,3))%,levels,'Showtext','off')
hold on
xl = xline(-UP_win,'-o');
xl.LineWidth = 4;
xl.Color = [0.75, 0, 0.75];
xl = xline(-SB_win,'-o');
xl.LineWidth = 4;
xl.Color = [0.0039, 0.5294, 0.0745];
cmocean('balance',length(levels))
shading interp
caxis([cmin cmax])
cRange=caxis;
hold on
[C_W,h_W] = contour(fliplr(-X_plot),Z_plot,nanmean(W_win,3),levels,'Color','k','LineWidth',1,'LineStyle','-');
clabel(C_W,h_W,'FontSize',13,'Color','k')
caxis(cRange)
hold on
[C_T,h_T] = contour(fliplr(-X_plot),Z_plot,nanmean(TEMP_win,3),mytemps,'Color','r','LineWidth',1.5,'LineStyle','--');
clabel(C_T,h_T,'FontSize',13,'Color','red')
caxis(cRange)
hold on
[C_R,h_R] = contour(fliplr(-X_plot),Z_plot,nanmean(RHO_win,3),myrho,'Color','g','LineWidth',1.5,'LineStyle','--');
clabel(C_R,h_R,'FontSize',13,'Color','green')
caxis(cRange)
hold on
plot(fliplr(-X_plot(1,:)),nanmean(MLD_win,2),'Color','y','LineWidth',1.5,'LineStyle','-');
hold on 
plot(fliplr(-X_plot(1,:)),Z_plot(1,:),'-k','LineWidth',2);
caxis(cRange)
hold on
ax = gca;
set(gca,'color',[0.8 0.8 0.8])
ax.FontSize = 13;
set(gca,'XTick',[])
title("(d) St Helena Bay: 32^{\circ} S [JJA]",'FontSize',12)
ylim([-surface 0])
pos = get(h1,'Position');
pos(2) = pos(2) - 0.05 ;                         % shift down
set( h1, 'Position', pos );

clear C_T h_T C_R h_R C_W h_W

clear MLD_T qe imf rho bvf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mid shelf region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load Mid_shelf-31.mat

SB_sum = 104;
SB_win = 118;

UP_sum = 31;
UP_win = 41;

% Get the density and MLD data

for i = 1:length(TIME)
    for j = 1:size(TEMP,2)
        z=flipud(squeeze(Z(:,j,i)));
        t=flipud(squeeze(TEMP(:,j,i)));
        [MLD_T(j,i),qe,imf]=get_mld(z,t);
    end
end

for k = 1:length(TIME)
    [rho(:,:,k),bvf(:,:,k)] = rho_eos(squeeze(TEMP(:,:,k)),squeeze(SALT(:,:,k)),squeeze(Z(:,:,k)));
end

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

U_sum = U(:,:,ind_sum);
U_win = U(:,:,ind_win);

V_sum = V(:,:,ind_sum);
V_win = V(:,:,ind_win);

W_sum = W(:,:,ind_sum);
W_win = W(:,:,ind_win);

TEMP_sum = TEMP(:,:,ind_sum);
TEMP_win = TEMP(:,:,ind_win);

RHO_sum = rho(:,:,ind_sum);
RHO_win = rho(:,:,ind_win);

MLD_sum = MLD_T(:,ind_sum);
MLD_win = MLD_T(:,ind_win);

X_sum = X(:,:,ind_sum);
X_win = X(:,:,ind_win);

Z_sum = Z(:,:,ind_sum);
Z_win = Z(:,:, ind_win);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(6,3,2)
X_plot = nanmean(X_sum,3);
Z_plot = nanmean(Z_sum,3);
pcolor(fliplr(-X_plot),Z_plot,nanmean(W_sum,3))%,levels,'Showtext','off')
hold on
xl = xline(-UP_sum,'-o');
xl.LineWidth = 4;
xl.Color = [0.75, 0, 0.75];
xl = xline(-SB_sum,'-o');
xl.LineWidth = 4;
xl.Color = [0.0039, 0.5294, 0.0745];
cmocean('balance',length(levels))
shading interp
caxis([cmin cmax])
cRange=caxis;
hold on
[C_W,h_W] = contour(fliplr(-X_plot),Z_plot,nanmean(W_sum,3),levels,'Color','k','LineWidth',1,'LineStyle','-');
clabel(C_W,h_W,'FontSize',13,'Color','k')
caxis(cRange)
hold on
[C_T,h_T] = contour(fliplr(-X_plot),Z_plot,nanmean(TEMP_sum,3),mytemps,'Color','r','LineWidth',1.5,'LineStyle','--');
clabel(C_T,h_T,'FontSize',13,'Color','red')
caxis(cRange)
hold on
[C_R,h_R] = contour(fliplr(-X_plot),Z_plot,nanmean(RHO_sum,3),myrho,'Color','g','LineWidth',1.5,'LineStyle','--');
clabel(C_R,h_R,'FontSize',13,'Color','green')
caxis(cRange)
hold on
plot(fliplr(-X_plot(1,:)),nanmean(MLD_sum,2),'Color','y','LineWidth',1.5,'LineStyle','-');
hold on 
plot(fliplr(-X_plot(1,:)),Z_plot(1,:),'-k','LineWidth',2);
caxis(cRange)
hold on
ax = gca;
set(gca,'color',[0.8 0.8 0.8])
ax.FontSize = 13;
set(gca,'XTick',[])
title("(b) Mid-shelf: 31^{\circ} S [DJF]",'FontSize',12)
%xlabel('Distance along section [km]')
%ylabel('Depth [m]')
ylim([-surface 0])

clear C_T h_T C_R h_R C_W h_W

h1 = subplot(6,3,[5 8]);
X_plot = nanmean(X_sum,3);
Z_plot = nanmean(Z_sum,3);
pcolor(fliplr(-X_plot),Z_plot,nanmean(W_sum,3))%,levels,'Showtext','off')
cmocean('balance',length(levels))
shading interp
caxis([cmin cmax])
cRange=caxis;
hold on
[C_W,h_W] = contour(fliplr(-X_plot),Z_plot,nanmean(W_sum,3),levels,'Color','k','LineWidth',1,'LineStyle','-');
clabel(C_W,h_W,'FontSize',13,'Color','k')
caxis(cRange)
hold on
[C_T,h_T] = contour(fliplr(-X_plot),Z_plot,nanmean(TEMP_sum,3),mytemps,'Color','r','LineWidth',1.5,'LineStyle','--');
clabel(C_T,h_T,'FontSize',13,'Color','red')
caxis(cRange)
hold on
[C_R,h_R] = contour(fliplr(-X_plot),Z_plot,nanmean(RHO_sum,3),myrho,'Color','g','LineWidth',1.5,'LineStyle','--');
clabel(C_R,h_R,'FontSize',13,'Color','green')
caxis(cRange)
hold on
plot(fliplr(-X_plot(1,:)),nanmean(MLD_sum,2),'Color','y','LineWidth',1.5,'LineStyle','-');
hold on 
plot(fliplr(-X_plot(1,:)),Z_plot(1,:),'-k','LineWidth',2);
caxis(cRange)
hold on
ax = gca;
set(gca,'color',[0.8 0.8 0.8])
ax.FontSize = 13;
%xlabel('Distance along section [km]')
%ylab = ylabel('Depth [m]');
%ylab.Position(2) = 0; % change vertical position of ylabel
%ylab.Position(1) = ylab.Position(1) + 60;
ylim([-1000 -surface])
pos = get(h1,'Position');
pos(2) = pos(2) + 0.03 ;                         % shift up
set( h1, 'Position', pos ) ;

% % Create colorbar
hold on
ca = colorbar('southOutside');
ca.Position = ca.Position + 1e-10;
ca.Position(2) = ca.Position(2) - 0.12;
ca.Label.String = 'Vertical velocities [m/s]';
ca.FontSize = 12;
caxis([cmin cmax]);

h1 = subplot(6,3,[14 17]);
X_plot = nanmean(X_win,3);
Z_plot = nanmean(Z_win,3);
pcolor(fliplr(-X_plot),Z_plot,nanmean(W_win,3))%,levels,'Showtext','off')
cmocean('balance',length(levels))
shading interp
caxis([cmin cmax])
cRange=caxis;
hold on
[C_W,h_W] = contour(fliplr(-X_plot),Z_plot,nanmean(W_win,3),levels,'Color','k','LineWidth',1,'LineStyle','-');
clabel(C_W,h_W,'FontSize',13,'Color','k')
caxis(cRange)
hold on
[C_T,h_T] = contour(fliplr(-X_plot),Z_plot,nanmean(TEMP_win,3),mytemps,'Color','r','LineWidth',1.5,'LineStyle','--');
clabel(C_T,h_T,'FontSize',13,'Color','red')
caxis(cRange)
hold on
[C_R,h_R] = contour(fliplr(-X_plot),Z_plot,nanmean(RHO_win,3),myrho,'Color','g','LineWidth',1.5,'LineStyle','--');
clabel(C_R,h_R,'FontSize',13,'Color','green')
caxis(cRange)
hold on
plot(fliplr(-X_plot(1,:)),nanmean(MLD_win,2),'Color','y','LineWidth',1.5,'LineStyle','-');
hold on 
plot(fliplr(-X_plot(1,:)),Z_plot(1,:),'-k','LineWidth',2);
caxis(cRange)
hold on
ax = gca;
set(gca,'color',[0.8 0.8 0.8])
ax.FontSize = 13;
xlabel('Distance along section [km]')
%ylab = ylabel('Depth [m]');
%ylab.Position(2) = 0; % change vertical position of ylabel
%ylab.Position(1) = ylab.Position(1) + 60;
ylim([-1000 -surface])
pos = get(h1,'Position');
pos(2) = pos(2) - 0.02 ;                         % shift down
set( h1, 'Position', pos ) ;

h1 = subplot(6,3,11);
X_plot = nanmean(X_win,3);
Z_plot = nanmean(Z_win,3);
pcolor(fliplr(-X_plot),Z_plot,nanmean(W_win,3))%,levels,'Showtext','off')
hold on
xl = xline(-UP_win,'-o');
xl.LineWidth = 4;
xl.Color = [0.75, 0, 0.75];
xl = xline(-SB_win,'-o');
xl.LineWidth = 4;
xl.Color = [0.0039, 0.5294, 0.0745];
cmocean('balance',length(levels))
shading interp
caxis([cmin cmax])
cRange=caxis;
hold on
[C_W,h_W] = contour(fliplr(-X_plot),Z_plot,nanmean(W_win,3),levels,'Color','k','LineWidth',1,'LineStyle','-');
clabel(C_W,h_W,'FontSize',13,'Color','k')
caxis(cRange)
hold on
[C_T,h_T] = contour(fliplr(-X_plot),Z_plot,nanmean(TEMP_win,3),mytemps,'Color','r','LineWidth',1.5,'LineStyle','--');
clabel(C_T,h_T,'FontSize',13,'Color','red')
caxis(cRange)
hold on
[C_R,h_R] = contour(fliplr(-X_plot),Z_plot,nanmean(RHO_win,3),myrho,'Color','g','LineWidth',1.5,'LineStyle','--');
clabel(C_R,h_R,'FontSize',13,'Color','green')
caxis(cRange)
hold on
plot(fliplr(-X_plot(1,:)),nanmean(MLD_win,2),'Color','y','LineWidth',1.5,'LineStyle','-');
hold on 
plot(fliplr(-X_plot(1,:)),Z_plot(1,:),'-k','LineWidth',2);
caxis(cRange)
hold on
ax = gca;
set(gca,'color',[0.8 0.8 0.8])
ax.FontSize = 13;
set(gca,'XTick',[])
title("(e) Mid-shelf: 31^{\circ} S [JJA]",'FontSize',12)
ylim([-surface 0])
pos = get(h1,'Position');
pos(2) = pos(2) - 0.05 ;                         % shift down
set( h1, 'Position', pos ) ;

clear C_T h_T C_R h_R C_W h_W

clear MLD_T qe imf rho bvf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Namaqua
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Namaqua: 133, 155; 58, 80
load Namaqua-29.mat

SB_sum = 133;
SB_win = 155;

UP_sum = 58;
UP_win = 80;

% Get the density and MLD data

for i = 1:length(TIME)
    for j = 1:size(TEMP,2)
        z=flipud(squeeze(Z(:,j,i)));
        t=flipud(squeeze(TEMP(:,j,i)));
        [MLD_T(j,i),qe,imf]=get_mld(z,t);
    end
end

for k = 1:length(TIME)
    [rho(:,:,k),bvf(:,:,k)] = rho_eos(squeeze(TEMP(:,:,k)),squeeze(SALT(:,:,k)),squeeze(Z(:,:,k)));
end

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

U_sum = U(:,:,ind_sum);
U_win = U(:,:,ind_win);

V_sum = V(:,:,ind_sum);
V_win = V(:,:,ind_win);

W_sum = W(:,:,ind_sum);
W_win = W(:,:,ind_win);

TEMP_sum = TEMP(:,:,ind_sum);
TEMP_win = TEMP(:,:,ind_win);

RHO_sum = rho(:,:,ind_sum);
RHO_win = rho(:,:,ind_win);

MLD_sum = MLD_T(:,ind_sum);
MLD_win = MLD_T(:,ind_win);

X_sum = X(:,:,ind_sum);
X_win = X(:,:,ind_win);

Z_sum = Z(:,:,ind_sum);
Z_win = Z(:,:, ind_win);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(6,3,3)
X_plot = nanmean(X_sum,3);
Z_plot = nanmean(Z_sum,3);
pcolor(fliplr(-X_plot),Z_plot,nanmean(W_sum,3))%,levels,'Showtext','off')
hold on
xl = xline(-UP_sum,'-o');
xl.LineWidth = 4;
xl.Color = [0.75, 0, 0.75];
xl = xline(-SB_sum,'-o');
xl.LineWidth = 4;
xl.Color = [0.0039, 0.5294, 0.0745];
cmocean('balance',length(levels))
shading interp
caxis([cmin cmax])
cRange=caxis;
hold on
[C_W,h_W] = contour(fliplr(-X_plot),Z_plot,nanmean(W_sum,3),levels,'Color','k','LineWidth',1,'LineStyle','-');
clabel(C_W,h_W,'FontSize',13,'Color','k')
caxis(cRange)
hold on
[C_T,h_T] = contour(fliplr(-X_plot),Z_plot,nanmean(TEMP_sum,3),mytemps,'Color','r','LineWidth',1.5,'LineStyle','--');
clabel(C_T,h_T,'FontSize',13,'Color','red')
caxis(cRange)
hold on
[C_R,h_R] = contour(fliplr(-X_plot),Z_plot,nanmean(RHO_sum,3),myrho,'Color','g','LineWidth',1.5,'LineStyle','--');
clabel(C_R,h_R,'FontSize',13,'Color','green')
caxis(cRange)
hold on
plot(fliplr(-X_plot(1,:)),nanmean(MLD_sum,2),'Color','y','LineWidth',1.5,'LineStyle','-');
hold on 
plot(fliplr(-X_plot(1,:)),Z_plot(1,:),'-k','LineWidth',2);
caxis(cRange)
hold on
ax = gca;
set(gca,'color',[0.8 0.8 0.8])
ax.FontSize = 13;
set(gca,'XTick',[])
title("(c) Namaqua: 29^{\circ} [DJF]",'FontSize',12)
%xlabel('Distance along section [km]')
%ylabel('Depth [m]')
ylim([-surface 0])

clear C_T h_T C_R h_R C_W h_W

h1 = subplot(6,3,[6 9]);
X_plot = nanmean(X_sum,3);
Z_plot = nanmean(Z_sum,3);
pcolor(fliplr(-X_plot),Z_plot,nanmean(W_sum,3))%,levels,'Showtext','off')
cmocean('balance',length(levels))
shading interp
caxis([cmin cmax])
cRange=caxis;
hold on
[C_W,h_W] = contour(fliplr(-X_plot),Z_plot,nanmean(W_sum,3),levels,'Color','k','LineWidth',1,'LineStyle','-');
clabel(C_W,h_W,'FontSize',13,'Color','k')
caxis(cRange)
hold on
[C_T,h_T] = contour(fliplr(-X_plot),Z_plot,nanmean(TEMP_sum,3),mytemps,'Color','r','LineWidth',1.5,'LineStyle','--');
clabel(C_T,h_T,'FontSize',13,'Color','red')
caxis(cRange)
hold on
[C_R,h_R] = contour(fliplr(-X_plot),Z_plot,nanmean(RHO_sum,3),myrho,'Color','g','LineWidth',1.5,'LineStyle','--');
clabel(C_R,h_R,'FontSize',13,'Color','green')
caxis(cRange)
hold on
plot(fliplr(-X_plot(1,:)),nanmean(MLD_sum,2),'Color','y','LineWidth',1.5,'LineStyle','-');
hold on 
plot(fliplr(-X_plot(1,:)),Z_plot(1,:),'-k','LineWidth',2);
caxis(cRange)
hold on
ax = gca;
set(gca,'color',[0.8 0.8 0.8])
ax.FontSize = 13;
%xlabel('Distance along section [km]')
%ylab = ylabel('Depth [m]');
%ylab.Position(2) = 0; % change vertical position of ylabel
%ylab.Position(1) = ylab.Position(1) + 60;
ylim([-1000 -surface])
pos = get(h1,'Position');
pos(2) = pos(2) + 0.03 ;                         % shift up
set( h1, 'Position', pos ) ;

h1 = subplot(6,3,[15 18]);
X_plot = nanmean(X_win,3);
Z_plot = nanmean(Z_win,3);
pcolor(fliplr(-X_plot),Z_plot,nanmean(W_win,3))%,levels,'Showtext','off')
cmocean('balance',length(levels))
shading interp
caxis([cmin cmax])
cRange=caxis;
hold on
[C_W,h_W] = contour(fliplr(-X_plot),Z_plot,nanmean(W_win,3),levels,'Color','k','LineWidth',1,'LineStyle','-');
clabel(C_W,h_W,'FontSize',13,'Color','k')
caxis(cRange)
hold on
[C_T,h_T] = contour(fliplr(-X_plot),Z_plot,nanmean(TEMP_win,3),mytemps,'Color','r','LineWidth',1.5,'LineStyle','--');
clabel(C_T,h_T,'FontSize',13,'Color','red')
caxis(cRange)
hold on
[C_R,h_R] = contour(fliplr(-X_plot),Z_plot,nanmean(RHO_win,3),myrho,'Color','g','LineWidth',1.5,'LineStyle','--');
clabel(C_R,h_R,'FontSize',13,'Color','green')
caxis(cRange)
hold on
plot(fliplr(-X_plot(1,:)),nanmean(MLD_win,2),'Color','y','LineWidth',1.5,'LineStyle','-');
hold on 
plot(fliplr(-X_plot(1,:)),Z_plot(1,:),'-k','LineWidth',2);
caxis(cRange)
hold on
ax = gca;
set(gca,'color',[0.8 0.8 0.8])
ax.FontSize = 13;
xlabel('Distance along section [km]')
%ylab = ylabel('Depth [m]');
%ylab.Position(2) = 0; % change vertical position of ylabel
%ylab.Position(1) = ylab.Position(1) + 60;
ylim([-1000 -surface])
pos = get(h1,'Position');
pos(2) = pos(2) - 0.02 ;                         % shift down
set( h1, 'Position', pos ) ;

h1 = subplot(6,3,12);
X_plot = nanmean(X_win,3);
Z_plot = nanmean(Z_win,3);
pcolor(fliplr(-X_plot),Z_plot,nanmean(W_win,3))%,levels,'Showtext','off')
hold on
xl = xline(-UP_win,'-o');
xl.LineWidth = 4;
xl.Color = [0.75, 0, 0.75];
xl = xline(-SB_win,'-o');
xl.LineWidth = 4;
xl.Color = [0.0039, 0.5294, 0.0745];
cmocean('balance',length(levels))
shading interp
caxis([cmin cmax])
cRange=caxis;
hold on
[C_W,h_W] = contour(fliplr(-X_plot),Z_plot,nanmean(W_win,3),levels,'Color','k','LineWidth',1,'LineStyle','-');
clabel(C_W,h_W,'FontSize',13,'Color','k')
caxis(cRange)
hold on
[C_T,h_T] = contour(fliplr(-X_plot),Z_plot,nanmean(TEMP_win,3),mytemps,'Color','r','LineWidth',1.5,'LineStyle','--');
clabel(C_T,h_T,'FontSize',13,'Color','red')
caxis(cRange)
hold on
[C_R,h_R] = contour(fliplr(-X_plot),Z_plot,nanmean(RHO_win,3),myrho,'Color','g','LineWidth',1.5,'LineStyle','--');
clabel(C_R,h_R,'FontSize',13,'Color','green')
caxis(cRange)
hold on
plot(fliplr(-X_plot(1,:)),nanmean(MLD_win,2),'Color','y','LineWidth',1.5,'LineStyle','-');
hold on 
plot(fliplr(-X_plot(1,:)),Z_plot(1,:),'-k','LineWidth',2);
caxis(cRange)
hold on
ax = gca;
set(gca,'color',[0.8 0.8 0.8])
ax.FontSize = 13;
set(gca,'XTick',[])
title("(f) Namaqua: 29^{\circ} S [JJA]",'FontSize',12)
ylim([-surface 0])
pos = get(h1,'Position');
pos(2) = pos(2) - 0.05 ;                         % shift down
set( h1, 'Position', pos ) ;

clear C_T h_T C_R h_R C_W h_W
clear MLD_T qe imf rho bvf
%% Save the figure

set(gcf, 'InvertHardcopy', 'off')
print('-f1','SBUS_cross_W','-dpng','-r600');