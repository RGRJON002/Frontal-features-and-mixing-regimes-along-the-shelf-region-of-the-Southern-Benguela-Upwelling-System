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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W = [];
TIME = [];
for i = Ymin:Ymax
    for j = 1:12
    file = strcat(CROCO_path,'/','croco_avg_Y',string(i),'M',string(j),'.nc');
    if exist(file, 'file')
        disp(file)
        disp('Reading data')
        % Read in the time array
        time = ncread(file,'time');
        % Get the vertical grid for that month
        grd = getmydepth(convertStringsToChars(file),10);
        disp('Processing grid')
        % Read in the vertical velocity data
        addpath('/usr/local/MATLAB/R2020a/toolbox/matlab/imagesci/')
        w = ncread(file,'w');
        w = bsxfun(@rdivide, squeeze(nansum(w.*grd,3)), h);
        %  Store arrays
        disp('Storing w and time')
        TIME = cat(1,TIME,time);
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

%%%% First we plot the season

cmin = -1e-6;
cmax = 1e-6;
mydepths = [200,300,500];

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
m_contour(CROCO_lon,CROCO_lat,h,mydepths,'Color','g','LineWidth',1,'LineStyle','--');
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
m_contour(CROCO_lon,CROCO_lat,h,mydepths,'Color','g','LineWidth',1,'LineStyle','--');
caxis(cRange)

% Create colorbar
hold on
ca = colorbar('eastOutside');
ca.Position = ca.Position + 1e-10;
ca.Label.String ='Vertical velocities (m/s)';
ca.FontSize = 12;
caxis([cmin cmax]);

% Save the figure

set(gcf, 'InvertHardcopy', 'off')
print('-f1','W_200','-dpng','-r600');

var=get_hslice(fname,fname,'w',1,-10,'r');

var = vinterp(var_sigma,z,level);  % to code in getdepths

cmin= -1e-6
cmax = 1e-6
pcolor(var)
colorbar
caxis([cmin cmax])
