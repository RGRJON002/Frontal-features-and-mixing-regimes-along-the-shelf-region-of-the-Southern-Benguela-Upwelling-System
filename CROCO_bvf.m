%%%% This script is intended to calculate the Brunt-Vaisala frequency for
%%%% the near shore domain. Focus will be given to the upper 50 m extending
%%%% 20 km (7 grid cells) offshore. 

%%%%%%%%%%%%%%%%%%
% PROCESS CROCO
%%%%%%%%%%%%%%%%%%

addpath /home/jono/CROCO/croco_tools/UTILITIES/m_map1.4h
addpath /media/data/DATASETS_CROCOTOOLS/m_map1.4f
addpath /home/jono/Documents/MATLAB/CMOCEAN
addpath '/home/jono/Documents/MATLAB/CDT-master/cdt'
addpath '/home/jono/Documents/MATLAB/CDT-master/cdt/cdt_data'
addpath '/home/jono/Documents/MATLAB/tight_subplot'
addpath /home/jono/CROCO/croco_tools/Preprocessing_tools
addpath /home/jono/CROCO/croco_tools/Visualization_tools

CROCO_path = '/media/jono/SBUS/SBUS_3km/CHPC_OUTPUT';
Ymin = 2004;
Ymax = 2018;
Yorig = 1990;
make_plot = 1;   % Check the domain

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE coastal region to focus on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define my region, will focus on the coastal domain

CROCO_lat=ncread(strcat(CROCO_path,'/','croco_avg_Y',string(Ymin),'M',string(1),'.nc'),'lat_rho');
CROCO_lon=ncread(strcat(CROCO_path,'/','croco_avg_Y',string(Ymin),'M',string(1),'.nc'),'lon_rho');
mask=ncread(strcat(CROCO_path,'/','croco_avg_Y',string(Ymin),'M',string(1),'.nc'),'mask_rho');
mask(mask==0)=nan;  
h = ncread(strcat(CROCO_path,'/','croco_avg_Y',string(Ymin),'M',string(1),'.nc'),'h');
h = h.*mask;

% Co-ordinates

lat_min = -36;
lat_max = -28;
lon_min = 15;
lon_max = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,idx_lat_min]=min(abs(CROCO_lat(1,:)-lat_min));
[~,idx_lat_max]=min(abs(CROCO_lat(1,:)-lat_max));
[~,idx_lon_min]=min(abs(CROCO_lon(:,1)-lon_min));
[~,idx_lon_max]=min(abs(CROCO_lon(:,1)-lon_max));
% 
CROCO_lat = CROCO_lat(idx_lon_min:idx_lon_max,idx_lat_min:idx_lat_max);
CROCO_lon = CROCO_lon(idx_lon_min:idx_lon_max,idx_lat_min:idx_lat_max);
mask = mask(idx_lon_min:idx_lon_max,idx_lat_min:idx_lat_max);
h = double(h(idx_lon_min:idx_lon_max,idx_lat_min:idx_lat_max));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUICK PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if make_plot == 1 
    figure;
    m_proj('miller','long',[double(min(CROCO_lon(:,1))) double(max(CROCO_lon(:,1)))]...
        ,'lat',[double(min(CROCO_lat(1,:))) double(max(CROCO_lat(1,:)))]);
    m_gshhs_h('patch',[.8 .8 .8],'edgecolor','none');
    m_grid('box','fancy','linest','none','tickdir','out','fontsize',13);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRID INFO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 60;
theta_s    =  5.;
theta_b    =  0.;
hc         = 10.;
vtransform =  2.;
type = 'r';

% Vertical levels to save for the BVF calculation 
g = 9.81;
rho0 = 1025;
zlevels = [30];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To save time for the computation, space will be prealocated but need to
% calculate the length of the time-domain

TIME = [];
for i = Ymin:Ymax
    for j = 1:12
    file = strcat(CROCO_path,'/','croco_avg_Y',string(i),'M',string(j),'.nc');
    if exist(file, 'file')
        disp(file)
        disp('Reading data')
        time = ncread(file,'time');
        % Store arrays
        disp('Storing time')
        TIME = cat(1,TIME,time);
    else
        disp(strcat('No data for',file))
    end
    end
end

[M,L] = size(CROCO_lat);

BVF = [];
for i = Ymin:Ymax
    for j = 1:12
    file = strcat(CROCO_path,'/','croco_avg_Y',string(i),'M',string(j),'.nc');
    if exist(file, 'file')
        disp(file)
        disp('Reading data')    
        % Read in the neccessary info
        zeta = ncread(file,'zeta',[idx_lon_min idx_lat_min 1],[M L inf]);
        zeta = zeta.*mask;
        for k = 1:size(zeta,3) 
            zeta_temp = squeeze(zeta(:,:,k));
            temp = squeeze(ncread(file,'temp',[idx_lon_min idx_lat_min 1 k],[M L inf 1]));
            salt = squeeze(ncread(file,'salt',[idx_lon_min idx_lat_min 1 k],[M L inf 1]));
            temp = permute(temp,[3 1 2]);
            salt = permute(salt,[3 1 2]);
            z_rho = zlevs(h,zeta_temp,theta_s,theta_b,hc,N,type,vtransform);
            z_w = zlevs(h,zeta_temp,theta_s,theta_b,hc,N,'w',vtransform);
            bvf=bvf_eos(temp,salt,z_rho,z_w,g,rho0);
            for l = 1:length(zlevels)
                bvf_lev(:,:,l) = vinterp(bvf,z_w,-zlevels(l));
            end
             BVF = cat(4,BVF,bvf_lev);
        end
%         Store arrays
        disp('Storing BVF complete')
    else
        disp(strcat('No data for',file))
    end
    end
end

%% Get some quick plots

BVF = squeeze(BVF);
%TIME = TIME(1:size(BVF,4));
[~, ia, ~] = unique(TIME);
BVF = BVF(:,:,ia);
TIME = TIME(ia);

Chunk = length(TIME)/4;

BVF_1 = BVF(:,:,1:Chunk);
BVF_2 = BVF(:,:,Chunk+1:Chunk*2);
BVF_3 = BVF(:,:,((Chunk*2)+1):Chunk*3);
BVF_4 = BVF(:,:,((Chunk*3)+1):end);

for ii = 1:4
    save(strcat('BVF_',string(ii),'.mat'),strcat('BVF_',string(ii)),'TIME',...
        'CROCO_lat','CROCO_lon')
end

%% Save data

CROCO_time = datetime(Yorig,1,1) + seconds(TIME);
[Y,MO,D] = datevec(CROCO_time);  

figure
for i = 1:12
    subplot(3,4,i)
    m_proj('miller','long',[double(min(CROCO_lon(:,1))) double(max(CROCO_lon(:,1)))]...
        ,'lat',[double(min(CROCO_lat(1,:))) double(max(CROCO_lat(1,:)))]);
    m_pcolor(CROCO_lon,CROCO_lat,squeeze(BVF(:,:,4,i*30)))
    m_gshhs_h('patch',[.8 .8 .8],'edgecolor','none');
    m_grid('box','fancy','linest','none','tickdir','out','fontsize',13);
    shading interp
    cmocean('dense')
    title(string(CROCO_time(i*30)))
end
hold on
ca = colorbar('eastOutside');
ca.Position = ca.Position + 1e-10;
ca.Label.String = 'N^2 (1/s^2)';
ca.FontSize = 12;

set(gcf, 'InvertHardcopy', 'off')
print('-f1','BVF_30m','-dpng','-r600');
%% Load in a file

% We need to load in the variable ZETA

file = 'croco_avg_Y2005M6.nc';
zeta = ncread(file,'zeta');
zeta = zeta(:,:,1);
%zeta = permute(zeta,[3 1 2]);
h = ncread(file,'h');
temp = ncread(file,'temp');
temp = temp(:,:,:,1);
salt = ncread(file,'salt');
salt = salt(:,:,:,1);

temp = permute(temp,[3 1 2]);
salt = permute(salt,[3 1 2]);

g = 9.81;
rho0 = 1025;
z_rho = zlevs(h,zeta,theta_s,theta_b,hc,N,type,vtransform);
z_w = zlevs(h,zeta,theta_s,theta_b,hc,N,'w',vtransform);

bvf=bvf_eos(temp,salt,z_rho,z_w,g,rho0);

var = vinterp(bvf,z_w,-30);

pcolor(var')
colorbar

figure
plot(squeeze(bvf(:,337,179)),squeeze(z_w(:,337,179)))