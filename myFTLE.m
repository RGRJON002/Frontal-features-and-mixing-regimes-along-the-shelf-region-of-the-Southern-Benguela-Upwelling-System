%%%% This script plots the FTLE on top of the hyperbolic attractors for
%%%% CROCO model output

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
addpath '/home/jono/CROCO/croco_tools/Preprocessing_tools'

CROCO_path = '/media/jono/SBUS/SBUS_3km/CHPC_OUTPUT';
Ymin = 2014;
Ymax = 2015;
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
%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORTANT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We need to have our data in the correct format for the integration
% to work

% Long and lat as Nx1 and Mx1 vectors in units of degrees
% Time: Any unit of time is allowed but days will be used
% Velocity: must be arranged such that lon is first and lat is second (x,y)
% and is in units of deg/day. Reason is that now the spatial and temporal
% coordinates have the same units as the velocity components.

% Convert m/s to m/day

raduis_earth = 6371*10^3;  % meters of earth raduis
day_seconds = 86400;     % Seconds in day
lat_dis = 111*10^3;      % Distance between latitudes

% At a specific latitude, the longitude length will differ, influencing the
% dx component. dy is constant

dlong = (2*pi*raduis_earth*cosd(CROCO_lat))/360; 
u = double((u.*day_seconds./dlong));
v = double((v.*day_seconds./lat_dis));

time_CROCO = datetime(Yorig,1,1) + seconds(TIME);

% Create a days since initialisation vector

res_time = 1;
days = 1:res_time:length(TIME)*res_time;
time= days';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subset the domain of interest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
u = u(idx_lon_min:idx_lon_max,idx_lat_min:idx_lat_max,:);
v = v(idx_lon_min:idx_lon_max,idx_lat_min:idx_lat_max,:);

% Get Lon and lat as Nx1 vectors

lon = double(CROCO_lon(:,1));
lat = double(CROCO_lat(1,:)');

u = permute(u, [3 2 1]);
v = permute(v,[3 2 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save as a data structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = struct('lon',lon,'lat',lat,'time',time,'u',u,'v',v);  % Create data structure
save('SBUS_2014_2015','-struct','data'); 