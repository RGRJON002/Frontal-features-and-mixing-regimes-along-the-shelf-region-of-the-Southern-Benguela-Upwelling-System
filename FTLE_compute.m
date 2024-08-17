%%%% This script is intended to store and compute FTLE
%%%% The objective is to compute FTLE on a higher resolution
%%%% grid and then interpolate the high res grid back to 
%%%% the resolution of CROCO for plotting and analysis
%%%% Written: Jonathan Rogerson
%%%% University of Cape Town
%%%% July 2024

addpath '/home/jono/Documents/MATLAB/LCS-Tool-master'
addpath '/home/jono/Documents/PhD/CHAPTER_2'
addpath /home/jono/Documents/MATLAB/CMOCEAN
addpath /home/jono/CROCO/croco_tools/UTILITIES/m_map1.4h
addpath /media/data/DATASETS_CROCOTOOLS/m_map1.4f
addpath '/usr/local/MATLAB/R2020a/toolbox/matlab/imagesci/'
addpath '/media/data/DIAGNOSTICS/FRONTS'
addpath '/home/jono/CROCO/croco_tools/Preprocessing_tools'

% Declare some croco variables for the grid process

%%%%%%%%%%%%%%%%%%%%%%%%%%
% CROCO
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Same temporal domain overr which the FTLE was computed
CROCO_path = '/media/jono/SBUS/SBUS_3km/CHPC_OUTPUT';
Ymin = 2014;
Ymax = 2015;
Yorig = 1990;  % This will change for the origin for FTLE but is fine for now

%%%%%%%%%
% CROCO
%%%%%%%%%

% Once-off gridding data and dimensions
file = strcat(CROCO_path,'/','croco_avg_Y',string(Ymin),'M',string(1),'.nc');
% Define my region, will focus on the coastal domain
CROCO_lat=ncread(strcat(CROCO_path,'/','croco_avg_Y',string(Ymin),'M',string(1),'.nc'),'lat_rho');
CROCO_lon=ncread(strcat(CROCO_path,'/','croco_avg_Y',string(Ymin),'M',string(1),'.nc'),'lon_rho');

% Domain of FTLE computation

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

% Create domain equal in resolution to size(FTLE) noting the zoom

domain = [15,20;-33.5,-28];
resolutionX = size(CROCO_lon,1)*5;
[resolutionY,deltaX] = equal_resolution(domain,resolutionX);
resolution = [resolutionX,resolutionY];

% Increase resolution

lat_tmp = min(CROCO_lat(1,:)):(max(CROCO_lat(1,:))-min(CROCO_lat(1,:)))/resolutionY:max(CROCO_lat(1,:));
lon_tmp = min(CROCO_lon(:,1)):(max(CROCO_lon(:,1))-min(CROCO_lon(:,1)))/resolutionX:max(CROCO_lon(:,1));

% Clean

lat_tmp(end) = [];
lon_tmp(end) = [];

[XX,YY] = meshgrid(lon_tmp',lat_tmp'); % This will now serve as the template to interpolate from

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FTLE Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load 'SBUS_2014_2015.mat'   % Do need the current data

Yorig = 2014;   % Date of initial save

u(isnan(u)) = 0;
v(isnan(v)) = 0; 

% Number of times to interate base on an integration of 30 days and a
% window of 1 day

num_times = round((time(end)-30)/1);   % Total number. We only want to run for one year
num_times = num_times/2;

shift = 350;
for i = 351:365
    timespan = [1+shift,30+shift];
    disp(timespan)
    domain = [15,20;-33.5,-28];
    resolutionX = size(u,3)*5;
    [resolutionY,deltaX] = equal_resolution(domain,resolutionX);
    resolution = [resolutionX,resolutionY];

    interpMethod = 'nearest';
    vLonInterpolant = griddedInterpolant({time,lat,lon},u,interpMethod);
    vLatInterpolant = griddedInterpolant({time,lat,lon},v,interpMethod);
    lDerivative = @(t,x,~)derivative(t,x,vLonInterpolant,vLatInterpolant);
    incompressible = true;

    % Cauchy-Green strain
    cgEigenvalueFromMainGrid = false;
    cgAuxGridRelDelta = .1;
    [cgEigenvector,cgEigenvalue] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible,'eigenvalueFromMainGrid',cgEigenvalueFromMainGrid,'auxGridRelDelta',cgAuxGridRelDelta);

    % Compute FTLE
    cgEigenvalue2 = reshape(cgEigenvalue(:,2),fliplr(resolution));
    ftle_ = ftle(cgEigenvalue2,diff(timespan));
    
    % Normalize by max FTLE of field
    %ftle_ = ftle_./max(max(ftle_));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Interpolate FTLE to CROCO grid
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Interpolate FTLE to native grid')
    ftle_ = interp2(XX,YY,ftle_,CROCO_lon',CROCO_lat');

    % store
    disp('Storing FTLE')
    myFTLE(:,:,i) = ftle_;
    
    % Store a time-stamp pf the data
    
    CROCO_time = datetime(Yorig,1,1) + days(timespan);
    date = datenum(CROCO_time);
    %date = mean(date);
    date = datetime(date, 'ConvertFrom', 'datenum', 'Format', 'dd-MM-yyyy');
    
    % store
    disp('Storing time')
    mydate(:,i) = date;
    
    shift = shift+1;
end

save('FTLE_zoom_nonorm_2014.mat','lon','lat','mydate','myFTLE');