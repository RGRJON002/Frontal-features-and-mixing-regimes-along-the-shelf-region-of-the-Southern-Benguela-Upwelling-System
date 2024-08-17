%%%% This script is intended to be used to generate probability maps
%%%% of lagrangian floats released along a region

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /home/jono/CROCO/croco_tools/UTILITIES/m_map1.4h
addpath /media/data/DATASETS_CROCOTOOLS/m_map1.4f
addpath /home/jono/Documents/MATLAB/CMOCEAN
addpath /home/jono/CROCO/Roff
addpath /home/jono/Documents/MATLAB/PIFF-master
addpath /home/jono/Documents/PhD/CHAPTER_2/SCRIPTS

% Inital file
CROCO_path = '/media/data/CHPC_SBUS_3km'; 
CROCO_file = 'croco_N01.nc';

% Read in the grid dimension of one of the CROCO_files

croco_lon = ncread(strcat(CROCO_path,'/',CROCO_file),'lon_rho');
croco_lat = ncread(strcat(CROCO_path,'/',CROCO_file),'lat_rho');
 
% Read in the topography data

croco_top = ncread(strcat(CROCO_path,'/',CROCO_file),'h');
mask = ncread(strcat(CROCO_path,'/',CROCO_file),'mask_rho');
mask(mask==0)=nan;               % Land is = 0, make nan
croco_top=croco_top.*mask;

% Read in the surface velocities

% Specify the days i will use for the plot, size is arbitrary

myDays = sort(randperm(31,3),'ascend');

% How many floats per release point

Fcount = 31;

% Subset the Grid domain to smaller domain used in paper 1

lat_min = -36;
lat_max = -28;
lon_min = 15;
lon_max = 20;

[~,idx_lat_min]=min(abs(croco_lat(1,:)-lat_min));
[~,idx_lat_max]=min(abs(croco_lat(1,:)-lat_max));
[~,idx_lon_min]=min(abs(croco_lon(:,1)-lon_min));
[~,idx_lon_max]=min(abs(croco_lon(:,1)-lon_max));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOWER RESOLUTION GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To avoind ailiasing, we need a lower resolution grid on which to have the
% positions interpolated or counted to

res = 1/12;  % Set the resolution in units of degrees

lat_low = lat_min:res:lat_max; 
lon_low = lon_min:res:lon_max;

[X_low, Y_low] = meshgrid(lon_low,lat_low);

X_low = X_low';
Y_low = Y_low';

% Interpolate the topography to the same grid resolution as the probability
% map for it to be overlayed

croco_lon = croco_lon(idx_lon_min:idx_lon_max,idx_lat_min:idx_lat_max);
croco_lat = croco_lat(idx_lon_min:idx_lon_max,idx_lat_min:idx_lat_max);
croco_top =  croco_top(idx_lon_min:idx_lon_max,idx_lat_min:idx_lat_max);
mask = mask(idx_lon_min:idx_lon_max,idx_lat_min:idx_lat_max);

croco_top_low = interp2(croco_lon',croco_lat',croco_top',X_low,Y_low); % Interp topography

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pick which floats you wish to process, ie
% 0(-1), -50, -100, -150, -200, -250

my_depth = -1;
freq = 4; % 6 hourly saved

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process zfiles in loop         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ymin = 2004;
Ymax = 2018;
myseason = 'Jul';
float_path = '/home/jono/Documents/PhD/CHAPTER_2/DATA/RESIDENCE/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MAIN = zeros(length(lon_low),length(lat_low),90,length(Ymin:Ymax));

for t = Ymin:Ymax
% Read in the grid and time dimension from the float file
    float_file = strcat('SBUS_floats_3km_',myseason,string(t),'.nc');
    disp(float_file)
    disp("Processing file")
    float_lon = Fto3D(ncread(strcat(float_path,'/',float_file),'lon'),Fcount);
    float_lat = Fto3D(ncread(strcat(float_path,'/',float_file),'lat'),Fcount);
    float_depth = Fto3D(ncread(strcat(float_path,'/',float_file),'depth'),Fcount);
    float_temp = Fto3D(ncread(strcat(float_path,'/',float_file),'temp'),Fcount);
    float_time = ncread(strcat(float_path,'/',float_file),'scrum_time');

    % Process time if it needs to be stacked
    float_time = Ftime(float_time,1200,Fcount);

    % Convert the time array
    time_float = datetime(1990,1,1) + seconds(float_time);

    % Calculate the count of each point on the model of a float vs a grid point
    % realtive to the lower resolution grid

    % Rules
    % 1) If a float is north of -28 S or west of 15 E, then the particle is 
    %    not considered in the domain
    % 2) The sum of the probabilities must always be equal to 1. This means we
    %    must account for point lost at the boundaries
    % 3) The lagrangain points are relative to the low resolution nodal points
    % 4) Count the points on the prob grid and divide by total num floats 

    OUTPUT = zeros(length(lon_low),length(lat_low),(length(float_time)-1)/freq,31); 
    for x = 1:size(float_temp,3)
        prob_grd = zeros(length(lon_low),length(lat_low),length(float_time)-1); % Low resolution probability grid
        for i = 2:length(time_float)     
            % STEP 1: find which floats are in the domain
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tmp_depth = squeeze(float_depth(:,:,x));
            ind_mydepth = find(tmp_depth == my_depth);
            tmp_lon = float_lon(ind_mydepth,i,x); % Assign to a temporary variable
            tmp_lat = float_lat(ind_mydepth,i,x); % Assign to a temporary variable
            ind_floats = find(round(tmp_lat,2) > max(lat_low) |...
                round(tmp_lon,2) < min(lon_low)); % Locate the indexes of fronts not in domain
            tmp_lon(ind_floats) = [];  % Remove missing
            tmp_lat(ind_floats) = [];  % Remove missing
            mycount = size(tmp_lat,1) - length(ind_floats); % Total number of floats in the domain to use for the probability
            % STEP 2: interpolate grid points onto prob_grd
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for j = 1:mycount
                [xx,yy] = mll2grid(tmp_lat(j),tmp_lon(j),Y_low,X_low);
                prob_grd(yy,xx,i-1) = prob_grd(yy,xx,i-1) + 1;
            end
            prob_grd(:,:,i-1) = prob_grd(:,:,i-1)./mycount;   % Probability field
            %disp(string(sum(prob_grd(:,:,i-1),'all'))); % Should be one
        end
        %    STEP 3: Average the prob_grd according to frequency of saved output to
        %    convert to a daily average
        prob_day = zeros(length(lon_low),length(lat_low),size(prob_grd,3)/freq);
        for k = 1:size(prob_grd,3)/freq
            prob_day(:,:,k) = mean(prob_grd(:,:,((k*freq)-freq)+1:k*freq),3,'omitnan');
        end
        OUTPUT(:,:,:,x) = prob_day;
    end

    MAIN(:,:,:,t+1-Ymin) = sum(OUTPUT,4,'omitnan');
end

% Quick note: due to the error on the first time-step it is not considered
% in the averaging process

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%
% PLOT
%%%%%%%%%%%%%%%%%

% The objectove is to show the evolution of the probability field in time
% as well as an average just for some completeness. Days 5, 20, 30, mean
% will be shown. Along with that, the probability field will be overlayed
% with the topography to see how particles are retained relative to the
% shelf

plt_var = mean(MAIN,4,'omitnan');
pos = [830     1   881   961 ];
mydepths = [200,300,500];
cmin = 0;
cmax = 0.2;%round(max(prob_day,[],'all'),3)/5;
levels = [cmin:0.02:cmax];
mydays = {'Day 5','Day 30','Day 60','Day 90'};
myinds = [5,30,60,90];
domain = [15,19;-33.5,-28];

f = figure;
for kk = 1:length(myinds)
    subplot(2,2,kk)
    m_proj('miller','long',[domain(1,:)]...
    ,'lat',[domain(2,:)]);
    m_contourf(X_low,Y_low,plt_var(:,:,myinds(kk)),levels,'ShowText','off');
    shading interp
    %colormap('jet')
    cmocean('matte',length(levels))
    m_gshhs_f('patch',[.7 .7 .7],'edgecolor','none');
    m_grid('box','fancy','linest','none','tickdir','out','fontsize',13);
    caxis([cmin cmax])
    cRange=caxis;
    hold on
    [C,h] = m_contour(X_low,Y_low,croco_top_low,mydepths,'ShowText','on');
    caxis(cRange)
    h.LineWidth = 3;
    h.Color = 'w';
    title(mydays{kk},'fontsize',18)
end
hold on
% Create colorbar
ca = colorbar('eastOutside');
ca.Position = ca.Position + 1e-10;
ca.Label.String = 'Mean cumulative particle density (%)';
ca.FontSize = 12;
caxis([cmin cmax]);

f.Position = pos;
set(gcf, 'InvertHardcopy', 'off')

print('-f1','Prob_map_Jul_sur','-dpng','-r600');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
