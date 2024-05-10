%%%% Script to extract depth of a chosen variable

addpath /home/jono/CROCO/croco_tools/UTILITIES/m_map1.4h
addpath /media/data/DATASETS_CROCOTOOLS/m_map1.4f
addpath /home/jono/Documents/MATLAB/CMOCEAN
addpath '/home/jono/Documents/MATLAB/CDT-master/cdt'
addpath '/home/jono/Documents/MATLAB/CDT-master/cdt/cdt_data'
addpath '/home/jono/Documents/MATLAB/tight_subplot'
addpath('/usr/local/MATLAB/R2020a/toolbox/matlab/imagesci/')

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
% Read in topo
h = ncread(file,'h');
h = h.*mask;
% Define my region, will focus on the coastal domain
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

CROCO_lat = CROCO_lat(idx_lon_min:idx_lon_max,idx_lat_min:idx_lat_max);
CROCO_lon = CROCO_lon(idx_lon_min:idx_lon_max,idx_lat_min:idx_lat_max);

%%%%%%%%%%%%%%
% VARIABLES
%%%%%%%%%%%%%%

var = 'w';
type = 'r'; % rho grid
depth = -10;  % Depth in [m]

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
        % Read in the vertical velocity data
        addpath('/usr/local/MATLAB/R2020a/toolbox/matlab/imagesci/')
        w = ncread(file,'w');
        for k = 1:length(time)
            w_slice(:,:,k) = get_hslice(convertStringsToChars(file),...
                convertStringsToChars(file),var,1,depth,type);
        end
        w_slice = w_slice(idx_lat_min:idx_lat_max,idx_lon_min:idx_lon_max,:);
        %  Store arrays
        disp('Storing w and time')
        TIME = cat(1,TIME,time);
        W = cat(3,W,w_slice);
    else
        disp(strcat('No data for',file))
    end
    clear w w_slice
    end
end

% Format time

[~, ia, ~] = unique(TIME);
W = W(:,:,ia);
TIME = TIME(ia);

% Save to a structure

% Divide W into chunks

Chunk = length(TIME)/4;

W_1 = W(:,:,1:Chunk);
W_2 = W(:,:,Chunk+1:Chunk*2);
W_3 = W(:,:,((Chunk*2)+1):Chunk*3);
W_4 = W(:,:,((Chunk*3)+1):end);

for ii = 1:4
    save(strcat('W_',string(ii),'.mat'),strcat('W_',string(ii)),'TIME',...
        'CROCO_lat','CROCO_lon')
end


