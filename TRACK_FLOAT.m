%%%% This function is designed to calculate the residence time of floats in
%%%% the SBUS from the ROFF outputs. Originally, the input was the Roff file
%%%% but on revision, the best is to input the lon, lat and depth data to
%%%% easily facilitate the function running in a loop. The CROCO_file is
%%%% required to set the limitations on the file and set the isobath
%%%% limits. third is a boolean input of 1 to make a plot and any
%%%% value thereafter to not plot the fits. 

function OUTPUT = TRACK_FLOAT(lon,lat,depth,CROCO_file,makeplots)

addpath /home/jono/Documents/MATLAB/ezyfit2.44/ezyfit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read in the grid and time dimension from the float file

float_lon = lon;
float_lat = lat;
float_depth = depth;

% float_lon = ncread(float_file,'lon');
% float_lat = ncread(float_file,'lat');
% float_depth = ncread(float_file,'depth');
% float_time = ncread(float_file,'scrum_time');
% 
% % Convert the time array
% time_float = datetime(1990,1,1) + seconds(float_time);

% Read in the grid dimension of one of the CROCO_files

croco_lon = ncread(CROCO_file,'lon_rho');
croco_lat = ncread(CROCO_file,'lat_rho');

% Read in the topography data

croco_top = ncread(CROCO_file,'h');
mask = ncread(CROCO_file,'mask_rho');
mask(mask==0)=nan;               % Land is = 0, make nan
croco_top=croco_top.*mask;

% Find the long and lat coordinate pairs for the desired isobath

mybath = 300;   % Set the depth for the desired isobath delimeter

lon_tmp = croco_lon(:,1);
lat_tmp = croco_lat(1,:);

for ii = 1:size(croco_top,2)
 ind = min(find(croco_top(:,ii)<=mybath));
 if isempty(ind)
    data(:,ii) = NaN;
 else
    data(:,ii) = lon_tmp(ind);
 end
end

loc_300 = [data' lat_tmp'];

%% Algorithm to calculate residence times for the different floats initiated at the different depth

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RULES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To calculate the residence time, we need to create criteria for which
% floats are considered to be on the SBUS shelf. These outputs need to be
% saved after every time-step and then plotted agianst time to achieve some
% sort of exponential decay function.

% 1) If a float goes south of Cape Columbine (32.82 S)
% 2) If a float goes north of 28 S 
% 3) If a float spends one day beyond the 300 m isobath

% Do some pre-processing for the criteria

%%%%%%%%%%%%%%%
% Criteria 1
%%%%%%%%%%%%%%%

lat_min = -32.82;
LTmin_index = min(find(lat_tmp>=lat_min));
lat_min = lat_tmp(LTmin_index);

%%%%%%%%%%%%%%%%
% Criteria 2
%%%%%%%%%%%%%%%%

% Particles move outside the boundary region of the model when their
% coordinate value is 1.0e15

mis_var = 1.0e15;
ind_mis = find(float_lat == mis_var);

% Make NaN

float_lon(ind_mis) = NaN;
float_lat(ind_mis)= NaN;
float_depth(ind_mis) = NaN;

%%%%%%%%%%%%%%%
% Criteria 3
%%%%%%%%%%%%%%%

% We saved data every 6 hours, so...

time_freq = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pick which floats you wish to process, ie
% 0(-1), -50, -100, -150, -200, -250

my_depth = flipud(unique(float_depth(:,1)))'; % All depths that are used 
NPD = size(float_depth,1)/length(my_depth); % Once off use as will use ind_depth

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

float_ID = NaN(NPD,size(float_depth,2),length(my_depth));
for ii = 1:length(my_depth)
    ind_mydepth = find(float_depth == my_depth(ii));
    for jj = 1:size(float_depth,2)
        for kk = 1:length(ind_mydepth)
            LT_index = min(find(loc_300(:,2)>=float_lat(ind_mydepth(kk),jj)));
            if float_lat(ind_mydepth(kk),jj) < lat_min 
               float_ID(kk,jj:end,ii) = ind_mydepth(kk);
            elseif isnan(float_lat(ind_mydepth(kk),jj))
               float_ID(kk,jj:end,ii) = ind_mydepth(kk);
            elseif float_lon(ind_mydepth(kk),jj) < loc_300(LT_index,1)
                ind_count = jj + time_freq;
                if ind_count > size(float_depth,2)
                    float_ID(kk,jj:end,ii) = ind_mydepth(kk);
                else
                    for ll=1:time_freq
                        LT_index = min(find(loc_300(:,2)>=float_lat(ind_mydepth(kk),jj+ll-1)));
                        if isempty(LT_index)
                            check(:,ll) = NaN;
                        else
                            check(:,ll) = float_lon(ind_mydepth(kk),jj+ll-1) - loc_300(LT_index,1);
                        end
                    end
                    check_ind = find(check < 0);
                    clear check
                    if length(check_ind == time_freq)
                       float_ID(kk,jj:end,ii) = ind_mydepth(kk);
                    end
                end
            else
               %float_ID(kk,jj,ii) = NaN;
            end 
        end
    end
end
           
%%%%%%%%%%%%%%%%%%%%%
% POST-PROCESS
%%%%%%%%%%%%%%%%%%%%%

% Get rid of the first time-step
float_ID(:,1,:) = [];   

% Reassign the main float_ID array

ALL_ID = float_ID;

% Loop over all the floats 
% We now want to calculate how many floats were in the domain
% Remeber, data is saved every six hours 

ALL_data = NaN(length(my_depth),size(ALL_ID,2));
for pp = 1:length(my_depth)
    tot = length(ind_mydepth);
    record = NaN(1,length(ind_mydepth));
    float_ID = squeeze(ALL_ID(:,:,pp));
    for xx=1:size(float_ID,2)
        ind_num = ~isnan(float_ID(:,xx));
        if isempty(float_ID(ind_num,xx))
           record(xx) = 0;
        else
            myfloat = float_ID(ind_num,xx);
            record(xx) = length(myfloat);
        end
    end

    for yy = 1:length(record)-1
        mytotal(yy) = record(yy+1) - record(yy);
    end
    mytotal = [0 mytotal];

    for nn =1:length(mytotal)
        tot = tot-mytotal(nn);
        mydata(:,nn) = tot;
    end
    ALL_data(pp,:) = mydata;
    clear mytotal mydata
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot data to aquire curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Represent function as a piecewise function, assuming that we have a
% straight line corresponding to the period of no particles being lost and
% then an exponential decay function

ALL_floats = NaN(length(my_depth),size(float_ID,2)/time_freq);
for gg = 1:length(my_depth)
    mydata = ALL_data(gg,:);
    a = mydata(1,:);
    n = time_freq;
    b = arrayfun(@(i) mean(a(i:i+n-1)),1:n:length(a)-n+1)'; % the averaged vector
    ALL_floats(gg,:) = b;
end

clear b
% Just round the value 

ALL_floats = round(ALL_floats);

% Determine for how many days there is no change in the number of floats 

Lambda = NaN(1,length(my_depth));
Rsquare = NaN(1,length(my_depth));
residence = NaN(1,length(my_depth));  % Units of days

for bb = 1:size(ALL_floats,1)
    b = ALL_floats(bb,:);
    ind_tot = max(find(b == length(ind_mydepth)));
    
    if ~isempty(ind_tot) 
        plot(ind_tot:length(b),b(ind_tot:end)./length(ind_mydepth),'o')
    else
        plot(1:length(b),b./length(ind_mydepth),'o')
    end

    %plot(1:length(b),b./length(ind_mydepth),'o')   % Original
    % showfit exp
    undofit  % deletes the previous fit
    curve = showfit('f(t)=a*exp(lambda*t)','fitlinewidth',2,'fitcolor','red'); 
    Lambda(bb) = curve.m(2);
    Rsquare(bb) = curve.r; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RESULTS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Our residence time is computed for the value 1/e

    %res_time = 1/abs(curve.m(2));   % Disregard intercept and piecwise function
    %res_time = ind_tot + (1+log(curve.m(1)))/abs(curve.m(2));  % Acknowledge intercept

    if ~isempty(ind_tot)
        res_time = ind_tot + 1/abs(curve.m(2));  % Disregard intercept
    else
        res_time = 1/abs(curve.m(2));
    end
    residence(bb) = res_time; 
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combine all info for a run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OUTPUT = [my_depth; Lambda; Rsquare; residence];  % This is the summary of all the result

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot desired curves from the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if makeplots == 1
    for ff = 1:length(my_depth)
        b = ALL_floats(ff,:);
        ind_tot = max(find(b == length(ind_mydepth)));
        figure
        if ~isempty(ind_tot) 
            plot(ind_tot:length(b),b(ind_tot:end)./length(ind_mydepth),'o')
        else
            plot(1:length(b),b./length(ind_mydepth),'o')
        end
    %plot(1:length(b),b./length(ind_mydepth),'o')   % Original
    % showfit exp
        undofit  % deletes the previous fit
        curve = showfit('f(t)=a*exp(lambda*t)','fitlinewidth',2,'fitcolor','red');
        title(strcat('Best fit for',' ',string(my_depth(ff)),'m'));
        xlabel('Days')
        ylabel('Ratio')
    end
end

end