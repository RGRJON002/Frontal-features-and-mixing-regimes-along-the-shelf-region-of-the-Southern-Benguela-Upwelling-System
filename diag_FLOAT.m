%%%% Diagnostic function to derive frontal variables such as component
%%%% velocities, temperature gradients as well as days since initialisation
%%%% Takes the float output file of ROFF and the time over which files are
%%%% saved in hours for time_step. Note that Fto3D runs inside this script

function [U,V,W,Td,Days] = diag_FLOAT(file,time_step,Fcount)

    float_lon = Fto3D(ncread(file,'lon'),Fcount);
    float_lat = Fto3D(ncread(file,'lat'),Fcount);
    float_depth = Fto3D(ncread(file,'depth'),Fcount);
    float_temp = Fto3D(ncread(file,'temp'),Fcount);
    float_time = ncread(file,'scrum_time');

    % Convert the time array
    time_float = datetime(1990,1,1) + seconds(float_time);
    
    % Compute the zonal and meridional components for velocity and
    % temperarure gradient
 
    for k = 1:size(float_lat,3)  
        for i = 1:size(float_lat,2)-1
            for j = 1:size(float_lat,1)
                distv = spheric_dist(float_lat(j,i,k),float_lat(j,i+1,k),float_lon(j,i,k),float_lon(j,i,k)); % Meridional
                distu = spheric_dist(float_lat(j,i,k),float_lat(j,i,k),float_lon(j,i,k),float_lon(j,i+1,k)); % Zonal
                % Meridional velocity
                if float_lat(j,i,k) < float_lat(j,i+1,k)
                    v(j,i,k) = distv/(time_step*60*60);
                else
                    v(j,i,k) = -distv/(time_step*60*60);
                end
                % Zonal velocity
                if float_lon(j,i,k) < float_lon(j,i+1,k)
                    u(j,i,k) = distu/(time_step*60*60);
                else
                    u(j,i,k) = -distu/(time_step*60*60);
                end
                % Calculate the temporal temperature gradients
                diff = float_temp(j,i+1,k)-float_temp(j,i,k);
                Td(j,i,k) = diff/((time_step*60*60)/(86400)); % deg C/day
                % This is a mistake. You cannot compute the spatial
                % gradient from taking the time derivative of temperature. 
%                 Tx(j,i,k) = diff/distu;  
%                 Ty(j,i,k) = diff/distv;
                % Do the vertical velocity
                tmp = float_depth(j,i+1,k)-float_depth(j,i,k);
                w(j,i,k) = tmp/((time_step*60*60)/(86400));  % m/day 
            end
        end
    end
    
    % Convert temperature gradient to deg/km
    
%     Tx = Tx.*1000;
%     Ty = Ty.*1000;
    
    % Pad dimensions
    
    u = [u(:,1,:), u];
    v = [v(:,1,:), v];
    w = [w(:,1,:), w];
    Td = [Td(:,1,:), Td];
%     Tx = [Tx(:,1,:), Tx];
%     Ty = [Ty(:,1,:), Ty];
        
    % Assign. velocity are in m/s and m/day (vertical)
    
    U = u;
    V = v;
    W = w;
    
    % Days since initialisation
    % We can be fancy and compute it using date formats but knowing the
    % freq of the saved variable and length will give and indication.
    
    Duration = (size(float_lat,2) -1)/(24/time_step); % Number of days floats are tracked 
    freq = 24/time_step;
    Days = repmat(1:Duration,[freq 1]);
    Days = Days(:)';
    
    % Pad 
    Days = [Days(1),Days];
    
end
    
  