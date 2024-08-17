%%%% This function is intended to wrap the output of Roff from 2D to 3D. If
%%%% you release float at an interval, the 2D array is stacked with NaN. We
%%%% want to convert the two-D array such that the thrid dimension is the
%%%% date. myvar is the variable you want to process and Fcount is the
%%%% number of floats. Default settings is a time stepping of 6 hours with
%%%% a Dt = 1 day (Adjust if needed)

function MYOUTPUT = Fto3D(myvar,Fcount)

% Release parameters

Dt = 1;     % Interval between releases
%Fcount = Fcount; % Number of floats released per location
Fstep = 4;    % Frequency of saved outputs per day (6 hour time-step)
%disp(strcat('The number of release dates is... ',' ',string(Fcount),' at an interval of...',...
%    ' ', string(Dt),' ',' days')); 

% % Idea 1:
% % To make this easy, we are going to convert the known extent of the
% % release range to NaN values to account for the assignment of missing
% % values within ROFF
% 
% tmp = myvar(:,1:Fcount*Dt*Fstep);
% 
% % Convert to NaN missing data
% mis_var = 1.0e15;
% ind_mis = find(tmp == mis_var);
% 
% % Make NaN
% tmp(ind_mis) = NaN;
% 
% % Count the number of NaN values in each row
% %nancount = sum(isnan(tmp),2);

% IDEA 2: Using Nan values holds for small values of Fcount but does not
% hold thereafter for long runs. We therefore take a very hard coded
% apporach. Lets assume that we know Dt and Fstep such that we know where
% the indices will be
nancount = [1:Fstep*Dt:Fcount*Fstep*Dt];  % Cretae indice length
nancount = repmat(nancount',1,size(myvar,1)/Fcount);  %Stack array  
nancount = nancount(:); % Fromat
grr = (nancount./Fstep)+1;  
[B,I] = sort(grr,'ascend'); % Get indexs

% Number of floats
Numfloats = size(myvar,1)/length(unique(grr));

% Calculate the maximum tracking length such that all floats are tracked
% for the same length of time
tlen = size(myvar,2)-1;
% Add in case Fcount is > 1
if Fcount > 1
    nanlen = Fcount*Fstep;
else
    nanlen = 0;
end
track_length = tlen - nanlen;

%disp(['The number of day to track is',string(track_length/Fstep), 'days']);

% Format
a = myvar(I,:);
[r,c] = size(a);
nlay  = length(unique(grr));
out   = permute(reshape(a',[c,r/nlay,nlay]),[2,1,3]);

% Reformat the output of 3D array such that each day (third dimensions),
% corresponds to the length of the max tracking

MYOUTPUT = zeros(size(out,1),track_length+1,Fcount);
skp = 0;
for i = 1:Fcount
    tff = squeeze(out(:,[1,2+(skp*Fstep):track_length+1+(skp*Fstep)],i));
    MYOUTPUT(:,:,i) = tff;
    skp = skp+1;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
