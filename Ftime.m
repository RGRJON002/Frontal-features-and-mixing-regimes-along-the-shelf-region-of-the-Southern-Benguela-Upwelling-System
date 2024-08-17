%%%% Function to transform the time variable so that it is stacked if
%%%% Fcount > 1. Will need to know the number of floats as well as the
%%%% Fcount variable in order to do the processing

function [float_time] = Ftime(time,numfloats,Fcount)
    A = repmat(time',numfloats*Fcount,1);
    tmp = Fto3D(A,Fcount);
    float_time = squeeze(mean(tmp,1));
end