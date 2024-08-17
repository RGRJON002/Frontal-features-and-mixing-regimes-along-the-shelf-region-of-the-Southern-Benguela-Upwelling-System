%%%% The intention of this function is to take daily averages of a 
%%%% float variable using the number of time-steps given by N. Need to preserve
%%%% the first index as this is the intialisation array info needed for
%%%% other processes


function [myoutput] = FU2day(myvar,N)
    % Average data every N-time-steps 
    szA = size(myvar(:,2:end)) ;
    % Now need to process lon, lat, depth, temperature and Days
    A = myvar(:,2:end);
    tmp = myvar(:,1);
    B = arrayfun(@(k) mean(A(:,k:min(szA(2),k+N-1)),2), 1:N:szA(2), 'un', 0) ;
    B = [B{:}] ;
    myoutput = [tmp, B];
end
