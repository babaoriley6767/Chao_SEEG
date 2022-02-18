function y = Calc_RMS(a,time_s,time_e)
% Function to calculate root mean square

% Input: A

%        time_s: start time
%        time_e: end time 

% Output: Raw RMS value
    x=0;
    for i=time_s:time_e
        x = x+a(i)^2;
    end
    y=x/(time_e-time_s+1);
    y=sqrt(y);   
end

