function [index] = calcIndex(w, w_min, w_max, N)
% returns index correspoding to a serial
% number of section in which argument w gets into if space between w_min
% and w_max is divided into N sections.
    
% input arguments:
% 1) w - original frequency
% 2) w_min - minimal frequency
% 3) w_max - maximal frequency
% 4) N - number of segments into which section [w_min, w_max] will be
%        splitted 

    dw = (w_max - w_min)/N;
    w_n = w - w_min;
    n = floor(w_n/dw);
    index = n;
end

