function [l,r] = convertToIndices(interval, inputSpectrum)
% this function calculates of indices interval bounds within which spectrum
% will be fitted.

% Input arguments:
% 1) interval - an array representiing [w_min, w_max], w - frequency (ppm)


    l_bound = interval(1);
    r_bound = interval(2);
    dW = inputSpectrum(2,1) - inputSpectrum(1,1);
    min_delta = inputSpectrum(1,1);
    l = floor((l_bound-min_delta)/dW);
    r = ceil((r_bound-min_delta)/dW);
end

