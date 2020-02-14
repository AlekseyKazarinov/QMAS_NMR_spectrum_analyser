function [eps] = calcDiscrepancy(calcSpectrum, expSpectrum, interval)
% returns an error of calculation between theoretical spectrum and
% experimental one using LS method (Least Squares)

% Input arguments:
% 1) calcSpectrum - a theoretically expected spectrum line
% 2) expSpectrum - a specrum obtained from experimental data
% 3) interval - the interval on which calculation is performed

    % Deprecated:
    % %calcSpectrum(:,2) - expSpectrum(:,2)
    % dW = inputSpectrum(2,1) - inputSpectrum(1,1);
    % min_delta = inputSpectrum(1,1);
    % l = floor((l_b-min_delta)/dW);
    % r = ceil((r_b-min_delta)/dW);
    
    left = interval(1);
    right = interval(2);
    eps = sqrt(sum((calcSpectrum(left:right,2) - expSpectrum(left:right,2)).^2));
end

