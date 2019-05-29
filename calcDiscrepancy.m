function [eps] = calcDiscrepancy(calcSpectrum, expSpectrum, interval)
%UNTITLED14 считаем расхождение между двумя спектрами по МНК
%   Detailed explanation goes here
% %calcSpectrum(:,2) - expSpectrum(:,2)
% dW = inputSpectrum(2,1) - inputSpectrum(1,1);
% min_delta = inputSpectrum(1,1);
% l = floor((l_b-min_delta)/dW);
% r = ceil((r_b-min_delta)/dW);
left = interval(1);
right = interval(2);
eps = sqrt(sum((calcSpectrum(left:right,2) - expSpectrum(left:right,2)).^2));
end

