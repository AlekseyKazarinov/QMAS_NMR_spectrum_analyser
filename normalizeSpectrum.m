function [outputSpectrum] = normalizeSpectrum(inputSpectrum, interval)
% returns normilized intensity of inputSpectrum on interval

% Input arguments:
% 1) inputSpectum - a specrum line (represented by array) to be normalized
% 2) interval = [a, b]


    outputSpectrum = zeros(length(inputSpectrum), 2);
    outputSpectrum(:,1) = inputSpectrum(:,1);
    
    % normalization by area
    area = integrate(inputSpectrum, interval);
    outputSpectrum(:,2) = inputSpectrum(:,2) ./ area;
    
    % normalization by intencity (amplitude)  % temporary unused
    %maxIntensity = max(inputSpectrum(:,2)); 
    %outputSpectrum(:,2) = inputSpectrum(:,2) ./ maxIntensity;
end

