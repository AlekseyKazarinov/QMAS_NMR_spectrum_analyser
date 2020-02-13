function [outputSpectrum] = normalizeSpectrum(inputSpectrum, interval)
%normalizeSpectrum - returns normilized intensity of inputSpectrum
    outputSpectrum = zeros(length(inputSpectrum), 2);
    outputSpectrum(:,1) = inputSpectrum(:,1);
    
    % ���������� �� �������
    area = integrate(inputSpectrum, interval);
    outputSpectrum(:,2) = inputSpectrum(:,2) ./ area;
    
    % ���������� �� ���������
    %maxIntensity = max(inputSpectrum(:,2)); 
    %outputSpectrum(:,2) = inputSpectrum(:,2) ./ maxIntensity;
end

