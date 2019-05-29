function [outputSpectrum] = convolution(inputSpectrum, sigma, gamma, alpha)
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here
    dW = inputSpectrum(2,1) - inputSpectrum(1,1);
    %leftW = inputSpectrum(1,1);
    %rightW = inputSpectrum(length(inputSpectrum),1);
    wRef = 105.84; % MHz
    sigma = (sigma)./(wRef);  % перевод Hz to ppm
    gamma = (gamma)./(wRef);
    k = ceil(sigma/dW);
    n_width = 100; % должно быть чётным и не меньше 6
    norm = linspace(1, n_width*k, n_width*k);
    for i = 1:length(norm)
       norm(i) = alpha*normpdf(norm(i)*dW, 0.5*n_width*k*dW, sigma) + ...
           (1-alpha)*lorenzian(norm(i)*dW, 0.5*n_width*k*dW, gamma);
    end
    norm = norm ./ max(norm);
    outputSpectrum = zeros(length(inputSpectrum), 2);
    outputSpectrum(:,1) = inputSpectrum(:,1);
    % т.к. операция conv (свёртка) увеличивает размер массива: n+k-1:
    %inputSpectrum(round(n_width*k*0.5):length(inputSpectrum)-round(n_width*k*0.5), 2)'
    outputSpectrum(:,2) = conv(inputSpectrum(round(n_width*k*0.5):length(inputSpectrum)-round(n_width*k*0.5), 2), norm);
end

