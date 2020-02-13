function [outputSpectrum] = convolution(inputSpectrum, sigma, gamma, alpha)
% returns outputSpectrum which is a convolution result of inputSpectrum and 
% bell-shaped line expressed by Lorenzian and Gaussian fuctions.
% It is need for account of broadening phanomena for real spectum lines
% from each monocrystal system composing polycrystal to be researched.

% Input arguments:
% 1) inputSpectum - original spectrum line to be fixed
% 2) sigma - exponenta parameter for Gaussian
% 3) gamma - broadening parameter for Lorenzian
% 4) alpha - ratio between Gaussian and Lorenzian


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
    % because convolunion operation increases size of the array: n+k-1:
    % inputSpectrum(round(n_width*k*0.5):length(inputSpectrum)-round(n_width*k*0.5), 2)'
    outputSpectrum(:,2) = conv(inputSpectrum(round(n_width*k*0.5):length(inputSpectrum)-round(n_width*k*0.5), 2), norm);
end

