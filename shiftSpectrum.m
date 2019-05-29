function [outputSpectrum] = shiftSpectrum(inputSpectrum, shift)
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here
    outputSpectrum = zeros(length(inputSpectrum),2);
    dW = inputSpectrum(2,1) - inputSpectrum(1,1);
    outputSpectrum(:,1) = inputSpectrum(:,1);
    N = length(inputSpectrum);
    if (shift >=0)
        n_shift = round(shift/dW);
        outputSpectrum(1:n_shift,2) = zeros(n_shift,1); 
        outputSpectrum(n_shift+1:N,2) = inputSpectrum(1:N-n_shift,2);
%         outputSpectrum(:,1) = outputSpectrum(:,1) + shift; % здесь косяк
    else
        n_shift = round(-shift/dW);
        outputSpectrum(N-n_shift+1:N,2) = zeros(n_shift,1);
        outputSpectrum(1:N-n_shift,2) = inputSpectrum(n_shift+1:N,2);
    end
end

