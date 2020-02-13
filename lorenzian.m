function [value] = lorenzian(x, x0, gamma)
% returns a value of lorenzian (another defition: Koshi distribution)

% Input arguments:
% 1) x0 - parametere of shift 
% 2) gamma - scale parameter


    value = (1/pi) * gamma./((x-x0).^2 + gamma*gamma);
end

