function [value] = lorenzian(x,x0,gamma)
%returns a value of lorenzian (another: Koshi distribution)
%   x0 - parametere of shift, gamma - scale parameter
    value = (1/pi) * gamma./((x-x0).^2 + gamma*gamma);
end

