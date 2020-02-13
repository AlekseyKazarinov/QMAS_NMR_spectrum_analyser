function [integral] = integrate(f, interval)
% just returns area under line given from function f on the interval.

% Input arguments:
% 1) f - input spectrum 
% 2) interval = [a,b], where a and b are just numeric


%[L,R] = convertToIndices(interval, f);  % deprecated operation
L = interval(1);
R = interval(2);
dW = f(2,1) - f(1,1);
integral = (f(L,2) + f(R,2))*0.5 + sum(f(L+1:R-1,2));
integral = integral*dW;
end

