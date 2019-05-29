function [integral] = integrate(f, interval)
%UNTITLED3 Summary of this function goes here
%   f - input spectrum, interval = [a,b]
%[L,R] = convertToIndices(interval, f);
L = interval(1);
R = interval(2);
dW = f(2,1) - f(1,1);
integral = (f(L,2) + f(R,2))*0.5 + sum(f(L+1:R-1,2));
integral = integral*dW;
end

