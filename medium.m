function [m] = medium(f, interval)
% определяет медианное значение спектра (точка, интеграл справа = инт. слева)
%[L,R] = convertToIndices(interval, f);
L = interval(1);
R = interval(2);
%fullArea = (b-a)/(2*n)*(f(L) + f(R) + 2*sum(f(L+1:R-1,2)));
fullArea = (f(L,2) + f(R,2))*0.5 + sum(f(L+1:R-1,2));

area = 0;
k = L;
while (2*area < fullArea)
    area = (f(k,2)+f(k+1,2))/2 + area;
    k = k + 1;
end
m = f(k,1);
end

