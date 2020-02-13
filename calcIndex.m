function [index] = calcIndex(w, w_min, w_max, N)
%getIndex - возвращает соответствующий индекс в массиве G_w
% значения y будут попадать в [Ymin, Ymax]
% Ymin, Ymax рассчитываются заранее для конкретного случая
% N - размер массива-гистограммы
dw = (w_max - w_min)/N;
w_n = w - w_min;
n = floor(w_n/dw);
index = n;
    %Yn = y - Ymin;
    %n = round(Yn/dY);
    %index = n + 1;
end

