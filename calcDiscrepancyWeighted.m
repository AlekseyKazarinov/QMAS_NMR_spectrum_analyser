function [eps] = calcDiscrepancyWeighted(calcSpectrum, expSpectrum, interval)
% Считает расхождение (ошибку) между двумя спектрами,
% что вес пропорционален установленной величине (интенсивности второго спектра)
% m - медианое значение спектра в ppm
left = interval(1);
right = interval(2);

variant = 2; 
% 0 - среднеквадратичное отклонение
% 1 - пользовательский вариант
% 2 - вес прямо пропорционален интенсивности эксп. спектра

if (variant == 0)
    eps = sqrt(sum((calcSpectrum(left:right,2) - expSpectrum(left:right,2)).^2));
%     plot((calcSpectrum(left:right,2) - expSpectrum(left:right,2)).^2);
%     stop
end

% вес прямо пропорционален расстоянию от отсчёта до медианы 
if (variant ==1)
    m = medium(calcSpectrum, interval);

    N = length(expSpectrum(:,1));
    w_min = expSpectrum(1,1);
    w_max = expSpectrum(N,1);
    dw = (w_max - w_min)/N;

    m_index = floor((m - w_min)/dw);
    %abs(calcSpectrum(left:right,1) - m) * 0.2
    weight = 0.5;
    eps = sum( abs( (calcSpectrum(left:right,2) - expSpectrum(left:right,2)).*  abs(calcSpectrum(left:right,1) - m) * weight ));

% plot(( (calcSpectrum(left:right,2) - expSpectrum(left:right,2)).*  abs(calcSpectrum(left:right,1) - m) * weight ).^2)
% stop
end

% вес прямо пропорционален высоте отсчёта
if (variant == 2)
    eps = sqrt(sum( ((calcSpectrum(left:right,2) - expSpectrum(left:right,2)).^2) .* abs(expSpectrum(left:right,2)) ));
    %eps = sum(abs( (calcSpectrum(left:right,2) - expSpectrum(left:right,2)).*  abs(maxIntensity./calcSpectrum(left:right,2)) ));
    % визуализация веса для подсчёта ошибки разных отсчётов
    %plot( ((calcSpectrum(left:right,2) - expSpectrum(left:right,2)).^2) .* expSpectrum(left:right,2));
    %stop
end
end


