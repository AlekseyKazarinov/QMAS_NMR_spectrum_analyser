function [parameters, eps, G_w] = stepLineParameter(number_line, number, dValue, hValue, N, A, M, interval, parameters, ratios, j_coupling, J)
%stepLineParameter - градиентный спуск вдоль данного параметра отдельной
%   линии
%   number_line - номер линии (только 1 или 2)
%   number - номер параметра в списке parameters соотв. линии:
%           eta -> 1, Chi -> 2, delta -> 3, sigma -> 4, gamma -> 5
%   sp_mesh - сетка для спектра (отсчёты хим. сдвига)
%   dValue, hValue - приращение и шаг аргумента
%   interval = [l,r] - задаются крайние индексы массива, где будет
%           считаться ошибка
%   parameters - наборы параметров для всех линий в спектре
%   ratios - соотношения между интенсивностями линий
    
    fast_mode = 1;  % (1/0) - вкл/выкл радикальный спуск (смотрим только на знак производной)
    
    [G_ws, G_w] = calcSpectrum(N, A(:,1), M, interval, parameters, ratios, j_coupling, J);
    eps = calcDiscrepancyWeighted(G_w, A, interval);
    if (ratios(number_line) > 0)  % для экономии ресурсов
        parameters(number_line,number) = parameters(number_line,number) + dValue;

        [G_ws_n, G_w_new] = calcSpectrum(N, A(:,1), M, interval, parameters, ratios, j_coupling, J);
        eps_new = calcDiscrepancyWeighted(G_w_new, A, interval);
        dEps_dValue = (eps_new - eps)/dValue;
        parameters(number_line,number) = parameters(number_line,number) - dValue;
        
        if (fast_mode == 1)
            decrement = hValue*sign(dEps_dValue);
        else
            decrement = hValue*dEps_dValue;
        end 
        
        if (parameters(number_line,number) - decrement >= 0)
            parameters(number_line,number) = parameters(number_line,number) - decrement;
        end
    end
end
