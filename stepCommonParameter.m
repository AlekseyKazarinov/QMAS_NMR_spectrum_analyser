function [parameters, eps, G_w] = stepCommonParameter(number, dValue, hValue, N, A, M, interval, parameters, ratios, j_coupling, J)
%stepCommonParameter - градиентный спуск по общему параметру всех линий
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
    
    % приращение по общему параметру всех линий
    parameters(:,number) = parameters(:,number) + dValue; % x <- x0 + dx
    
    [G_ws_n, G_w_new] = calcSpectrum(N, A(:,1), M, interval, parameters, ratios, j_coupling, J);
    eps_new = calcDiscrepancyWeighted(G_w_new, A, interval);
    dEps_dValue = (eps_new - eps)/dValue;
    parameters(:,number) = parameters(:,number) - dValue;  % x0 <- x - dx
    
    if (fast_mode == 1)
        decrement = hValue*sign(dEps_dValue);
    else
        decrement = hValue*dEps_dValue;
    end 
    
    % второе условие специфично только для gamma
    if (parameters(1,number) - decrement > 0 && parameters(1,number) - decrement <=100)  
        parameters(:,number) = parameters(:,number) - decrement;
    end
end

