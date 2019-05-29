function [ratios, eps, G_w] = StepRatios(dRatio, hRatio, N, A, M, interval, parameters, ratios, j_coupling, J)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    fast_mode = 1;  % (1/0) - вкл/выкл радикальный спуск (смотрим только на знак производной)

    for i=1:length(ratios)
        if (ratios(i) == 0)
            continue
        end
        
        [G_ws, G_w] = calcSpectrum(N, A(:,1), M, interval, parameters, ratios, j_coupling, J);
        eps = calcDiscrepancyWeighted(G_w, A, interval);

        
        [ratios] = changeRatio(ratios, dRatio, i);
    
        [G_ws_n, G_w_new] = calcSpectrum(N, A(:,1), M, interval, parameters, ratios, j_coupling, J);
        eps_new = calcDiscrepancyWeighted(G_w_new, A, interval);
        
        dEps_dRatio = (eps_new - eps)/dRatio;
        
        [ratios] = changeRatio(ratios, -dRatio, i);
        if (fast_mode == 1)
            decrement = hRatio*sign(dEps_dRatio);
        else
            decrement = hRatio*dEps_dRatio;
        end
        
        if (ratios(i) - decrement >= 0 && ratios(i) - decrement <=1)  % второе условие специфично только для gamma
            [ratios] = changeRatio(ratios, -decrement, i);
        end
    end
end

function [ratios] = changeRatio(ratios, dRatio, i)
%changeRatio - изменяет i-ый вес на величину dRatio, чтобы сохранялась
%нормировка: sum(ratios)=1
    
    for k = 1:length(ratios)
        if (k == i)
            continue
        end
        % меняем веса на пропорциональную им величину, сохраняем нормировку
        ratios(k) = ratios(k) - dRatio*ratios(k)/(1-ratios(i));
    end
    ratios(i) = ratios(i) + dRatio; % x <- x0 + dx
    
    % учёт неточности вычислений
    delta_r = 1 - sum(ratios);
    ratios(1) = ratios(1) + delta_r;
end
