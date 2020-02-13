function [ratios, eps, G_w] = StepRatios(dRatio, hRatio, N, A, M, interval, parameters, ratios, j_coupling, J)
% this function makes fitting spectrum lines that form the complex spectrum
% line by regulation of ratio coefficiens of each line.
% The function realizes just one step of fitting per each call.
% The fitting implements gradient descent by ratios of each component.

% Input arguments:
% 1) dRatio - increment of ration for calculating the derivative
% 2) hRatio - value of ratio step
% 3) N - 
% 4) A - 
% 5) M - 
% 6) interval - 
% 7) parameters - 
% 8) ratios - 
% 9) j_coupling - 
% 10) J - 

% Special implicit agrument:
fast_mode = 1;  % (1/0) - on/of mode of fast gradient descent
% (in this case we just takes into accout sign of derivative)


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
% changes i-th ratio coefficient by dRatio value so that
% the normalization condition is performed
% normalization: sum(ratios)=1
    
% Input arguments:
% 1) ratios - an array of ratio coefficients of each component (spectrum
%             line)
% 2) dRatio - increase of ratio for calcucating derivative


    for k = 1:length(ratios)
        if (k == i)
            continue
        end
        % changes the weights веса by a proportional value with saving 
        % normalization:
        ratios(k) = ratios(k) - dRatio*ratios(k)/(1-ratios(i));
    end
    ratios(i) = ratios(i) + dRatio; % x <- x0 + dx
    
    % calculate error:
    delta_r = 1 - sum(ratios);
    ratios(1) = ratios(1) + delta_r;
end
