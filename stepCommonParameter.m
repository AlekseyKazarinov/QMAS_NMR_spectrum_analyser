function [parameters, eps, G_w] = stepCommonParameter(number, dValue, hValue, N, A, M, interval, parameters, ratios, j_coupling, J)
% stepCommonParameter performs one gradient descent step by common 
% parameter shared by all the spectrum components

% Input arguments:
% 1) number_line - number of a line (deprecated: accessible only 1 or 2) 
% 2) number - номер параметра в списке parameters соотв. линии:
%           eta -> 1, Chi -> 2, delta -> 3, sigma -> 4, gamma -> 5
% 3) sp_mesh - grid (mesh) for spectrum counts (counts of chemical shift)
% 4) dValue, hValue - increasing and step of related value argument
% 5) interval = [l,r] - outer indices of the array where error will be 
%                       calculated
% 6) parameters - sets of parameters for all the components consisting 
%                       the whole spectrum
% 7) ratios - ratios between intencities of spectrum components
% 8) j_coupling - 
% 9) J - 
  
% Special implicit agrument:
fast_mode = 1;  % (1/0) - on/of mode of fast gradient descent (
                % (in this case we just takes into accout sign of derivative)
    

    [G_ws, G_w] = calcSpectrum(N, A(:,1), M, interval, parameters, ratios, j_coupling, J);
    eps = calcDiscrepancyWeighted(G_w, A, interval);
    
    % increasing of shared parameter
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
    
    % secont condition is specified for gamma parameter only
    if (parameters(1,number) - decrement > 0 && parameters(1,number) - decrement <=100)  
        parameters(:,number) = parameters(:,number) - decrement;
    end
end

