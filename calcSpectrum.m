function [G_ws, G_w] = calcSpectrum(N, sp_mesh, M, interval, params, ratios, j_coupling, J)
% calculates a complex spectrum which can be represented as composinion of
% several different specrum lines (to be summed together). Parameters of
% each spectrum are listed in two-dimensional array in params, that
% represents [parameters1, parameters2, ..., parametersN], where
% parameters# are setting up a spectrum line with # ordinal number.
% Each of line has relative integral intensity defined from ratios:
% ratios = [r1, r2, r3, ..., rN] so that:
% r1 + r2 + ... + rk + ... + rN == 1, where rk - an element from ratios.

% Input arguments:
% 1) N - a number of points on X-axis which will be used for calculating
%        the whole spectrum
% 2) sp_mesh - 
% 3) M - 
% 4) interval - 
% 5) params - 
% 6) ratios - 
% 7) j_coupling - 
% 8) J - 

    
    % Assertion that the condition of normalization is performed:
    if (sum(ratios) < 0.999)
        stop       % in case of the violation
    end
    
    n_lines = length(ratios);  % count score number of spectrum lines
    G_ws = zeros(N, 2, n_lines);
        
    % calculating the first one (It is just specific case)
    G_w1 = calcLine(N,sp_mesh,M,interval, params(1,:), j_coupling, J);
    G_w1(:,2) = ratios(1)*G_w1(:,2);
    G_w = G_w1;
    G_ws(:,:,1) = G_w1;
    % Then add rest lines using normalization coefficients in ratios:
    for k = 2:n_lines
        G_wk = calcLine(N,sp_mesh,M,interval, params(k,:), 0, J);
        G_wk(:,2) = ratios(k)*G_wk(:,2);
        G_w(:,2) = G_w(:,2) + G_wk(:,2);
        G_ws(:,:,k) = G_wk;
    end
end

