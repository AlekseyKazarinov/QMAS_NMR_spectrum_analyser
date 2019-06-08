function [G_w] = calcLine(N, sp_mesh, M,interval, parameters, j_coupling, J)
%calcLine - численное построение одной спектральной линии по заданным
%параметрам
%   parameters - набор параметров для данной линии, это массив вида:
%       [eta, Chi, delta, sigma, gamma, alpha];

% Параметры экспериментальной системы:
    wL = 105.84; % MHz (Ларморова частота)
    I = 3/2;     % спин ядра (для 23Na)


%parameters
    eta = parameters(1);
    Chi = parameters(2);
    delta = parameters(3);
    sigma = parameters(4);
    gamma = parameters(5);
    alpha = parameters(6);

    G_w = zeros(N,2);  % гистограмма-спектр, рассчитывается
    %[Ymin, Ymax] = getYMinMax();
    %dY = (Ymax - Ymin)/N;
    %G_y(:,1) = linspace(Ymin+dY*0.5, Ymax-dY*0.5, N);
    G_w(:,1) = sp_mesh;
    w_min = sp_mesh(1,1);
    w_max = sp_mesh(N,1);
    my_pi = 3.14159265358979323846;

    const = -1/(6*wL*1e6)*(3*Chi*1e6/(2*I*(2*I-1))).^2 .* (I*(I+1) - 3/4) ;
    %const = -2186; %3386.858187815424 - *около*точное значение
    dw = (w_max - w_min)/N;
    for p = 0:M
        for q = 0:M
            mu = p/M;
            phi = 2*my_pi*q/M;
            y = calcY(phi,mu,eta);
            w = const*y;
            w_ppm = w/wL;
            %index = calcIndex(w, w_min, w_max, N);
            index = floor((w_ppm - w_min)/dw);
            if (index <= 0)
                index
            end
            G_w(index,2) = G_w(index,2) + 1;  % добавляем точку в спектр
        end
    end
    
    %hist(a)
%   stop
%    figure
% %    plot(G_w(:,1), G_w(:,2));
    if (j_coupling == 1)
        J = J/wL;      % перевод в ppm
        G_w1 = G_w;
        G_w1(:,2) = G_w1(:,2)./2;
        G_w1 = shiftSpectrum(G_w1, -J/2);      % расщепляем в обе стороны
        G_w2 = shiftSpectrum(G_w, +J/2);
        G_w(:,2) = G_w1(:,2) + G_w2(:,2);
    end
    G_w = convolution(G_w, sigma, gamma, alpha);
    G_w = normalizeSpectrum(G_w, interval);
    
    %m = medium(G_w, interval);
    G_w = shiftSpectrum(G_w, delta);  
    %m = medium(G_w, interval)
end

