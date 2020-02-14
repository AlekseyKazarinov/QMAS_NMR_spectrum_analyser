function [eps] = calcDiscrepancyWeighted(calcSpectrum, expSpectrum, interval)
% returns an weighted error of calculation between theoretical spectrum and
% experimental one using LS method (Least Squares).
% The weigth is proportional to specified value. That value may be distance
% between an frequency point and a median point of the whole spectrum
% or in another hand it can be proportional to absolute height of spectrum
% line

% Input arguments:
% 1) calcSpectrum - a theoretically expected spectrum line
% 2) expSpectrum - a specrum obtained from experimental data
% 3) interval - the interval on which calculation is performed


    left = interval(1);
    right = interval(2);

    variant = 2;
    % Possible values of the parameter are listed below:
    % 0 - standard deviation
    % 1 - user (custom) variant
    % 2 - the case of variant is  directly proportional to intencity of experimental spectrum

    if (variant == 0)
        eps = sqrt(sum((calcSpectrum(left:right,2) - expSpectrum(left:right,2)).^2));
    %     plot((calcSpectrum(left:right,2) - expSpectrum(left:right,2)).^2);
    %     stop
    end

    % The weight is proportional to a distance between an frequency point 
    % and a median point of the whole spectrum:
    if (variant ==1)
        % m - median value of a spectrum (ppm)
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

    % weight is proportional to absolute height of spectrum:
    if (variant == 2)
        eps = sqrt(sum( ((calcSpectrum(left:right,2) - expSpectrum(left:right,2)).^2) .* abs(expSpectrum(left:right,2)) ));
        %eps = sum(abs( (calcSpectrum(left:right,2) - expSpectrum(left:right,2)).*  abs(maxIntensity./calcSpectrum(left:right,2)) ));
        % visual representation of weight for calculating of error depending on different counts
        %plot( ((calcSpectrum(left:right,2) - expSpectrum(left:right,2)).^2) .* expSpectrum(left:right,2));
        %stop
    end
end


