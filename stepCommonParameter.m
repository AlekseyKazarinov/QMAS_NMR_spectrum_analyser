function [parameters, eps, G_w] = stepCommonParameter(number, dValue, hValue, N, A, M, interval, parameters, ratios, j_coupling, J)
%stepCommonParameter - ����������� ����� �� ������ ��������� ���� �����
%   number_line - ����� ����� (������ 1 ��� 2)
%   number - ����� ��������� � ������ parameters �����. �����:
%           eta -> 1, Chi -> 2, delta -> 3, sigma -> 4, gamma -> 5
%   sp_mesh - ����� ��� ������� (������� ���. ������)
%   dValue, hValue - ���������� � ��� ���������
%   interval = [l,r] - �������� ������� ������� �������, ��� �����
%           ��������� ������
%   parameters - ������ ���������� ��� ���� ����� � �������
%   ratios - ����������� ����� ��������������� �����
    
    fast_mode = 1;  % (1/0) - ���/���� ����������� ����� (������� ������ �� ���� �����������)

    [G_ws, G_w] = calcSpectrum(N, A(:,1), M, interval, parameters, ratios, j_coupling, J);
    eps = calcDiscrepancyWeighted(G_w, A, interval);
    
    % ���������� �� ������ ��������� ���� �����
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
    
    % ������ ������� ���������� ������ ��� gamma
    if (parameters(1,number) - decrement > 0 && parameters(1,number) - decrement <=100)  
        parameters(:,number) = parameters(:,number) - decrement;
    end
end

