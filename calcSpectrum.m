function [G_ws, G_w] = calcSpectrum(N, sp_mesh, M, interval, params, ratios, j_coupling, J)
%UNTITLED ������ ���������� ������� �� ����� � ����������� � parameters,
%   �� ������������ ������������� ������������ �������������� �� ratios,
%   ��� r1 + r2 + .. + rk + .. + rn = 1 , rk - ������� �� ratios

    % �������� �� ������ ����������:
    if (sum(ratios) < 0.999)
        stop       % ������ ����������
    end
    
    n_lines = length(ratios);  % ���������� ����� �����
    G_ws = zeros(N, 2, n_lines);
        
    % ������ ������ �����
    G_w1 = calcLine(N,sp_mesh,M,interval, params(1,:), j_coupling, J);
    G_w1(:,2) = ratios(1)*G_w1(:,2);
    G_w = G_w1;
    G_ws(:,:,1) = G_w1;
    % ��������� ��� ��������� ����� c ������ ����� ���������� �� 1
    for k = 2:n_lines
        G_wk = calcLine(N,sp_mesh,M,interval, params(k,:), 0, J);
        G_wk(:,2) = ratios(k)*G_wk(:,2);
        G_w(:,2) = G_w(:,2) + G_wk(:,2);
        G_ws(:,:,k) = G_wk;
    end
end

