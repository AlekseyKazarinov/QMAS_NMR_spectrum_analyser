function [index] = calcIndex(w, w_min, w_max, N)
%getIndex - ���������� ��������������� ������ � ������� G_w
% �������� y ����� �������� � [Ymin, Ymax]
% Ymin, Ymax �������������� ������� ��� ����������� ������
% N - ������ �������-�����������
dw = (w_max - w_min)/N;
w_n = w - w_min;
n = floor(w_n/dw);
index = n;
    %Yn = y - Ymin;
    %n = round(Yn/dY);
    %index = n + 1;
end

