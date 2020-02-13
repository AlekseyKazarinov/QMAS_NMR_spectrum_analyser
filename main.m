clear;
%clc;
% version 1.0.2
% ��������� ��� ������������� ��������������� MAS ������� �� ������� 23Na (I=3/2)
% ������������ ��������� �������.


%% ��������� ���������/��������� - ���������� ��� ������ ��������
file_name = 'id9117_93340_NaBiO3_23Na.023.001.1r_mod.txt'; % ����������������� ������
checkout = 0;    % ���/���� (1/0) ����� ������ �������

% ���������� ����� ������������� ��������������� ����� (� ����� = 1):
ratios = [0.5	0.5];

% ��������� �����:                        
%   eta   - �������� ���������� ������� ���
%   Chi   - �������������� ��������� (MHz)
%   delta - ����� �������, delta iso (ppm)
%   sigma - ������ (����.)������. ����� �� ���������� ������������� (Hz)
%   gamma - �������� �������� ���������� ����� (Hz)
%   alpha - ����������� ����� ��������� � ���������� ������ �����
%           0 - ����. ������ ���������� �����, 1 - ���������

% ��� ������� � ������ ��������� ��������� ��������� ����� �������� �� 
% �������. ������� ���������� ���������� ��� �����:
% [eta, Chi, delta, sigma, gamma, alpha]

params1 = [0    1.300   11.0   80.0000   75.000        0];
params2 = [0    1.5     12.8   80.0000   90.000        0];

params = [params1; params2];


% �������������� �������������� ��� 1-�� �����:
j_coupling = 0;         % J-coupling � 1H (0 - ���, 1 - ����)
J = 50;                 % �������� ����������� (Hz)


% ����� ������� ����� ���������:
common = 0; % (1/0) - �������� ����� ����� ����������
sigma = 80;
gamma = 70;
alpha = 0.5;

if (common == 1)
    for n_line = 1:length(ratios)
       params(n_line, 4) = sigma; 
       params(n_line, 5) = gamma;
       params(n_line, 6) = alpha;
    end
end
                        
%% ����������� ���������
epsilon = 1e-5;         % ����������� ��������
max_iter = 50;          % ����������� ����� ����������
M = 1000;               % ���������� ����� � �����, ����������� ��������� �����
left_bound = -25;       % ������ �������� �������������
right_bound = 30;


%% ���������������� ����

% ������ � ���������� ��������� ������������������ ������� �� �����
fileID = fopen(file_name, 'r');
formatSpec = '%f';
sizeA = [2 Inf];  % �.�. 2 �������, � ���-�� �������� ������� ����������
A = fscanf(fileID, formatSpec, sizeA); % 1 - ��� �����, 2 - �������������
A = A';
A = A(end:-1:1,:);  % ��������������, �.�. ������� "+" -> "-"
delta_max = A(length(A),1);
%plot(A(:,1), A(:,2))
interval = [left_bound, right_bound];
% ����������� �������� ��������� �������� �������
[l,r] = convertToIndices(interval, A);
interval = [l, r];
A = normalizeSpectrum(A, interval);
% ����� ������� ��������� ���. ������ �������:
delta_min = A(1,1);           
N = length(A(:,1));  % ���������� ����� � ����� ��������������� �������
           % ���������� M>N (�� ���� �������)

           
%% ���ר� �������� �������     
% ������ ������ 
% 1 - �������� �����, 2 - ��������������
[G_ws, G_w] = calcSpectrum(N, A(:,1), M, interval, params, ratios, j_coupling, J);
eps = calcDiscrepancyWeighted(G_w, A, interval);

% ������ ������
figure
hold on
plot(G_ws(:,1,1), G_ws(:,2,1), 'c--', G_ws(:,1,2), G_ws(:,2,2), 'c-.', 'LineWidth', 1.5)
for n_line = 3:length(ratios)
    plot(G_ws(:,1,n_line), G_ws(:,2,n_line), 'b:', 'LineWidth', 2)
end

plot(G_w(:,1), G_w(:,2), 'b', A(:,1), A(:,2), 'r', 'LineWidth', 1.5)
plot(G_w(:,1), A(:,2)-G_w(:,2), 'g-', 'LineWidth', 1.5)

xlabel('chemical shift (ppm)');
ylabel('Intensity (relative units)');
xlim([left_bound, right_bound]);

if (length(ratios) == 3)
    legend('line1', 'line2', 'line3', 'calculated sum', 'experimental', 'difference');
else
    legend('line1', 'line2', 'calculated sum', 'experimental', 'difference');
end

%title('Test spectrum');

% ����� ���������� � �������:
intencity = max(max(G_w(:,2)), max(A(:,2)));
p = 'eta     Chi     delta     sigma     gamma';
text(left_bound+1,1.0*intencity, ['eps = ' num2str(eps)] );
text(left_bound+1,0.95*intencity, 'ratios:');
text(left_bound+1,0.9*intencity, num2str(ratios));
text(left_bound+1,0.8*intencity, p);
text(left_bound+1,0.7*intencity, num2str(params(1:length(ratios),1:6)));

hold off

if (checkout == 1)
    stop   % ������� ����� ������ �������
end



%% �������� �������
hEta = 1e-3;  % ���� �� ������� �� ����������
hChi = 1e-2;
hDelta = 1e-1;
hSigma = 1;
hGamma = 1;
hAlpha = 1e-2;
hRatio = 1e-3;

dEta = 1e-3;   % ���������� � ������� �����������
dChi = 1e-3;
dDelta = 1e-1;
dSigma = 1;
dGamma = 0.5;
dAlpha = 1e-2;
dRatio = 1e-4;


diffs = [dEta, dChi, dDelta, dSigma, dGamma, dAlpha];
increments = [hEta, hChi, hDelta, hSigma, hGamma, hAlpha];

figure
% ���������� ����� ������ ��� ��������
v = linspace(0.9,0,max_iter);
map = [v' v' v'];

for iter = 1:max_iter
    iter
    
    % ��������� Eta
%      [params, eps, G_w] = stepLineParameter(1, 1, dEta, hEta, N, A, M, interval, params, ratios, j_coupling, J);
%      [params, eps, G_w] = stepLineParameter(2, 1, dEta, hEta, N, A, M, interval, params, ratios, j_coupling, J);
%     [params, eps, G_w] = stepLineParameter(3, 1, dEta, hEta, N, A, M, interval, params, ratios, j_coupling, J);

    % ��������� Chi
    [params, eps, G_w] = stepLineParameter(1, 2, dChi, hChi, N, A, M, interval, params, ratios, j_coupling, J);
    [params, eps, G_w] = stepLineParameter(2, 2, dChi, hChi, N, A, M, interval, params, ratios, j_coupling, J);
%    [params, eps, G_w] = stepLineParameter(3, 2, dChi, hChi, N, A, M, interval, params, ratios, j_coupling, J);

    % ��������� delta (isotropic)
    [params, eps, G_w] = stepLineParameter(1, 3, dDelta, hDelta, N, A, M, interval, params, ratios, j_coupling, J);
    [params, eps, G_w] = stepLineParameter(2, 3, dDelta, hDelta, N, A, M, interval, params, ratios, j_coupling, J);
%    [params, eps, G_w] = stepLineParameter(3, 3, dDelta, hDelta, N, A, M, interval, params, ratios, j_coupling, J);

    % ��������� gamma
    if (common == 0)
        [params, eps, G_w] = stepLineParameter(1, 5, dGamma, hGamma, N, A, M, interval, params, ratios, j_coupling, J);
        [params, eps, G_w] = stepLineParameter(2, 5, dGamma, hGamma, N, A, M, interval, params, ratios, j_coupling, J);
%        [params, eps, G_w] = stepLineParameter(3, 5, dGamma, hGamma, N, A, M, interval, params, ratios, j_coupling, J);
    end

    % ����� ����� ���������� ��� ���� ����� - sigma, gamma, alpha
%     [params, eps, G_w] = stepCommonParameter(4, dSigma, hSigma, N, A, M, interval, params, ratios, j_coupling, J);
    if (common == 1)
        [params, eps, G_w] = stepCommonParameter(5, dGamma, hGamma, N, A, M, interval, params, ratios, j_coupling, J);
    end       
%     [params, eps, G_w] = stepCommonParameter(6, dAlpha, hAlpha, N, A, M, interval, params, ratios, j_coupling, J);
    
    % �������� ������������� ������������ ��������������
      [ratios, eps, G_w] = StepRatios(dRatio, hRatio, N, A, M, interval, params, ratios, j_coupling, J);

    % ������� ���� � �������:
    eps
    ratios
    params(1:length(ratios),:)
    
    
    % ������ ������:
    hold on
    COLOR = [0, 0, 0.5+ 0.5*iter/max_iter];
    plot(G_w(:,1), G_w(:,2), 'color', map(iter,: ));  
    plot(A(:,1), A(:,2), 'r');
    plot(G_w(:,1), A(:,2)-G_w(:,2), 'g-')
    xlabel('chemical shift (ppm)');
    ylabel('Intensity (relative units)');
    xlim([left_bound, right_bound]);
    hold off
    legend('calculated', 'experimental', 'difference');
    if (eps <= epsilon)
        break
    end
end

