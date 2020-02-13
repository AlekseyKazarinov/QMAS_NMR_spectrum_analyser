clear;
%clc;
% version 1.0.2
% Программа для аппроксимации квадриупольного MAS спектра на примере 23Na (I=3/2)
% Определяются параметры спектра.


%% НАЧАЛЬНЫЕ ПАРАМЕТРЫ/КОНСТАНТЫ - необходимы для старта подгонки
file_name = 'id9117_93340_NaBiO3_23Na.023.001.1r_mod.txt'; % экспериментальный спектр
checkout = 0;    % вкл/выкл (1/0) режим ручной отладки

% сотношения между интегральными интенсивностями линий (в сумме = 1):
ratios = [0.5	0.5];

% Параметры линий:                        
%   eta   - параметр асимметрии тензора ГЭП
%   Chi   - квадриупольная константа (MHz)
%   delta - сдвиг спектра, delta iso (ppm)
%   sigma - ширина (дисп.)спектр. линии от отдельного монокристалла (Hz)
%   gamma - параметр масштаба лоренцевой линии (Hz)
%   alpha - соотношение между гауссовым и лоренцевым типами линий
%           0 - соот. только лоренцевой линии, 1 - гауссовой

% Все функции в данной программе принимают параметры линий согласно их 
% позиции. Порядок следования параметров для линий:
% [eta, Chi, delta, sigma, gamma, alpha]

params1 = [0    1.300   11.0   80.0000   75.000        0];
params2 = [0    1.5     12.8   80.0000   90.000        0];

params = [params1; params2];


% Дополнительные характеристики для 1-ой линии:
j_coupling = 0;         % J-coupling с 1H (0 - нет, 1 - есть)
J = 50;                 % величина расщепления (Hz)


% Можно указать общие параметры:
common = 0; % (1/0) - включает режим общих параметров
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
                        
%% ТЕХНИЧЕСКИЕ КОНСТАНТЫ
epsilon = 1e-5;         % погрешность подгонки
max_iter = 50;          % ограничение числа повторений
M = 1000;               % количество узлов в сетке, покрывающей единичную сферу
left_bound = -25;       % задают интервал аппроксимации
right_bound = 30;


%% ПОДГОТОВИТЕЛЬНЫЙ ЭТАП

% Чтение и подготовка исходного экспериментального спектра из файла
fileID = fopen(file_name, 'r');
formatSpec = '%f';
sizeA = [2 Inf];  % т.к. 2 столбца, а кол-во значений заранее неизвестно
A = fscanf(fileID, formatSpec, sizeA); % 1 - хим сдвиг, 2 - интенсивность
A = A';
A = A(end:-1:1,:);  % переворачиваем, т.к. записан "+" -> "-"
delta_max = A(length(A),1);
%plot(A(:,1), A(:,2))
interval = [left_bound, right_bound];
% определение индексов интервала подгонки спектра
[l,r] = convertToIndices(interval, A);
interval = [l, r];
A = normalizeSpectrum(A, interval);
% задаём границы диапазона хим. сдвига спектра:
delta_min = A(1,1);           
N = length(A(:,1));  % количество точек в сетке рассчитываемого спектра
           % оптимально M>N (на один порядок)

           
%% РАСЧЁТ ПРОБНОГО СПЕКТРА     
% строим спектр 
% 1 - основная линия, 2 - дополнительная
[G_ws, G_w] = calcSpectrum(N, A(:,1), M, interval, params, ratios, j_coupling, J);
eps = calcDiscrepancyWeighted(G_w, A, interval);

% рисуем график
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

% вывод информации о спектре:
intencity = max(max(G_w(:,2)), max(A(:,2)));
p = 'eta     Chi     delta     sigma     gamma';
text(left_bound+1,1.0*intencity, ['eps = ' num2str(eps)] );
text(left_bound+1,0.95*intencity, 'ratios:');
text(left_bound+1,0.9*intencity, num2str(ratios));
text(left_bound+1,0.8*intencity, p);
text(left_bound+1,0.7*intencity, num2str(params(1:length(ratios),1:6)));

hold off

if (checkout == 1)
    stop   % ВКЛЮЧЕН РЕЖИМ РУЧНОЙ ОТЛАДКИ
end



%% ПОДГОНКА СПЕКТРА
hEta = 1e-3;  % шаги по каждому из параметров
hChi = 1e-2;
hDelta = 1e-1;
hSigma = 1;
hGamma = 1;
hAlpha = 1e-2;
hRatio = 1e-3;

dEta = 1e-3;   % приращения в частных производных
dChi = 1e-3;
dDelta = 1e-1;
dSigma = 1;
dGamma = 0.5;
dAlpha = 1e-2;
dRatio = 1e-4;


diffs = [dEta, dChi, dDelta, dSigma, dGamma, dAlpha];
increments = [hEta, hChi, hDelta, hSigma, hGamma, hAlpha];

figure
% подготовка карты цветов для графиков
v = linspace(0.9,0,max_iter);
map = [v' v' v'];

for iter = 1:max_iter
    iter
    
    % варьируем Eta
%      [params, eps, G_w] = stepLineParameter(1, 1, dEta, hEta, N, A, M, interval, params, ratios, j_coupling, J);
%      [params, eps, G_w] = stepLineParameter(2, 1, dEta, hEta, N, A, M, interval, params, ratios, j_coupling, J);
%     [params, eps, G_w] = stepLineParameter(3, 1, dEta, hEta, N, A, M, interval, params, ratios, j_coupling, J);

    % варьируем Chi
    [params, eps, G_w] = stepLineParameter(1, 2, dChi, hChi, N, A, M, interval, params, ratios, j_coupling, J);
    [params, eps, G_w] = stepLineParameter(2, 2, dChi, hChi, N, A, M, interval, params, ratios, j_coupling, J);
%    [params, eps, G_w] = stepLineParameter(3, 2, dChi, hChi, N, A, M, interval, params, ratios, j_coupling, J);

    % варьируем delta (isotropic)
    [params, eps, G_w] = stepLineParameter(1, 3, dDelta, hDelta, N, A, M, interval, params, ratios, j_coupling, J);
    [params, eps, G_w] = stepLineParameter(2, 3, dDelta, hDelta, N, A, M, interval, params, ratios, j_coupling, J);
%    [params, eps, G_w] = stepLineParameter(3, 3, dDelta, hDelta, N, A, M, interval, params, ratios, j_coupling, J);

    % варьируем gamma
    if (common == 0)
        [params, eps, G_w] = stepLineParameter(1, 5, dGamma, hGamma, N, A, M, interval, params, ratios, j_coupling, J);
        [params, eps, G_w] = stepLineParameter(2, 5, dGamma, hGamma, N, A, M, interval, params, ratios, j_coupling, J);
%        [params, eps, G_w] = stepLineParameter(3, 5, dGamma, hGamma, N, A, M, interval, params, ratios, j_coupling, J);
    end

    % вдоль общих параметров для всех линий - sigma, gamma, alpha
%     [params, eps, G_w] = stepCommonParameter(4, dSigma, hSigma, N, A, M, interval, params, ratios, j_coupling, J);
    if (common == 1)
        [params, eps, G_w] = stepCommonParameter(5, dGamma, hGamma, N, A, M, interval, params, ratios, j_coupling, J);
    end       
%     [params, eps, G_w] = stepCommonParameter(6, dAlpha, hAlpha, N, A, M, interval, params, ratios, j_coupling, J);
    
    % вариация относительных интегральных интенсивностей
      [ratios, eps, G_w] = StepRatios(dRatio, hRatio, N, A, M, interval, params, ratios, j_coupling, J);

    % пишутся логи в консоль:
    eps
    ratios
    params(1:length(ratios),:)
    
    
    % рисуем график:
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

