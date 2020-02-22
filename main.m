clear;
% clc;
% version 1.0.3
% This program is aimed to approximation of 2-nd order quadrupolar MAS NMR spectrum
% (for an example of 23Na I=3/2, but it may customize)
% By the result of fitting spectra main physical parameters of spectrum 
% will be calculated.


%% INITIAL PARAMETERS/CONSTANTS - are necessary for start of fitting process
% Name of file containing experimental spectrum (by default this file 
% expected to be in the same directory from which the program launch)
    file_name = 'id9117_93340_NaBiO3_23Na.023.001.1r_mod.txt'; 
    checkout = 1;    % on/off (1/0) mode of manual debugging (or viewing)

% Relation between integral intensities of spectum lines (sum of them = 1):
	ratios = [0.5	0.5];

% Parametes that forming spectrum lines:                        
%   eta   -asymmetry parameter of EFG tenzor
%   Chi   - quadrupolar constant (MHz)
%   delta - isotropic shift of spectrum, delta iso (ppm)
%   sigma - width (or dispersion) of spectral line from each monocryst отдельного монокристалла (Hz)
%   gamma - scale paramter of Lorenzian (Hz)
%   alpha - ratio coefficient between Gaussian and Lorencian types of lines
%           0 - Lorenzian only, 1 - Gaussian only

% All the functions of that program inputs paremeters of spectrum lines
% according theirs positions in params array.
% Here is the order of spectrum parameters:
% [eta, Chi, delta, sigma, gamma, alpha]

    params1 = [0    1.300   11.0   80.0000   75.000        0]; % parameters of the first component (line)
    params2 = [0    1.5     12.8   80.0000   90.000        0]; % 

    params = [params1; params2];


% Additional characteristics for 1-st line:
j_coupling = 0;         % J-coupling of 1H (0 - not existing, 1 - exists)
J = 50;                 % value of splitting (Hz)


% You may define that some of parameters are shared by all of lines
common = 0; % (1/0) - enable mode of shared parameters
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

%% TECHNICAL CONSTANTS
epsilon = 1e-5;         % error of fitting (approximation)
max_iter = 50;          % restriction of iteration number
M = 1000;               % number of points in mesh (or grid), placed around the sphere of radius = 1
left_bound = -25;       % bounds of approximation interval
right_bound = 30;


%% PREPARATION STAGE

% Reading and preparing original experimental spectrum from the file
fileID = fopen(file_name, 'r');
formatSpec = '%f';
sizeA = [2 Inf];  % because there are 2 columns, and number of values could not know apriori
A = fscanf(fileID, formatSpec, sizeA); % 1 - chem. shift, 2 - intensity
A = A';
A = A(end:-1:1,:);  % need to flip, as written "+" -> "-"
delta_max = A(length(A),1);
%plot(A(:,1), A(:,2))
interval = [left_bound, right_bound];
% define indices of interval on which spectrum to fit
[l,r] = convertToIndices(interval, A);
interval = [l, r];
A = normalizeSpectrum(A, interval);
% define bounds of chemical shift of spectrum:
delta_min = A(1,1);
N = length(A(:,1));  % number of points in mesh (grid) of calculating spectrum
           % optimal value M>N (by 10x)

           
%% CALCULATION OF SAMPLE SPECTRUM     
% Spectrum building 
% 1 - main line, 2 (3, etc) - additional one(s)
[G_ws, G_w] = calcSpectrum(N, A(:,1), M, interval, params, ratios, j_coupling, J);
eps = calcDiscrepancyWeighted(G_w, A, interval);

% plot the figure
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

% Output information about spectrum (character parameters):
intencity = max(max(G_w(:,2)), max(A(:,2)));
p = 'eta     Chi     delta     sigma     gamma';
text(left_bound+1,1.0*intencity, ['eps = ' num2str(eps)] );
text(left_bound+1,0.95*intencity, 'ratios:');
text(left_bound+1,0.9*intencity, num2str(ratios));
text(left_bound+1,0.8*intencity, p);
text(left_bound+1,0.7*intencity, num2str(params(1:length(ratios),1:6)));

hold off

if (checkout == 1)
    stop   % MODE OF MANUAL DEBUGGING (or VIEWING) IS ENABLED
end



%% SPECTRUM APPROXIMATION (fitting)
hEta = 1e-3;  % values of steps by each parameters
hChi = 1e-2;
hDelta = 1e-1;
hSigma = 1;
hGamma = 1;
hAlpha = 1e-2;
hRatio = 1e-3;

dEta = 1e-3;   % increases in partial derivatives
dChi = 1e-3;
dDelta = 1e-1;
dSigma = 1;
dGamma = 0.5;
dAlpha = 1e-2;
dRatio = 1e-4;


diffs = [dEta, dChi, dDelta, dSigma, dGamma, dAlpha];
increments = [hEta, hChi, hDelta, hSigma, hGamma, hAlpha];

figure
% preparation of colour map for figures
v = linspace(0.9,0,max_iter);
map = [v' v' v'];

for iter = 1:max_iter
    iter
    
    % variation Eta
%      [params, eps, G_w] = stepLineParameter(1, 1, dEta, hEta, N, A, M, interval, params, ratios, j_coupling, J);
%      [params, eps, G_w] = stepLineParameter(2, 1, dEta, hEta, N, A, M, interval, params, ratios, j_coupling, J);
%     [params, eps, G_w] = stepLineParameter(3, 1, dEta, hEta, N, A, M, interval, params, ratios, j_coupling, J);

    % variation Chi
    [params, eps, G_w] = stepLineParameter(1, 2, dChi, hChi, N, A, M, interval, params, ratios, j_coupling, J);
    [params, eps, G_w] = stepLineParameter(2, 2, dChi, hChi, N, A, M, interval, params, ratios, j_coupling, J);
%    [params, eps, G_w] = stepLineParameter(3, 2, dChi, hChi, N, A, M, interval, params, ratios, j_coupling, J);

    % variation delta (isotropic)
    [params, eps, G_w] = stepLineParameter(1, 3, dDelta, hDelta, N, A, M, interval, params, ratios, j_coupling, J);
    [params, eps, G_w] = stepLineParameter(2, 3, dDelta, hDelta, N, A, M, interval, params, ratios, j_coupling, J);
%    [params, eps, G_w] = stepLineParameter(3, 3, dDelta, hDelta, N, A, M, interval, params, ratios, j_coupling, J);

    % variation gamma
    if (common == 0)
        [params, eps, G_w] = stepLineParameter(1, 5, dGamma, hGamma, N, A, M, interval, params, ratios, j_coupling, J);
        [params, eps, G_w] = stepLineParameter(2, 5, dGamma, hGamma, N, A, M, interval, params, ratios, j_coupling, J);
%        [params, eps, G_w] = stepLineParameter(3, 5, dGamma, hGamma, N, A, M, interval, params, ratios, j_coupling, J);
    end

    % if shared parameters mode is enabled (sigma, gamma, alpha are the same for each spectrum line)
%     [params, eps, G_w] = stepCommonParameter(4, dSigma, hSigma, N, A, M, interval, params, ratios, j_coupling, J);
    if (common == 1)
        [params, eps, G_w] = stepCommonParameter(5, dGamma, hGamma, N, A, M, interval, params, ratios, j_coupling, J);
    end       
%     [params, eps, G_w] = stepCommonParameter(6, dAlpha, hAlpha, N, A, M, interval, params, ratios, j_coupling, J);
    
    % variation of relative intergal intensities
      [ratios, eps, G_w] = StepRatios(dRatio, hRatio, N, A, M, interval, params, ratios, j_coupling, J);

    % logs are output into Matlab console:
    eps
    ratios
    params(1:length(ratios),:)
    
    
    % plot the figure:
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

