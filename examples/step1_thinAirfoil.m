%% STEP 1 - Thin airfoil theory
%  Applies classical thin airfoil theory (Anderson 4.7-4.9) to several
%  NACA 4-digit airfoils and compares Cl-alpha, alpha_L0, Cm_c/4 and x_cp.

clear; clc; close all

%% Paths
root   = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root, 'src')));
figDir = fullfile(root, 'figures');

%% Inputs
codes    = {'0012', '2412', '4412'};   % airfoils to compare
alphaDeg = -6:0.5:12;                  % angle-of-attack sweep [deg]
nPanels  = 2000;                       % fine grid -> accurate integrals

%% Thin airfoil theory for each airfoil
nAirfoils  = numel(codes);
nAlpha     = numel(alphaDeg);
Cl         = zeros(nAirfoils, nAlpha);   % one row per airfoil, one column per alpha
xcp        = zeros(nAirfoils, nAlpha);
alphaL0    = zeros(nAirfoils, 1);        % zero-lift angle [deg]
alphaIdeal = zeros(nAirfoils, 1);        % ideal (design) angle [deg]
Cm_c4      = zeros(nAirfoils, 1);        % moment about c/4 (= about the a.c.)
xcAll      = cell(nAirfoils, 1);
zcAll      = cell(nAirfoils, 1);

for k = 1:nAirfoils
    [~, ~, xc, zc, dzc, theta] = naca4(codes{k}, nPanels);
    TA = thinAirfoil(theta, dzc, alphaDeg);

    Cl(k, :)      = TA.Cl;
    xcp(k, :)     = TA.xcp;
    alphaL0(k)    = TA.alphaL0;
    alphaIdeal(k) = TA.alphaIdeal;
    Cm_c4(k)      = TA.Cm_c4;
    xcAll{k}      = xc;
    zcAll{k}      = zc;
end

%% Results table
fprintf('\n  Airfoil    alpha_L0 [deg]   alpha_ideal [deg]    Cm_c/4    dCl/dalpha [1/rad]\n');
for k = 1:nAirfoils
    slope = (Cl(k, end) - Cl(k, 1)) / deg2rad(alphaDeg(end) - alphaDeg(1));
    fprintf('  NACA %s    %8.3f          %8.3f         %8.4f      %8.4f\n', ...
            codes{k}, alphaL0(k), alphaIdeal(k), Cm_c4(k), slope);
end

%% Figure
T      = plotTheme();
colors = {T.accent1, T.accent2, T.accent3};
labels = strcat({'NACA '}, codes);

fig = figure('Color', T.bg, 'Position', [100 100 1400 450]);
tl  = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% ---- (a) camber lines
ax1 = nexttile(tl);  hold(ax1, 'on');
for k = 1:nAirfoils
    plot(ax1, xcAll{k}, zcAll{k}, 'Color', colors{k}, 'LineWidth', 1.8);
end
xlabel(ax1, 'x/c');  ylabel(ax1, 'z_c/c');
title(ax1, 'Mean camber lines');
styleAxes(ax1, T);

% ---- (b) lift curves
ax2 = nexttile(tl);  hold(ax2, 'on');
xline(ax2, 0, ':', 'Color', T.grid, 'HandleVisibility', 'off');
yline(ax2, 0, ':', 'Color', T.grid, 'HandleVisibility', 'off');
h = gobjects(nAirfoils, 1);
for k = 1:nAirfoils
    h(k) = plot(ax2, alphaDeg, Cl(k, :), 'Color', colors{k}, 'LineWidth', 1.8);
    plot(ax2, alphaL0(k), 0, 'o', 'Color', colors{k}, 'MarkerFaceColor', colors{k}, ...
         'HandleVisibility', 'off');       % zero-lift angle
end
xlabel(ax2, '\alpha [deg]');  ylabel(ax2, 'C_l');
title(ax2, 'Lift curve  (slope = 2\pi for all)');
lg = legend(ax2, h, labels, 'Location', 'northwest');
set(lg, 'TextColor', T.fg, 'Color', T.bg, 'EdgeColor', T.grid);
styleAxes(ax2, T);

% ---- (c) centre of pressure
ax3 = nexttile(tl);  hold(ax3, 'on');
xcpPlot = xcp;
xcpPlot(abs(Cl) < 0.1) = NaN;       % hide the jump to infinity where Cl -> 0
for k = 1:nAirfoils
    plot(ax3, alphaDeg, xcpPlot(k, :), 'Color', colors{k}, 'LineWidth', 1.8);
end
yline(ax3, 0.25, '--', 'Color', T.fg, 'HandleVisibility', 'off');
ylim(ax3, [0 1]);
xlabel(ax3, '\alpha [deg]');  ylabel(ax3, 'x_{cp}/c');
title(ax3, 'Centre of pressure  (a.c. at c/4)');
styleAxes(ax3, T);

%% Save
if ~exist(figDir, 'dir'), mkdir(figDir); end
outFile = fullfile(figDir, 'step1_thinAirfoil.png');
exportgraphics(fig, outFile, 'Resolution', 200, 'BackgroundColor', 'current');
fprintf('\nFigure saved to %s\n', outFile);