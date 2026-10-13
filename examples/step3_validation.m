%% STEP 3 - Comparison and validation
%  (A) Second-order panel method vs thin airfoil theory: Cl-alpha and Cm_c/4
%  (C) Validation against an exact solution: 15%-thick van de Vooren airfoil,
%      the test case of Katz & Plotkin Fig. 11.31

clear; clc; close all

%% Paths
root   = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root, 'src')));
figDir = fullfile(root, 'figures');

%% Inputs
codes    = {'0012', '2412', '4412'};
alphaDeg = -4:1:10;                    % angle-of-attack sweep [deg]
nPanels  = 200;

%% (A) Panel method vs thin airfoil theory
nAirfoils = numel(codes);
nAlpha    = numel(alphaDeg);
ClPanel   = zeros(nAirfoils, nAlpha);
CmPanel   = zeros(nAirfoils, nAlpha);
ClThin    = zeros(nAirfoils, nAlpha);
CmThin    = zeros(nAirfoils, 1);

for i = 1:nAirfoils
    [x, z]                   = naca4(codes{i}, nPanels);
    [~, ~, ~, ~, dzc, theta] = naca4(codes{i}, 2000);   % fine grid for thin airfoil theory
    TA = thinAirfoil(theta, dzc, alphaDeg);
    ClThin(i, :) = TA.Cl;
    CmThin(i)    = TA.Cm_c4;
    
    for k = 1:nAlpha
        VP = vortexPanel(x, z, alphaDeg(k));
        ClPanel(i, k) = VP.Cl;
        CmPanel(i, k) = VP.Cm_c4;
    end
end

% lift slope and zero-lift angle (straight-line fit Cl = slope*alpha + Cl0)
fprintf('\n  Airfoil    dCl/da panel   dCl/da thin   alpha_L0 panel   alpha_L0 thin\n');

for i = 1:nAirfoils
    pPanel = polyfit(deg2rad(alphaDeg), ClPanel(i, :), 1);
    pThin  = polyfit(deg2rad(alphaDeg), ClThin(i, :),  1);
    fprintf('  NACA %s    %6.3f         %6.3f          %6.3f           %6.3f\n', codes{i}, ...
            pPanel(1), pThin(1), -rad2deg(pPanel(2)/pPanel(1)), -rad2deg(pThin(2)/pThin(1)));
end

%% (C) Validation: van de Vooren airfoil, alpha = 5 deg (K&P Fig. 11.31)
epsVdV = 0.10;   kVdV = 1.97;   alphaVdV = 5;
[xv, zv]           = vanDeVooren(epsVdV, kVdV, 90, alphaVdV);     % 90 panels, as K&P
[xe, ze, CpE, ClE] = vanDeVooren(epsVdV, kVdV, 2000, alphaVdV);   % fine grid: exact curve
VPv = vortexPanel(xv, zv, alphaVdV);

fprintf('\nvan de Vooren (eps = %.2f, k = %.2f), alpha = %g deg, 90 panels\n', epsVdV, kVdV, alphaVdV);
fprintf('  Cl panel = %.4f   Cl exact = %.4f   error = %.2f %%\n', VPv.Cl, ClE, 100*(VPv.Cl - ClE)/ClE);

%% Figure
T      = plotTheme();
colors = {T.accent1, T.accent2, T.accent3};
labels = strcat({'NACA '}, codes);

fig = figure('Color', T.bg, 'Position', [100 100 1450 470]);
tl  = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% ---- (a) lift curves: panel (solid) vs thin airfoil theory (dashed)
ax1 = nexttile(tl);  hold(ax1, 'on');
h = gobjects(nAirfoils, 1);

for i = 1:nAirfoils
    h(i) = plot(ax1, alphaDeg, ClPanel(i, :), '-',  'Color', colors{i}, 'LineWidth', 1.8);
    plot(ax1, alphaDeg, ClThin(i, :), '--', 'Color', colors{i}, 'LineWidth', 1.2, ...
         'HandleVisibility', 'off');
end

xlabel(ax1, '\alpha [deg]');  ylabel(ax1, 'C_l');
title(ax1, 'Lift: panel (solid) vs thin airfoil (dashed)');
lg = legend(ax1, h, labels, 'Location', 'northwest');
set(lg, 'TextColor', T.fg, 'Color', T.bg, 'EdgeColor', T.grid);
styleAxes(ax1, T);

% ---- (b) moment about c/4
ax2 = nexttile(tl);  hold(ax2, 'on');
for i = 1:nAirfoils
    plot(ax2, alphaDeg, CmPanel(i, :), '-', 'Color', colors{i}, 'LineWidth', 1.8);
    yline(ax2, CmThin(i), '--', 'Color', colors{i}, 'LineWidth', 1.2);
end

xlabel(ax2, '\alpha [deg]');  ylabel(ax2, 'C_{m,c/4}');
title(ax2, 'Moment about c/4: panel vs thin airfoil');
styleAxes(ax2, T);

% ---- (c) validation against the exact solution (cf. K&P Fig. 11.31)
ax3 = nexttile(tl);  hold(ax3, 'on');
fill(ax3, xv, -zv, T.body, 'EdgeColor', T.grid, 'HandleVisibility', 'off');   % airfoil outline
plot(ax3, xe, CpE, '-', 'Color', T.fg, 'LineWidth', 1.4);
plot(ax3, VPv.P.xm, VPv.Cp, '^', 'Color', T.accent2, 'MarkerSize', 4, 'MarkerFaceColor', T.accent2);
set(ax3, 'YDir', 'reverse');
xlabel(ax3, 'x/c');  ylabel(ax3, 'C_p');

title(ax3, sprintf('van de Vooren 15%%, \\alpha = %g°:  C_l = %.4f (exact %.4f)', ...
      alphaVdV, VPv.Cl, ClE));

lg = legend(ax3, {'exact solution', 'panel method (90 panels)'}, 'Location', 'southeast');
set(lg, 'TextColor', T.fg, 'Color', T.bg, 'EdgeColor', T.grid);
styleAxes(ax3, T);

%% Save
if ~exist(figDir, 'dir'), mkdir(figDir); end

outFile = fullfile(figDir, 'step3_validation.png');
exportgraphics(fig, outFile, 'Resolution', 200, 'BackgroundColor', 'current');
fprintf('\nFigure saved to %s\n', outFile);