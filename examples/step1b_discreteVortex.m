%% STEP 1b - Discrete vortex method
%  Solves thin airfoil theory numerically with N point vortices
%  (vortex at 1/4, control point at 3/4 of each segment) and shows
%  convergence to the analytical thin airfoil results.

clear; clc; close all

%% Paths
root   = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root, 'src')));
figDir = fullfile(root, 'figures');

%% Inputs
codes    = {'0012', '2412', '4412'};
alphaDeg = 4;
Nlist    = [1 2 3 5 8 10 15 20 30 50 100];   % numbers of vortices to test
Ngamma   = 40;                               % N used for the gamma(x) plot

%% Preallocate
nAirfoils = numel(codes);
nN        = numel(Nlist);
ClDV      = zeros(nAirfoils, nN);    % one row per airfoil, one column per N
CmDV      = zeros(nAirfoils, nN);
ClTA      = zeros(nAirfoils, 1);     % thin airfoil reference values
CmTA      = zeros(nAirfoils, 1);
xGamma        = cell(nAirfoils, 1);      % gamma(x) curves for the plot
gammaThin      = cell(nAirfoils, 1);
xVortex      = cell(nAirfoils, 1);
gammaDiscrete       = cell(nAirfoils, 1);

%% Loop over airfoils
nF  = 40;                                    % Fourier terms for gamma(theta)
thetaGamma = linspace(0.02, pi, 500)';              % theta > 0: gamma -> infinity at the LE

for i = 1:nAirfoils
    % camber line and analytical reference
    [~, ~, xc, ~, dzc, theta] = naca4(codes{i}, 2000);
    TA = thinAirfoil(theta, dzc, alphaDeg);
    ClTA(i) = TA.Cl;
    CmTA(i) = TA.Cm_c4;

    % discrete vortex method for increasing N
    for k = 1:nN
        DV = discreteVortex(xc, dzc, alphaDeg, Nlist(k));
        ClDV(i, k) = DV.Cl;
        CmDV(i, k) = DV.Cm_c4;
    end

    % gamma(theta) from thin airfoil theory, V_inf = 1
    An = zeros(nF, 1);
    for n = 1:nF
        An(n) = 2/pi * trapz(theta, dzc .* cos(n*theta));
    end
    xGamma{i}   = (1 - cos(thetaGamma)) / 2;
    gammaThin{i} = 2 * ( TA.A0 * (1 + cos(thetaGamma)) ./ sin(thetaGamma) + sin(thetaGamma * (1:nF)) * An );

    % gamma from discrete vortices: Gamma_j / segment length
    DV = discreteVortex(xc, dzc, alphaDeg, Ngamma);
    xVortex{i} = DV.xv;
    gammaDiscrete{i}  = DV.Gamma * Ngamma;
end

%% Convergence table
for i = 1:nAirfoils
    fprintf('\nNACA %s at alpha = %g deg  |  thin airfoil theory: Cl = %.4f, Cm_c4 = %.4f\n', ...
            codes{i}, alphaDeg, ClTA(i), CmTA(i));
    fprintf('    N        Cl     error Cl     Cm_c4    |error Cm|\n');
    for k = 1:nN
        fprintf('  %3d    %.4f    %6.2f %%    %7.4f     %.4f\n', Nlist(k), ...
                ClDV(i, k), 100*(ClDV(i, k) - ClTA(i))/ClTA(i), ...
                CmDV(i, k), abs(CmDV(i, k) - CmTA(i)));
    end
end

%% Figure
T      = plotTheme();
colors = {T.accent1, T.accent2, T.accent3};
labels = strcat({'NACA '}, codes);

fig = figure('Color', T.bg, 'Position', [100 100 1400 450]);
tl  = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% ---- (a) Cl convergence
ax1 = nexttile(tl);  hold(ax1, 'on');
h = gobjects(nAirfoils, 1);
for i = 1:nAirfoils
    yline(ax1, ClTA(i), '--', 'Color', colors{i}, 'LineWidth', 1.2, 'HandleVisibility', 'off');
    h(i) = plot(ax1, Nlist, ClDV(i, :), '-o', 'Color', colors{i}, 'LineWidth', 1.5, ...
                'MarkerFaceColor', colors{i}, 'MarkerSize', 5);
end
set(ax1, 'XScale', 'log');
xlabel(ax1, 'number of vortices N');  ylabel(ax1, 'C_l');
title(ax1, 'Lift coefficient  (dashed: thin airfoil theory)');
lg = legend(ax1, h, labels, 'Location', 'east');
set(lg, 'TextColor', T.fg, 'Color', T.bg, 'EdgeColor', T.grid);
styleAxes(ax1, T);

% ---- (b) Cm convergence
ax2 = nexttile(tl);  hold(ax2, 'on');
for i = 1:nAirfoils
    yline(ax2, CmTA(i), '--', 'Color', colors{i}, 'LineWidth', 1.2);
    plot(ax2, Nlist, CmDV(i, :), '-o', 'Color', colors{i}, 'LineWidth', 1.5, ...
         'MarkerFaceColor', colors{i}, 'MarkerSize', 5);
end
set(ax2, 'XScale', 'log');
xlabel(ax2, 'number of vortices N');  ylabel(ax2, 'C_{m,c/4}');
title(ax2, 'Moment about the quarter chord');
styleAxes(ax2, T);

% ---- (c) vortex-sheet strength
ax3 = nexttile(tl);  hold(ax3, 'on');
for i = 1:nAirfoils
    plot(ax3, xGamma{i}, gammaThin{i}, '-', 'Color', colors{i}, 'LineWidth', 1.6);
    plot(ax3, xVortex{i}, gammaDiscrete{i}, 'o', 'Color', colors{i}, 'MarkerSize', 4, ...
         'HandleVisibility', 'off');
end
ylim(ax3, [0 2.5]);
xlabel(ax3, 'x/c');  ylabel(ax3, '\gamma / V_\infty');
title(ax3, sprintf('Vortex sheet strength  (lines: theory, dots: N = %d)', Ngamma));
lg = legend(ax3, labels);
set(lg, 'TextColor', T.fg, 'Color', T.bg, 'EdgeColor', T.grid);
styleAxes(ax3, T);

title(tl, sprintf('Discrete vortex method vs thin airfoil theory,  \\alpha = %g°', alphaDeg), ...
      'Color', T.fg, 'FontSize', 14);

%% Save
if ~exist(figDir, 'dir'), mkdir(figDir); end
outFile = fullfile(figDir, 'step1b_discreteVortex.png');
exportgraphics(fig, outFile, 'Resolution', 200, 'BackgroundColor', 'current');
fprintf('\nFigure saved to %s\n', outFile);