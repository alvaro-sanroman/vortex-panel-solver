%% STEP 0 - Airfoil geometry and panelling
%  Generates a NACA 4-digit airfoil with cosine spacing and plots it.
%  Run with F5: the paths are resolved from this file's location.

clear; clc; close all

%% Paths
root   = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root, 'src')));
figDir = fullfile(root, 'figures');

%% Inputs
code    = '2412';   % try '0012', '4412', '2612'...
nPanels = 60;       % few panels so the nodes are visible

%% Geometry
[x, z, xc, zc] = naca4(code, nPanels);

%% Sanity checks
fprintf('NACA %s with %d panels\n', code, nPanels);
fprintf('  nodes      : %d (should be %d)\n', numel(x), nPanels + 1);
fprintf('  TE gap     : %.2e\n', hypot(x(1) - x(end), z(1) - z(end)));
[zcMax, iMax] = max(zc);
fprintf('  max camber : %.4f at x/c = %.3f\n', zcMax, xc(iMax));

%% Figure
T   = plotTheme();
fig = figure('Color', T.bg, 'Position', [100 100 1200 480]);
tl  = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% ---- (a) full airfoil, spanning 2 of the 3 tiles
ax1 = nexttile(tl, [1 2]);  hold(ax1, 'on');
fill(ax1, x, z, T.body, 'EdgeColor', 'none', 'HandleVisibility', 'off');
hP = plot(ax1, x, z, '-o', 'Color', T.accent1, 'LineWidth', 1.2, ...
          'MarkerSize', 3.5, 'MarkerFaceColor', T.accent1);
hC = plot(ax1, xc, zc, '--', 'Color', T.accent2, 'LineWidth', 1.2);
axis(ax1, 'equal');  xlim(ax1, [-0.05 1.12]);
xlabel(ax1, 'x/c');  ylabel(ax1, 'z/c');
title(ax1, sprintf('NACA %s  -  %d panels, cosine spacing', code, nPanels));
lg = legend(ax1, [hP hC], {'panels / nodes', 'mean camber line'});
set(lg, 'TextColor', T.fg, 'Color', T.bg, 'EdgeColor', T.grid);
styleAxes(ax1, T);

% ---- (b) leading-edge zoom
ax2 = nexttile(tl);  hold(ax2, 'on');
fill(ax2, x, z, T.body, 'EdgeColor', 'none');
plot(ax2, x, z, '-o', 'Color', T.accent1, 'LineWidth', 1.2, ...
     'MarkerSize', 4, 'MarkerFaceColor', T.accent1);
plot(ax2, xc, zc, '--', 'Color', T.accent2, 'LineWidth', 1.2);
axis(ax2, 'equal');  xlim(ax2, [-0.01 0.10]);  ylim(ax2, [-0.05 0.07]);
xlabel(ax2, 'x/c');  ylabel(ax2, 'z/c');
title(ax2, 'Leading-edge zoom');
styleAxes(ax2, T);

%% Save
if ~exist(figDir, 'dir'), mkdir(figDir); end
outFile = fullfile(figDir, sprintf('step0_geometry_NACA%s.png', code));
exportgraphics(fig, outFile, 'Resolution', 200, 'BackgroundColor', 'current');
fprintf('Figure saved to %s\n', outFile)