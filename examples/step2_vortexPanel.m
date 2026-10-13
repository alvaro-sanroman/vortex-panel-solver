%% STEP 2 - Vortex panel method: first-order vs second-order
%  First-order  : constant gamma per panel (Anderson 4.10)
%                 -> panel-to-panel oscillations in gamma and Cp
%  Second-order : linear gamma per panel (Anderson Fig. 4.34, Katz & Plotkin 11.4.2)
%                 -> smooth Cp, Cl from circulation and from pressure agree

clear; clc; close all

%% Paths
root   = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root, 'src')));
figDir = fullfile(root, 'figures');

%% Inputs
code     = '2412';
alphaDeg = 4;
nPanels  = 100;                       % for the Cp comparison
Nlist    = [20 40 80 160 320];        % for the convergence study

%% First-order vs second-order, same airfoil and panels
[x, z] = naca4(code, nPanels);
VP1 = vortexPanelConstant(x, z, alphaDeg);    % first-order
VP2 = vortexPanel(x, z, alphaDeg);            % second-order

fprintf('\nNACA %s, alpha = %g deg, %d panels\n', code, alphaDeg, nPanels);
fprintf('  first-order : Cl = %.4f   Cp from %7.2f to %.2f\n', VP1.Cl, min(VP1.Cp), max(VP1.Cp));
fprintf('  second-order: Cl = %.4f   Cp from %7.2f to %.2f\n', VP2.Cl, min(VP2.Cp), max(VP2.Cp));

%% Convergence of the second-order method
nN    = numel(Nlist);
ClKJ  = zeros(1, nN);
ClP   = zeros(1, nN);
CdP   = zeros(1, nN);
for k = 1:nN
    [xk, zk] = naca4(code, Nlist(k));
    VP = vortexPanel(xk, zk, alphaDeg);
    ClKJ(k) = VP.Cl;
    ClP(k)  = VP.ClP;
    CdP(k)  = VP.CdP;
end

fprintf('\nSecond-order convergence\n');
fprintf('  panels   Cl (KJ)   Cl (pressure)   Cd (pressure)\n');
for k = 1:nN
    fprintf('  %5d    %.4f      %.4f         %.5f\n', Nlist(k), ClKJ(k), ClP(k), CdP(k));
end

%% Figure
T   = plotTheme();
fig = figure('Color', T.bg, 'Position', [100 100 1400 450]);
tl  = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% nodes are numbered TE -> lower -> LE -> upper -> TE, so the first half
% of the panels is the lower surface and the second half the upper surface
lower = 1:nPanels/2;
upper = nPanels/2+1:nPanels;

% ---- (a) first-order Cp
ax1 = nexttile(tl);  hold(ax1, 'on');
plot(ax1, VP1.P.xm, VP1.Cp, '-o', 'Color', T.accent3, 'LineWidth', 1, 'MarkerSize', 3);
set(ax1, 'YDir', 'reverse');
xlabel(ax1, 'x/c');  ylabel(ax1, 'C_p');
title(ax1, sprintf('First-order (constant \\gamma):  C_l = %.3f', VP1.Cl));
styleAxes(ax1, T);

% ---- (b) second-order Cp
ax2 = nexttile(tl);  hold(ax2, 'on');
plot(ax2, VP2.P.xm(upper), VP2.Cp(upper), '-', 'Color', T.accent1, 'LineWidth', 1.8);
plot(ax2, VP2.P.xm(lower), VP2.Cp(lower), '-', 'Color', T.accent2, 'LineWidth', 1.8);
set(ax2, 'YDir', 'reverse');
xlabel(ax2, 'x/c');  ylabel(ax2, 'C_p');
title(ax2, sprintf('Second-order (linear \\gamma):  C_l = %.3f', VP2.Cl));
lg = legend(ax2, {'upper surface', 'lower surface'}, 'Location', 'southeast');
set(lg, 'TextColor', T.fg, 'Color', T.bg, 'EdgeColor', T.grid);
styleAxes(ax2, T);

% ---- (c) convergence
ax3 = nexttile(tl);  hold(ax3, 'on');
plot(ax3, Nlist, ClKJ, '-o', 'Color', T.accent1, 'LineWidth', 1.6, 'MarkerFaceColor', T.accent1);
plot(ax3, Nlist, ClP,  '-s', 'Color', T.accent2, 'LineWidth', 1.6, 'MarkerFaceColor', T.accent2);
set(ax3, 'XScale', 'log');
xlabel(ax3, 'number of panels');  ylabel(ax3, 'C_l');
title(ax3, 'Second-order: convergence');
lg = legend(ax3, {'Kutta-Joukowski', 'pressure integration'}, 'Location', 'southeast');
set(lg, 'TextColor', T.fg, 'Color', T.bg, 'EdgeColor', T.grid);
styleAxes(ax3, T);

title(tl, sprintf('Vortex panel method, NACA %s, \\alpha = %g°, %d panels', ...
      code, alphaDeg, nPanels), 'Color', T.fg, 'FontSize', 14);

%% Save
if ~exist(figDir, 'dir'), mkdir(figDir); end
outFile = fullfile(figDir, 'step2_vortexPanel.png');
exportgraphics(fig, outFile, 'Resolution', 200, 'BackgroundColor', 'current');
fprintf('\nFigure saved to %s\n', outFile);