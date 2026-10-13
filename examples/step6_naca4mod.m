%% STEP 6 - NACA modified 4-digit sections (F-86F "6-3" wing)
%  Generates NACA 0009-64 (root) and 0008-64 (tip) with naca4mod,
%  checks them against NASA TM X-3284, compares 0009-64 with a
%  standard 0009 and exports the coordinates.

clear; clc; close all
root    = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root, 'src')));
figDir  = fullfile(root, 'figures');
dataDir = fullfile(root, 'data');

if ~exist(dataDir, 'dir'), mkdir(dataDir); end

codes    = {'0009-64', '0008-64'};
RleTable = [0.008925, 0.007052];        % TM X-3284 Table I, index 6
alphaDeg = 0;

% ------------------------------------------------------------ checks

for k = 1:2
    [x, z] = naca4mod(codes{k}, 2000);
    nS = 1000;                          % LE node is nS + 1
    [zMax, i] = max(z);
    Rle = z(nS+2)^2 / (2*x(nS+2));      % near LE: z ~ sqrt(2*R*x)
    fprintf('NACA %s\n', codes{k});
    fprintf('  max thickness : %.4f at x/c = %.3f\n', 2*zMax, x(i));
    fprintf('  LE radius     : %.6f (Table I: %.6f)\n', Rle, RleTable(k));

    % export: columns x/c, z/c, 0 (chord = 1, scale by the real chord in CAD)
    [xs, zs] = naca4mod(codes{k}, 160);
    fname = fullfile(dataDir, sprintf('NACA%s.txt', codes{k}));
    writematrix([xs zs zeros(size(xs))], fname, 'Delimiter', 'tab');
    fprintf('  coordinates   : %s\n\n', fname);
end

% -------------------------------- compare 0009-64 with standard 0009

nP = 200;
[xA, zA] = naca4('0009', nP);
[xB, zB] = naca4mod('0009-64', nP);
VA = vortexPanel(xA, zA, alphaDeg);
VB = vortexPanel(xB, zB, alphaDeg);
fprintf('alpha = %g deg:  Cl(0009) = %.4f   Cl(0009-64) = %.4f\n', ...
        alphaDeg, VA.Cl, VB.Cl);

% minimum Cp (highest local velocity) and where it occurs

[cpA, iA] = min(VA.Cp);   [cpB, iB] = min(VB.Cp);
fprintf('Cp min: 0009 = %.3f at x/c = %.3f | 0009-64 = %.3f at x/c = %.3f\n', ...
        cpA, VA.P.xm(iA), cpB, VB.P.xm(iB));
% ------------------------------------------------------------- plot

T   = plotTheme();
fig = figure('Color', T.bg, 'Position', [100 100 1100 700], 'WindowStyle', 'normal');
tl  = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% (a) geometry

ax1 = nexttile(tl); hold(ax1, 'on');
plot(ax1, xA, zA, '-', 'Color', T.accent1, 'LineWidth', 1.6);
plot(ax1, xB, zB, '-', 'Color', T.accent2, 'LineWidth', 1.6);
xline(ax1, 0.3, ':', 'Color', T.accent1, 'HandleVisibility', 'off');
xline(ax1, 0.4, ':', 'Color', T.accent2, 'HandleVisibility', 'off');
xlabel(ax1, 'x/c'); ylabel(ax1, 'z/c');
title(ax1, 'Geometry (vertical scale exaggerated) - max thickness at 30 % vs 40 %');
lg1 = legend(ax1, {'NACA 0009', 'NACA 0009-64'}, 'Location', 'northeast');
set(lg1, 'TextColor', T.fg, 'Color', T.bg, 'EdgeColor', T.grid);
styleAxes(ax1, T);

% (b) pressure distribution

ax2 = nexttile(tl); hold(ax2, 'on');
plot(ax2, VA.P.xm, VA.Cp, '-', 'Color', T.accent1, 'LineWidth', 1.6);
plot(ax2, VB.P.xm, VB.Cp, '-', 'Color', T.accent2, 'LineWidth', 1.6);
set(ax2, 'YDir', 'reverse');            % negative Cp (suction) upwards
xlabel(ax2, 'x/c'); ylabel(ax2, 'C_p');
title(ax2, sprintf('Pressure distribution, \\alpha = %g°', alphaDeg));
lg2 = legend(ax2, {'NACA 0009', 'NACA 0009-64'}, 'Location', 'southeast');
set(lg2, 'TextColor', T.fg, 'Color', T.bg, 'EdgeColor', T.grid);
styleAxes(ax2, T);

ax1.Toolbar.Visible = 'off';  ax2.Toolbar.Visible = 'off';
exportgraphics(fig, fullfile(figDir, 'step6_naca4mod.png'), ...
               'Resolution', 200, 'BackgroundColor', 'current');