%% STEP 4 - Flow field: streamlines and pressure field
%  Velocity anywhere = freestream + all vortex panels (superposition).
%  Check: the velocity inside the airfoil must be ~0 (Anderson Fig. 4.33).

clear; clc; close all

%% Paths

root   = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root, 'src')));
figDir = fullfile(root, 'figures');

%% Inputs

code     = '2412';
alphaDeg = 5;
nPanels  = 200;
xLim = [-0.5 1.5];   zLim = [-0.5 0.5];      % plotting window
nx   = 300;          nz   = 150;             % grid resolution

%% Solve the airfoil

[x, z] = naca4(code, nPanels);
VP = vortexPanel(x, z, alphaDeg);

%% Velocity on a grid

[X, Z] = meshgrid(linspace(xLim(1), xLim(2), nx), linspace(zLim(1), zLim(2), nz));
[U, W, inside] = flowField(VP, X, Z, alphaDeg);
Vmag = hypot(U, W);

fprintf('\nNACA %s, alpha = %g deg, Cl = %.4f\n', code, alphaDeg, VP.Cl);
fprintf('  mean speed INSIDE the airfoil = %.5f  (should be ~0)\n', mean(Vmag(inside)));

% hide the inside of the airfoil, then pressure from Bernoulli

U(inside) = NaN;   W(inside) = NaN;   Vmag(inside) = NaN;
CpField = 1 - Vmag.^2;

fprintf('  speed range around the airfoil: %.3f to %.3f\n', min(Vmag(:)), max(Vmag(:)));

%% Figure

T   = plotTheme();
fig = figure('Color', T.bg, 'Position', [100 50 1100 950], 'WindowStyle', 'normal');
tl  = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% ---- (a) streamlines coloured by speed

ax1 = nexttile(tl);  hold(ax1, 'on');
zSeed = linspace(zLim(1) + 0.02, zLim(2) - 0.02, 30);       % starting points, left edge
xSeed = xLim(1) * ones(size(zSeed));
lines = stream2(X, Z, U, W, xSeed, zSeed);                   % integrate the streamlines

for s = 1:numel(lines)
    v = lines{s};
    
    if isempty(v), continue; end
    sp = interp2(X, Z, Vmag, v(:,1), v(:,2));                % speed along the line
    patch(ax1, [v(:,1); NaN], [v(:,2); NaN], [sp; NaN], ...
          'EdgeColor', 'interp', 'FaceColor', 'none', 'LineWidth', 1.3);
end

fill(ax1, x, z, T.body, 'EdgeColor', T.fg, 'LineWidth', 0.8);
colormap(ax1, turbo);
cb1 = colorbar(ax1);  cb1.Color = T.fg;  cb1.Label.String = 'V / V_\infty';
axis(ax1, 'equal');  xlim(ax1, xLim);  ylim(ax1, zLim);
xlabel(ax1, 'x/c');  ylabel(ax1, 'z/c');
title(ax1, sprintf('Streamlines, NACA %s, \\alpha = %g°,  C_l = %.3f', code, alphaDeg, VP.Cl));
styleAxes(ax1, T);  grid(ax1, 'off');

% ---- (b) pressure coefficient field

ax2 = nexttile(tl);  hold(ax2, 'on');
contourf(ax2, X, Z, CpField, linspace(-1.5, 1.5, 61), 'LineStyle', 'none');
fill(ax2, x, z, T.body, 'EdgeColor', T.body, 'LineWidth', 3);   % thick edge hides the grid gaps
plot(ax2, x, z, '-', 'Color', T.fg, 'LineWidth', 0.8);          % thin outline on top
cmapCp = interp1([0 0.5 1], [0.23 0.30 0.75; 0.95 0.95 0.95; 0.70 0.02 0.15], ...
                 linspace(0, 1, 256));                       % blue - white - red
colormap(ax2, cmapCp);  clim(ax2, [-1.5 1.5]);
cb2 = colorbar(ax2);  cb2.Color = T.fg;  cb2.Label.String = 'C_p';
axis(ax2, 'equal');  xlim(ax2, xLim);  ylim(ax2, zLim);
xlabel(ax2, 'x/c');  ylabel(ax2, 'z/c');
title(ax2, 'Pressure coefficient field  (blue: suction, red: high pressure)');
styleAxes(ax2, T);  grid(ax2, 'off');

%% Save

if ~exist(figDir, 'dir'), mkdir(figDir); end
outFile = fullfile(figDir, sprintf('step4_flowField_NACA%s.png', code));
ax1.Toolbar.Visible = 'off';
ax2.Toolbar.Visible = 'off';
exportgraphics(fig, outFile, 'Resolution', 200, 'BackgroundColor', 'current');
fprintf('\nFigure saved to %s\n', outFile);