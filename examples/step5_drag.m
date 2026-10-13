%% STEP 5 - Skin-friction drag (flat-plate estimate, Anderson 4.12)
%  Inviscid panel method -> lift and Cp, but Cd = 0 (d'Alembert).
%  Here drag is estimated as flat-plate skin friction, cd = 2*Cf.

clear; clc; close all
root   = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root, 'src')));
figDir = fullfile(root, 'figures');

% ----------------------------------------------------------- inputs

code     = '2412';
alphaDeg = 4;
nPanels  = 200;
V        = 50;          % freestream speed       [m/s]
c        = 1.0;         % chord                  [m]
ReCrit   = 5e5;         % transition Reynolds    [-]

% ISA sea level

rho = 1.225;            % density          [kg/m^3]
mu  = 1.789e-5;         % viscosity        [Pa s]
a   = 340.3;            % speed of sound   [m/s]

% ------------------------------- validation: Anderson Examples 4.8-4.11
% Same NACA 2412 at Re = 3.1e6; measured cd = 0.0068 (Anderson Fig. 4.11)

ReA    = 3.1e6;
cdMeas = 0.0068;
[cL, cT, cM5] = skinFriction(ReA, 5e5);
[~,  ~,  cM1] = skinFriction(ReA, 1e6);
fprintf('Validation vs Anderson (NACA 2412, Re = 3.1e6):\n');
fprintf('  laminar          : %.5f   (Ex. 4.8 : 0.0015)\n',  cL);
fprintf('  turbulent        : %.5f   (Ex. 4.9 : 0.00744)\n', cT);
fprintf('  mixed, Recr=5e5  : %.5f   (Ex. 4.10: 0.0063)\n',  cM5);
fprintf('  mixed, Recr=1e6  : %.5f   (Ex. 4.11: 0.00528)\n', cM1);
fprintf('  measured         : %.4f    -> pressure drag ~ %.0f%% of total\n\n', ...
        cdMeas, 100*(cdMeas - cM5)/cdMeas);

% ------------------------------------------------- flight condition

Re = rho * V * c / mu;
M  = V / a;
q  = 0.5 * rho * V^2;   % dynamic pressure [Pa]
fprintf('V = %.1f m/s, c = %.2f m  ->  Re = %.3e, M = %.3f\n', V, c, Re, M);

if M > 0.3
    warning('M > 0.3: incompressible assumption no longer valid.');
end

% ------------------------------------------------ inviscid solution

[x, z] = naca4(code, nPanels);
VP = vortexPanel(x, z, alphaDeg);

% ----------------------------------------------------- friction drag

[cdL, cdT, cdM] = skinFriction(Re, ReCrit);
fprintf('\nSkin-friction cd (flat plate, 2 faces):\n');
fprintf('  laminar   : %.5f\n', cdL);
fprintf('  turbulent : %.5f\n', cdT);
fprintf('  mixed     : %.5f   (transition at x/c = %.3f)\n', cdM, min(ReCrit/Re, 1));

% CdP from the panel method is only discretisation error (should be 0),
% so it is NOT added to the drag; it is printed as an accuracy check.

fprintf('\nPanel CdP (numerical error, ideally 0): %.5f\n', VP.CdP);

% ------------------------------------------------------ final numbers

cd = cdM;                                   % chosen estimate
fprintf('\nNACA %s, alpha = %g deg:\n', code, alphaDeg);
fprintf('  Cl = %.4f   cd = %.5f   L/D = %.1f\n', VP.Cl, cd, VP.Cl/cd);
fprintf('  per unit span: L'' = %.0f N/m,  D'' = %.1f N/m\n', q*c*VP.Cl, q*c*cd);

% ------------------------------------------------------------- plot

ReV = logspace(5, 8, 300)';
[cdLv, cdTv, cdMv] = skinFriction(ReV, ReCrit);

T   = plotTheme();
fig = figure('Color', T.bg, 'Position', [100 100 900 520], 'WindowStyle', 'normal');
ax  = axes(fig); hold(ax, 'on');
plot(ax, ReV, cdLv, '-',  'Color', T.accent1, 'LineWidth', 1.8);
plot(ax, ReV, cdTv, '-',  'Color', T.accent3, 'LineWidth', 1.8);
plot(ax, ReV, cdMv, '--', 'Color', T.accent2, 'LineWidth', 1.8);
plot(ax, ReA, cdMeas, 's', 'MarkerSize', 9, ...
     'MarkerFaceColor', T.accent4, 'MarkerEdgeColor', T.fg);
set(ax, 'XScale', 'log', 'YScale', 'log');
xlabel(ax, 'Re_c'); ylabel(ax, 'c_d = 2 C_f');
title(ax, 'Skin-friction drag - flat-plate estimate vs NACA 2412 data');
lg = legend(ax, {'laminar', 'turbulent', ...
            sprintf('mixed (Re_{cr} = %.0e)', ReCrit), ...
            'NACA 2412 measured (Anderson Fig. 4.11)'}, 'Location', 'northeast');
set(lg, 'TextColor', T.fg, 'Color', T.bg, 'EdgeColor', T.grid);
styleAxes(ax, T);
ax.Toolbar.Visible = 'off';

exportgraphics(fig, fullfile(figDir, 'step5_skinFriction.png'), ...
               'Resolution', 200, 'BackgroundColor', 'current');