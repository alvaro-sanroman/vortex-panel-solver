%% STEP 7 - F-86F "6-3" wing sections at approach and touchdown
%  2D, inviscid lift (panel method) + flat-plate friction drag.
%  Limitations: no flaps, no sweep; CL = W/(q*S) is a 3D aircraft value,
%  so the 2D angle of attack obtained is only illustrative.

clear; clc; close all
root   = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root, 'src')));
figDir = fullfile(root, 'figures');

% ---------------------------------------------------- unit factors

kt = 0.514444;            % knot -> m/s
ft = 0.3048;              % foot -> m
lb = 4.44822;             % pound-force -> N

% ------------------------------------------------- aircraft data
% T.O. 1F-86F-1, Fig. 2-6 (p. 2-17): airplanes without slats, 13,800 lb

W     = 13800 * lb;                   % weight [N]
cases = {'Final approach', 'Touchdown'};
Vkt   = [145, 120];                   % IAS [kt]  (IAS ~ TAS at sea level)

% NACA RM A56G30, Table I and footnote p. 5 ("6-3" wing)

S     = 302 * ft^2;                   % wing area [m^2]
codes = {'0009-64', '0008-64'};       % root, tip (Lednicer, streamwise)
names = {'root', 'tip'};
chord = [10.77, 5.49] * ft;           % root, tip chord [m] (derived)

% ISA sea level

rho = 1.225;  mu = 1.789e-5;  a = 340.3;

nPanels  = 200;
alphaDeg = 0:1:10;

% ------------------------------- Cl-alpha of each section (2D)

Cl    = zeros(numel(alphaDeg), 2);
slope = zeros(1, 2);

for s = 1:2
    [x, z] = naca4mod(codes{s}, nPanels);

    for k = 1:numel(alphaDeg)
        VP = vortexPanel(x, z, alphaDeg(k));
        Cl(k, s) = VP.Cl;

    end

    pf       = polyfit(alphaDeg, Cl(:, s)', 1);   % Cl = slope*alpha + Cl0
    slope(s) = pf(1);                              % [1/deg]
    fprintf('NACA %s (%s): Cl_alpha = %.4f /deg = %.3f /rad\n', ...
            codes{s}, names{s}, slope(s), slope(s)*180/pi);
end

% -------------------------------------------- flight conditions

CLreq = zeros(1, 2);

for c = 1:2
    
    V = Vkt(c) * kt;
    M = V / a;
    q = 0.5 * rho * V^2;
    CLreq(c) = W / (q * S);                        % required aircraft CL

    fprintf('\n%s: %d kt = %.1f m/s, M = %.3f, q = %.0f Pa, CL required = %.3f\n', ...
            cases{c}, Vkt(c), V, M, q, CLreq(c));

    if M > 0.3
        
        warning('M > 0.3: incompressible assumption no longer valid.');
    
    end

    for s = 1:2
        
        Re        = rho * V * chord(s) / mu;
        [~, ~, cd] = skinFriction(Re);             % mixed, Re_cr = 5e5
        aReq      = CLreq(c) / slope(s);           % 2D alpha for that CL
        fprintf('  %-4s NACA %s: c = %.2f m, Re = %.2e, cd = %.5f, alpha_2D = %.1f deg\n', ...
                names{s}, codes{s}, chord(s), Re, cd, aReq);
    
    end
    
end

% ------------------------------------------------------------- plot
T   = plotTheme();
fig = figure('Color', T.bg, 'Position', [100 100 900 560], 'WindowStyle', 'normal');
ax  = axes(fig); hold(ax, 'on');
plot(ax, alphaDeg, Cl(:, 1), '-o', 'Color', T.accent1, 'LineWidth', 1.6, ...
     'MarkerFaceColor', T.accent1, 'MarkerSize', 4);
plot(ax, alphaDeg, Cl(:, 2), '-s', 'Color', T.accent2, 'LineWidth', 1.6, ...
     'MarkerFaceColor', T.accent2, 'MarkerSize', 4);
yline(ax, CLreq(1), '--', sprintf('Final approach (%d kt)', Vkt(1)), ...
      'Color', T.accent3, 'LineWidth', 1.2, 'LabelHorizontalAlignment', 'left');
yline(ax, CLreq(2), '--', sprintf('Touchdown (%d kt)', Vkt(2)), ...
      'Color', T.accent4, 'LineWidth', 1.2, 'LabelHorizontalAlignment', 'left');
xlabel(ax, '\alpha [deg]'); ylabel(ax, 'C_l');
title(ax, 'F-86F "6-3" wing sections - 2D lift vs required C_L');
lg = legend(ax, {'root, NACA 0009-64', 'tip, NACA 0008-64'}, 'Location', 'southeast');
set(lg, 'TextColor', T.fg, 'Color', T.bg, 'EdgeColor', T.grid);
styleAxes(ax, T);
ax.Toolbar.Visible = 'off';

exportgraphics(fig, fullfile(figDir, 'step7_f86Approach.png'), ...
               'Resolution', 200, 'BackgroundColor', 'current');