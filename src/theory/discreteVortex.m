function DV = discreteVortex(xc, dzc, alphaDeg, N)
%DISCRETEVORTEX  Discrete (lumped) vortex method on the mean camber line.
%   DV = discreteVortex(xc, dzc, alphaDeg, N)
%   xc, dzc  : camber line x/c and slope dz/dx (e.g. from naca4)
%   alphaDeg : angle(s) of attack [deg]
%   N        : number of discrete vortices
%
%   Each of the N equal segments gets a vortex at its 1/4 point and a
%   control point at its 3/4 point. Flow tangency at the control points
%   gives A*Gamma = b. Kutta is satisfied automatically. V_inf = 1, c = 1.

alpha = deg2rad(alphaDeg(:)');        % row vector of angles

% vortex and control point positions
edges = linspace(0, 1, N + 1)';       % segment end points
h     = diff(edges);                  % segment lengths
xv    = edges(1:end-1) + h/4;         % vortices at the 1/4 points
xCtrl = edges(1:end-1) + 3*h/4;       % control points at the 3/4 points
dzdx  = interp1(xc, dzc, xCtrl);      % camber slope at the control points

% influence matrix and right-hand side
% A(i,j) = velocity induced at control point i by a unit vortex at j
A = 1 ./ (2*pi * (xCtrl - xv'));      % N x N
b = alpha - dzdx;                     % V_inf * (alpha - dz/dx)

% solve for the vortex strengths
Gamma = A \ b;

    % ---------------- coefficients
DV.xv    = xv;
DV.xCtrl = xCtrl;
DV.Gamma = Gamma;
DV.Cl    = 2 * sum(Gamma, 1);                % Kutta-Joukowski
DV.Cm_LE = -2 * sum(Gamma .* xv, 1);
DV.Cm_c4 = DV.Cm_LE + DV.Cl/4;
end