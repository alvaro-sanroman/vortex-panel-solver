function [Vx, Vz] = vortexPanelInfluenceConstant(xp, zp, x1, z1, x2, z2)
%VORTEXPANELINFLUENCECONSTANT  Velocity induced by a constant-strength vortex panel.
%   [Vx, Vz] = vortexPanelInfluenceConstant(xp, zp, x1, z1, x2, z2)
%   xp, zp : points where the velocity is wanted (vectors)
%   x1, z1 : panel start      x2, z2 : panel end
%
%   Velocity for UNIT strength (gamma = 1, clockwise positive). Since velocity
%   is proportional to gamma, this gives the coefficients A_ij that multiply
%   the unknowns gamma_j in  V_inf cos(beta_i) + sum A_ij gamma_j = 0.

% ---------------- panel length and orientation
L = hypot(x2 - x1, z2 - z1);         % panel length
c = (x2 - x1) / L;                   % cos(phi)
s = (z2 - z1) / L;                   % sin(phi)

% ---------------- P in local panel coordinates (x' along, z' normal)
dx = xp - x1;
dz = zp - z1;
xl =  dx * c + dz * s;
zl = -dx * s + dz * c;

% ---------------- angles and distances from P to the panel ends
theta1 = atan2(zl, xl);              % to the start
theta2 = atan2(zl, xl - L);          % to the end
r1     = hypot(xl, zl);
r2     = hypot(xl - L, zl);
dTheta = theta2 - theta1;

% ---------------- points ON the panel (its own control point)
% atan2 can't tell above (+pi) from below (-pi) the sheet there,
% so we force +pi: the outer side, where the flow is.
onPanel = abs(zl) < 1e-12 & xl > 0 & xl < L;
dTheta(onPanel) = pi;

% ---------------- induced velocity, local coordinates
u =  dTheta / (2*pi);                % tangential
w = -log(r1 ./ r2) / (2*pi);         % normal

% ---------------- back to global axes
Vx = u * c - w * s;
Vz = u * s + w * c;
end