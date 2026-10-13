function [U, W, inside] = flowField(VP, X, Z, alphaDeg)
%FLOWFIELD  Velocity field around the airfoil from a vortexPanel solution.
%   [U, W, inside] = flowField(VP, X, Z, alphaDeg)
%   VP       : solution returned by vortexPanel
%   X, Z     : grid of points (e.g. from meshgrid)
%   alphaDeg : angle of attack used in vortexPanel [deg]
%
%   U, W     : velocity components at each grid point (V_inf = 1)
%   inside   : true for points inside the airfoil (no flow there)
%
%   Superposition (K&P Sec. 3.9): freestream + every linear-vortex panel,
%   each evaluated with K&P eq. (11.101) at the grid points.

alpha = deg2rad(alphaDeg);
P     = VP.P;
g     = VP.gamma;

% ---------------- work with columns of points

xp = X(:);
zp = Z(:);

% ---------------- freestream
U = cos(alpha) * ones(size(xp));
W = sin(alpha) * ones(size(xp));

% ---------------- add every panel with the solved gamma values

for j = 1:P.N
     [Vxa, Vza, Vxb, Vzb] = vortexPanelInfluenceLinear(xp, zp, ...
                                   P.x1(j), P.z1(j), P.x2(j), P.z2(j));
      U = U + g(j) * Vxa + g(j+1) * Vxb;
      W = W + g(j) * Vza + g(j+1) * Vzb;
end

% ---------------- points inside the airfoil contour

inside = inpolygon(xp, zp, [P.x1; P.x2(end)], [P.z1; P.z2(end)]);

% ---------------- back to the grid shape
U      = reshape(U, size(X));
W      = reshape(W, size(X));
inside = reshape(inside, size(X));
end