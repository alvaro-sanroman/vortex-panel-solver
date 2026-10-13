function VP = vortexPanelConstant(x, z, alphaDeg, dropRow)
%VORTEXPANELCONSTANT  First-order vortex panel method (Anderson 4.10).
%   VP = vortexPanelConstant(x, z, alphaDeg, dropRow)
%   x, z     : airfoil nodes, clockwise from the TE (from naca4)
%   alphaDeg : angle of attack [deg]
%   dropRow  : control point whose tangency equation is replaced by Kutta
%              (optional, default N/2 = near the LE)
%
%   Constant gamma on each panel. V_inf = 1, c = 1.
%   Solves  V_inf cos(beta_i) + sum_j A_ij gamma_j = 0  at every control point.

P = panelGeometry(x, z);                % panel geometry (Part 1)
N = P.N;                                % number of panels = number of unknowns
if nargin < 4, dropRow = round(N/2); end
alpha = deg2rad(alphaDeg);

% ---------------- influence matrix
% A(i,j) = normal velocity at control point i induced by panel j with gamma_j = 1.
% One call per panel gives its effect on ALL control points = column j.
A = zeros(N, N);
for j = 1:N
[Vx, Vz] = vortexPanelInfluenceConstant(P.xm, P.zm, ...
                        P.x1(j), P.z1(j), P.x2(j), P.z2(j));
A(:, j) = Vx .* P.nx + Vz .* P.nz;  % dot product with the normal
end

% ---------------- right-hand side: freestream normal component
% beta_i = angle between V_inf (at alpha) and the normal (at phi_i + 90 deg).
% V_inf cos(beta_i) is the flow the freestream pushes through panel i;
% the vortices must cancel it  ->  sum A_ij gamma_j = -V_inf cos(beta_i)
beta = P.phi + pi/2 - alpha;
b    = -cos(beta);

% ---------------- Kutta condition: gamma_1 + gamma_N = 0 (Anderson eq. 4.81)
% N tangency equations + Kutta = N+1 equations for N unknowns, so the
% tangency equation at control point 'dropRow' is replaced by Kutta.
A(dropRow, :)   = 0;
A(dropRow, 1)   = 1;
A(dropRow, end) = 1;
b(dropRow)      = 0;

% ---------------- solve for the vortex strengths
gamma = A \ b;

% ---------------- results
VP.P     = P;
VP.gamma = gamma;
VP.Cp    = 1 - gamma.^2;               % surface velocity = gamma (Anderson)
VP.Cl    = 2 * sum(gamma .* P.L);      % L' = rho V_inf sum(gamma_j s_j) (eq. 4.83)
end