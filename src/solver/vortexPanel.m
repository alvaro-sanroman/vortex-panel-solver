function VP = vortexPanel(x, z, alphaDeg)
%VORTEXPANEL  Second-order (linear-strength) vortex panel method.
%   VP = vortexPanel(x, z, alphaDeg)
%   x, z     : airfoil nodes, clockwise from the TE (from naca4)
%   alphaDeg : angle of attack [deg]
%
%   gamma varies linearly on each panel; the unknowns are the N+1 nodal
%   values. N tangency conditions + Kutta = N+1 equations. V_inf = 1, c = 1.
%
%   References: Anderson, Fundamentals of Aerodynamics, Sec. 1.5, 4.10 (Fig. 4.34);
%               Katz & Plotkin, Low-Speed Aerodynamics, 2nd ed., Sec. 11.4.2.

P = panelGeometry(x, z);
N = P.N;                              % panels (unknowns = N+1 nodes)
alpha = deg2rad(alphaDeg);

% ---------------- influence matrices (K&P step c, eq. 11.103)
% A(i,j)  = NORMAL velocity at control point i per unit gamma at node j
% At(i,j) = TANGENTIAL velocity, same idea (needed for Cp)
% Panel j touches nodes j and j+1: its (.)^a part goes to column j,
% its (.)^b part to column j+1. Middle nodes get contributions from two
% panels, hence the '+' (accumulate, don't overwrite).
A  = zeros(N+1, N+1);
At = zeros(N,   N+1);
for j = 1:N
[Vxa, Vza, Vxb, Vzb] = vortexPanelInfluenceLinear(P.xm, P.zm, ...
                                   P.x1(j), P.z1(j), P.x2(j), P.z2(j));
A(1:N, j)   = A(1:N, j)   + Vxa .* P.nx + Vza .* P.nz;   % normal comp.
A(1:N, j+1) = A(1:N, j+1) + Vxb .* P.nx + Vzb .* P.nz;
At(:, j)    = At(:, j)    + Vxa .* P.tx + Vza .* P.tz;   % tangential comp.
At(:, j+1)  = At(:, j+1)  + Vxb .* P.tx + Vzb .* P.tz;
end

% ---------------- right-hand side (K&P step d)
% V_inf cos(beta_i) + sum A_ij gamma_j = 0  ->  sum A_ij gamma_j = -cos(beta_i)
beta = P.phi + pi/2 - alpha;          % angle between V_inf and the normal
b    = [-cos(beta); 0];               % last entry: Kutta row

% ---------------- Kutta condition as an EXTRA row (K&P eq. 11.105)
% gamma_1 + gamma_N+1 = 0 -> square (N+1)x(N+1) system, no row dropped
A(N+1, 1)   = 1;
A(N+1, N+1) = 1;

% ---------------- solve (K&P step e)
gamma = A \ b;

% ---------------- surface velocity and Cp (K&P step f, eqs. 11.106, 11.18)
% tangential velocity = freestream along the panel (V_inf sin(beta_i))
%                     + velocity induced by all the nodes (At * gamma)
Vt = sin(beta) + At * gamma;
Cp = 1 - Vt.^2;                       % Bernoulli

% ---------------- lift (1): Kutta-Joukowski  [main result]
% L' = rho V_inf * Gamma, with Gamma = sum(mean gamma on panel * length)
% Anderson eq. (4.83), K&P eq. (11.107)
gammaMean = 0.5 * (gamma(1:N) + gamma(2:N+1));
ClKJ      = 2 * sum(gammaMean .* P.L);

% ---------------- lift (2): pressure integration  [independent check]
% Anderson Sec. 1.5: pressure pushes inward (-n) on each panel ->
% normal force N' and axial force A', then
%   L' = N' cos(alpha) - A' sin(alpha),   D' = N' sin(alpha) + A' cos(alpha)
dFx = -Cp .* P.L .* P.nx;             % force on each panel, x component
dFz = -Cp .* P.L .* P.nz;             % force on each panel, z component
Fx  = sum(dFx);                       % axial force A'
Fz  = sum(dFz);                       % normal force N'
ClP = Fz * cos(alpha) - Fx * sin(alpha);
CdP = Fz * sin(alpha) + Fx * cos(alpha);   % inviscid -> should tend to 0

% ---------------- pitching moment about the quarter chord (Anderson Sec. 1.5)
% moment of each panel force about (x_ref, 0), positive nose-up:
%   dM = z*dFx - (x - x_ref)*dFz
xRef  = 0.25;
Cm_c4 = sum(P.zm .* dFx - (P.xm - xRef) .* dFz);
% ---------------- results
VP.P     = P;
VP.gamma = gamma;           % N+1 nodal values
VP.Vt    = Vt;              % surface velocity at the control points
VP.Cp    = Cp;
VP.Cl    = ClKJ;
VP.ClP   = ClP;
VP.CdP   = CdP;
VP.Cm_c4 = Cm_c4;

end

