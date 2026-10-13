function [Vxa, Vza, Vxb, Vzb] = vortexPanelInfluenceLinear(xp, zp, x1, z1, x2, z2)
%VORTEXPANELINFLUENCELINEAR  Velocity induced by a linear-strength vortex panel.
%   [Vxa, Vza, Vxb, Vzb] = vortexPanelInfluenceLinear(xp, zp, x1, z1, x2, z2)
%   xp, zp : points where the velocity is wanted (vectors)
%   x1, z1 : panel start (node j)      x2, z2 : panel end (node j+1)
%
%   gamma varies linearly along the panel (Anderson Fig. 4.34):
%       gamma(x') = gamma_j (1 - x'/L) + gamma_j+1 (x'/L)
%   so the induced velocity is
%       V = gamma_j * (Vxa, Vza) + gamma_j+1 * (Vxb, Vzb)
%   Outputs are per unit gamma_j and gamma_j+1 (the (.)^a and (.)^b parts).
%
%   Reference: Katz & Plotkin, Low-Speed Aerodynamics, 2nd ed., Sec. 10.3.3,
%   eqs. (10.72)-(10.73), and Sec. 11.4.2, eqs. (11.97)-(11.101).

    % ---------------- panel length and orientation
    L = hypot(x2 - x1, z2 - z1);         % panel length
    c = (x2 - x1) / L;                   % cos(phi)
    s = (z2 - z1) / L;                   % sin(phi)

    % ---------------- P in local panel coordinates (x' along, z' normal)
    dx = xp - x1;
    dz = zp - z1;
    xl =  dx * c + dz * s;
    zl = -dx * s + dz * c;

    % ---------------- angle subtended by the panel and log of distances
    dTheta = atan2(zl, xl - L) - atan2(zl, xl);          % theta2 - theta1
    lnR    = log(hypot(xl, zl) ./ hypot(xl - L, zl));    % ln(r1/r2)

    % ---------------- points ON the panel: take the outer side (as before)
    onPanel = abs(zl) < 1e-12 & xl > 0 & xl < L;
    dTheta(onPanel) = pi;
    zl(onPanel)     = 0;                 % remove rounding noise (~1e-17)

    % ---------------- constant part: gamma = 1   (K&P eqs. 10.39-10.40)
    u0 =  dTheta / (2*pi);
    w0 = -lnR    / (2*pi);

    % ---------------- linear part: gamma = x'   (K&P eqs. 10.72-10.73)
    u1 =  (xl .* dTheta - zl .* lnR)     / (2*pi);
    w1 = -(xl .* lnR - L + zl .* dTheta) / (2*pi);

    % ---------------- split into the two nodes   (K&P eq. 11.101)
    % (.)^a: gamma = 1 - x'/L  (1 at node j,   0 at node j+1)
    % (.)^b: gamma = x'/L      (0 at node j,   1 at node j+1)
    ua = u0 - u1/L;    wa = w0 - w1/L;
    ub = u1/L;         wb = w1/L;

    % ---------------- back to global axes
    Vxa = ua * c - wa * s;    Vza = ua * s + wa * c;
    Vxb = ub * c - wb * s;    Vzb = ub * s + wb * c;
end