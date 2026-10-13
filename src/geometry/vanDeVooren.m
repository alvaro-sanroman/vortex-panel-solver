function [x, z, Cp, Cl] = vanDeVooren(epsilon, k, nPanels, alphaDeg)
%VANDEVOOREN  Van de Vooren airfoil and its exact potential-flow solution.
%   [x, z, Cp, Cl] = vanDeVooren(epsilon, k, nPanels, alphaDeg)
%   epsilon  : thickness parameter (0.10 with k = 1.97 -> 15% thick)
%   k        : trailing-edge angle parameter, tau = pi*(2 - k)
%   nPanels  : number of panels (nodes = nPanels + 1)
%   alphaDeg : angle of attack [deg]
%
%   x, z : nodes, chord 1, clockwise from the TE (same convention as naca4)
%   Cp   : exact pressure coefficient at the nodes
%   Cl   : exact lift coefficient
%
%   Reference: Katz & Plotkin, Low-Speed Aerodynamics, 2nd ed., Sec. 6.6-6.7,
%   eqs. (6.62)-(6.75), and Appendix D, Program No. 1.

ell   = 1;                                   % semichord (chord = 2*ell)
a     = 2*ell*(1 + epsilon)^(k-1) / 2^k;     % circle radius, eq. (6.65)
alpha = deg2rad(alphaDeg);

% ---------------- points on the circle: equal steps in theta
% theta = 0 is the TE, theta = pi the LE, 0 < theta < pi the upper surface.
% Going from 2*pi down to 0 gives clockwise numbering: TE -> lower -> LE -> upper -> TE
theta = linspace(2*pi, 0, nPanels + 1)';
th    = theta(2:end-1);                      % interior points (TE added at the end)

% ---------------- distances and angles in the circle plane (eq. 6.74)
r1  = sqrt((a*(cos(th) - 1)).^2       + (a*sin(th)).^2);
r2  = sqrt((a*(cos(th) - epsilon)).^2 + (a*sin(th)).^2);
th1 = atan(sin(th) ./ (cos(th) - 1)) + pi;
th2 = mod(atan2(sin(th), cos(th) - epsilon), 2*pi);   % quadrant logic of eq. (6.74)

% ---------------- airfoil coordinates (eq. 6.74, Program 1)
com1 = r1.^k ./ r2.^(k-1);
X = com1 .* (cos(k*th1).*cos((k-1)*th2) + sin(k*th1).*sin((k-1)*th2)) + ell;
Z = com1 .* (sin(k*th1).*cos((k-1)*th2) - cos(k*th1).*sin((k-1)*th2));

% ---------------- exact surface velocity (eq. 6.75, Program 1)
A  = cos((k-1)*th1).*cos(k*th2) + sin((k-1)*th1).*sin(k*th2);
B  = sin((k-1)*th1).*cos(k*th2) - cos((k-1)*th1).*sin(k*th2);
D0 = a*(1 - k + k*epsilon);
D1 = A.*(a*cos(th) - D0) - B.*a.*sin(th);
D2 = A.*a.*sin(th) + B.*(a*cos(th) - D0);
fac = 2 * (r2.^k ./ r1.^(k-1)) .* (sin(alpha) - sin(alpha - th)) ./ (D1.^2 + D2.^2);
u  =  fac .* (D1.*sin(th) + D2.*cos(th));
w  = -fac .* (D1.*cos(th) - D2.*sin(th));

% ---------------- assemble, adding the TE at both ends; scale to chord = 1
% the TE (finite angle) is a stagnation point: Cp = 1 there
x  = [1; (X + ell) / (2*ell); 1];
z  = [0;  Z        / (2*ell); 0];
Cp = [1; 1 - (u.^2 + w.^2); 1];               % Bernoulli, eq. (6.70)

% ---------------- exact lift: circulation that puts the rear stagnation
% point at the TE (Kutta), Gamma = 4 pi a V_inf sin(alpha); Cl = 2 Gamma / c
Cl = 4*pi*a*sin(alpha) / ell;
end