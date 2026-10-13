function [x, z, xc, zc, dzc, beta, K] = naca4mod(code, nPanels, closedTE)
%NACA4MOD  Coordinates of a NACA modified 4-digit airfoil (e.g. '0009-64').
%
%   [x, z] = naca4mod('0009-64', 200)
%   [x, z, xc, zc, dzc, beta, K] = naca4mod(code, nPanels, closedTE)
%
%   code     : 'MPTT-IX'
%                M, P : camber digits, same as NACA 4-digit
%                TT   : max thickness (% chord)
%                I    : leading-edge radius index (6 = normal 4-digit radius)
%                X    : position of max thickness, tenths of chord (2..6)
%   nPanels  : number of panels (even)
%   closedTE : true  -> d0 = 0, sharp TE (default)
%              false -> original definition, d0 = 0.002
%   K        : thickness constants (K.a = a0..a3, K.d = d0..d3)
%
%   All coordinates are non-dimensional (divided by the chord c), i.e. the
%   chord is c = 1: x means x/c and z means z/c. For real dimensions,
%   multiply x and z by the chord.
%
%   Same node ordering as naca4: TE -> lower -> LE -> upper -> TE.
%
%   Reference: Ladson & Brooks, NASA TM X-3284 (1975):
%              thickness equations and constraints pp. 4-6,
%              camber line pp. 6-7, LE radius values Table I.
%   Note: the report calls the max-thickness position "m"; here it is
%   "xt" (m is the max camber, as in naca4). The report's y is our z.

if nargin < 3, closedTE = true; end

% ----------------------------------------------------------- input

code = char(code);

if numel(code) ~= 7 || code(5) ~= '-' || any(~isstrprop(code([1:4 6 7]), 'digit'))
   
    error('naca4mod:badCode', 'The code must look like ''0009-64''.');

end

if mod(nPanels, 2) ~= 0
    
    error('naca4mod:oddPanels', 'nPanels must be even.');

end

m  = str2double(code(1))   / 100;   % max camber                (z/c)
p  = str2double(code(2))   / 10;    % position of max camber    (x/c)
t  = str2double(code(3:4)) / 100;   % max thickness             (t/c)
I  = str2double(code(6));           % leading-edge radius index
xt = str2double(code(7))   / 10;    % position of max thickness (x/c)

if xt < 0.2 || xt > 0.6
   error('naca4mod:badXt', 'Max-thickness digit must be 2 to 6.');
end

% -------------------------------------------------- cosine spacing
% xc = x/c, from 0 (LE) to 1 (TE)

nSide = nPanels / 2;
beta  = linspace(0, pi, nSide + 1)';
xc    = 0.5 * (1 - cos(beta));

% ------------------------------------- known thickness constants
% Equations are written for the 20 % thick base section (max y/c = 0.1)
% and scaled by (t/c)/0.20 at the end (p. 6).

% Constraint (2): LE radius R = a0^2/2.
% a0 is proportional to the index; index 6 -> a0 = 0.2969 (p. 6)

if I == 9

   a0 = 0.2969 * sqrt(3);      % index 9 = 3x normal LE radius (p. 6)

else

   a0 = 0.2969 * I / 6;

end

% Constraint (4): TE ordinate y/c = d0 = 0.002 (p. 5).
% d0 = 0 is our change to close the TE (needed for the Kutta condition).

if closedTE, d0 = 0; else, d0 = 0.002; end

% Constraint (5): TE slope d1 = f(xt), table on p. 5

d1 = interp1([0.2 0.3 0.4 0.5 0.6], ...
                 [0.200 0.234 0.315 0.465 0.700], xt);

% ------------------------------- aft part: solve for d2 and d3
% y/c = d0 + d1*u + d2*u^2 + d3*u^3,  u = 1 - x/c
% Constraint (1) at x/c = xt (u = U):  y/c = 0.1,  dy/dx = 0

U  = 1 - xt;
dd = [U^2,   U^3;
      2*U, 3*U^2] \ [0.1 - d0 - d1*U; -d1];
d2 = dd(1);  d3 = dd(2);

% Aft curvature at xt. Constraint (3): both parts have the same radius
% of curvature there (p. 5). The report's R = (1-xt)^2/(2*d1*(1-xt) - 0.588)
% is 1/|yppAft| with d0 = 0.002 (0.588 = 6*(0.1 - d0)).
    
yppAft = 2*d2 + 6*d3*U;

% -------------------------- forward part: solve for a1, a2, a3
% y/c = a0*sqrt(x/c) + a1*(x/c) + a2*(x/c)^2 + a3*(x/c)^3
% at x/c = xt:  y/c = 0.1,  dy/dx = 0   (constraint 1)
%               d2y/dx2 = yppAft        (constraint 3)

Aeq = [xt,  xt^2,  xt^3;
       1,   2*xt,  3*xt^2;
       0,   2,     6*xt];

beq = [0.1 - a0*sqrt(xt);
       -0.5*a0/sqrt(xt);
       yppAft + 0.25*a0*xt^(-1.5)];

aa = Aeq \ beq;
a1 = aa(1);  a2 = aa(2);  a3 = aa(3);

K.a = [a0 a1 a2 a3];
K.d = [d0 d1 d2 d3];

% --------------------------------------- thickness distribution

zt  = zeros(size(xc));          % half-thickness z_t/c
fwd = xc < xt;
aft = ~fwd;
zt(fwd) = a0*sqrt(xc(fwd)) + a1*xc(fwd) + a2*xc(fwd).^2 + a3*xc(fwd).^3;
u       = 1 - xc(aft);
zt(aft) = d0 + d1*u + d2*u.^2 + d3*u.^3;
zt      = zt * t / 0.2;         % scale base section to (t/c)

% ---------------------------- camber line (pp. 6-7, as in naca4)

zc  = zeros(size(xc));          % z_c/c
dzc = zeros(size(xc));          % slope dz_c/dx

if m > 0 && p > 0

    f = xc < p;   r = ~f;
    zc(f)  = m / p^2       * (2*p*xc(f) - xc(f).^2);
    zc(r)  = m / (1 - p)^2 * ((1 - 2*p) + 2*p*xc(r) - xc(r).^2);
    dzc(f) = 2*m / p^2       * (p - xc(f));
    dzc(r) = 2*m / (1 - p)^2 * (p - xc(r));

end

theta = atan(dzc);

% ------------------- thickness perpendicular to the camber line
xU = xc - zt.*sin(theta);   zU = zc + zt.*cos(theta);
xL = xc + zt.*sin(theta);   zL = zc - zt.*cos(theta);

% --------------------------- assemble clockwise, starting at TE
    
x = [flipud(xL); xU(2:end)];
z = [flipud(zL); zU(2:end)];

end