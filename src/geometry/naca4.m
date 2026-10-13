function [x, z, xc, zc, dzc, beta] = naca4(code, nPanels, closedTE)
%   [x, z, xc, zc] = naca4(code, nPanels, closedTE)
%   code     : the 4 digits as text, e.g. '2412'
%   nPanels  : number of panels (even number)
%   closedTE : true -> sharp trailing edge (default)

if nargin < 3 %is less than 3 inputs is given create closedTE
    closedTE = true;
end

code = char(code);

if numel(code) ~= 4 || any(~isstrprop(code, 'digit'))
    error('NACA code must have exactly 4 digits, e.g. ''0012''.');
%if the code doesn't have 4 characters, 
% OR any of its characters isn't a digit error

end

if mod(nPanels, 2) ~= 0  % "if the remainder is NOT 0" = "if it's odd"
    error('nPanels mus be even');
end

m = str2double(code(1))   / 100;   % max camber
p = str2double(code(2))   / 10;    % position of max camber
t = str2double(code(3:4)) / 100;   % max thickness

    % ---------------- cosine spacing
    nSide = nPanels / 2;
    beta  = linspace(0, pi, nSide + 1)';
    xc    = 0.5 * (1 - cos(beta));

if closedTE
    a4 = -0.1036; % coefficients sum to zero -> sharp TE
else
    a4 = -0.1015;   % original NACA value -> small TE gap
end

zt = 5 * t * ( 0.2969 * sqrt(xc) - 0.1260 * xc - 0.3516 * xc.^2 ...
                 + 0.2843 * xc.^3   + a4     * xc.^4 );



zc  = zeros(size(xc));        % start with a flat line (symmetric airfoil)
    dzc = zeros(size(xc));        % and zero slope
if m > 0 && p > 0
        fwd = xc <  p;            % points in front of max camber
        aft = xc >= p;            % points behind max camber

        zc(fwd)  = m / p^2       * (2*p*xc(fwd) - xc(fwd).^2);
        zc(aft)  = m / (1 - p)^2 * ((1 - 2*p) + 2*p*xc(aft) - xc(aft).^2);

        dzc(fwd) = 2*m / p^2       * (p - xc(fwd));
        dzc(aft) = 2*m / (1 - p)^2 * (p - xc(aft));
    end

% ---------------- upper and lower surfaces
theta = atan(dzc);                 % angle of the camber line at each point
xU = xc - zt .* sin(theta);        % upper surface
zU = zc + zt .* cos(theta);
xL = xc + zt .* sin(theta);        % lower surface
zL = zc - zt .* cos(theta);

% ---------------- assemble: TE -> lower -> LE -> upper -> TE (clockwise)
x = [flipud(xL); xU(2:end)];
z = [flipud(zL); zU(2:end)];

end