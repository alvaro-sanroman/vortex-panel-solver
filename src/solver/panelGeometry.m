function P = panelGeometry(x, z)
%PANELGEOMETRY  Geometry of the panels defined by the nodes (x, z).
%   P = panelGeometry(x, z)
%   x, z : node coordinates, clockwise from the TE (as returned by naca4)
%
%   Panel j goes from node j to node j+1. Output struct P with fields:
%   N        number of panels
%   x1,z1    start node of each panel      x2,z2  end node
%   xm,zm    control points (panel midpoints)
%   L        panel lengths (S_j in Anderson)
%   phi      panel angle from the x axis
%   tx,tz    unit tangent (direction of the numbering)
%   nx,nz    unit outward normal

x = x(:);  z = z(:);        % force column vectors, whatever shape came in
P.N = numel(x) - 1;         % N+1 nodes define N panels (fence posts)

% ---------------- panel end points
% panel j runs from node j to node j+1
P.x1 = x(1:end-1);          % start of each panel: every node except the last
P.z1 = z(1:end-1);
P.x2 = x(2:end);            % end of each panel: every node except the first
P.z2 = z(2:end);

% ---------------- control points: panel midpoints
% flow tangency (V_inf cos(beta_i) + sum A_ij gamma_j = 0) is imposed here
P.xm = 0.5 * (P.x1 + P.x2);
P.zm = 0.5 * (P.z1 + P.z2);

% ---------------- length and angle
dx    = P.x2 - P.x1;        % horizontal extent of each panel
dz    = P.z2 - P.z1;        % vertical extent of each panel
P.L   = hypot(dx, dz);      % panel length S_j = sqrt(dx^2 + dz^2)
P.phi = atan2(dz, dx);      % panel angle phi_j from the x axis, in (-pi, pi]
                            % atan2 (not atan) tells left-pointing panels
                            % (lower surface, phi ~ 180 deg) from right-pointing ones

% ---------------- unit tangent and outward normal
% tangent t = (cos phi, sin phi), pointing in the numbering direction
P.tx = cos(P.phi);
P.tz = sin(P.phi);
% normal n = tangent rotated 90 deg counter-clockwise = (-sin phi, cos phi)
% with clockwise numbering this points OUT of the airfoil, into the flow
% note: angle between n and V_inf is beta_i = phi_i + pi/2 - alpha
P.nx = -sin(P.phi);
P.nz =  cos(P.phi);
end

