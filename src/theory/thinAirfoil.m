function TA = thinAirfoil(theta, dzdx, alphaDeg)
%THINAIRFOIL  Thin airfoil theory for an arbitrary camber line.
%   theta    : angles from 0 (LE) to pi (TE), with x/c = (1 - cos(theta))/2
%   dzdx     : camber-line slope dz/dx at those angles
%   alphaDeg : angle(s) of attack in degrees (a single value or a vector)
%
%   Output struct TA with fields:
%   A0, A1, A2, alphaL0, alphaIdeal, Cl, Cm_c4, Cm_LE, xcp
%
%   Reference: Anderson, Fundamentals of Aerodynamics, ch. 4.

    alpha = deg2rad(alphaDeg);

% ---------------- Fourier coefficients
    %   A0 = alpha - (1/pi) * int_0^pi dz/dx dtheta
    %   An = (2/pi) * int_0^pi dz/dx cos(n theta) dtheta
A0 = alpha - 1/pi * trapz(theta, dzdx);
A1 = 2/pi * trapz(theta, dzdx .* cos(theta));
A2 = 2/pi * trapz(theta, dzdx .* cos(2*theta));

% ---------------- airfoil properties
TA.A0 = A0;
TA.A1 = A1;
TA.A2 = A2;

% zero-lift angle:  alpha_L0 = -(1/pi) * int_0^pi dz/dx (cos(theta) - 1) dtheta
TA.alphaL0    = rad2deg(-1/pi * trapz(theta, dzdx .* (cos(theta) - 1)));

% ideal angle (A0 = 0, smooth flow at the LE): alpha_i = (1/pi) * int_0^pi dz/dx dtheta
TA.alphaIdeal = rad2deg( 1/pi * trapz(theta, dzdx));

% lift and moment coefficients
TA.Cl    = pi * (2*A0 + A1);
TA.Cm_LE = -pi/2 * (A0 + A1 - A2/2);
TA.Cm_c4 = pi/4 * (A2 - A1);

% centre of pressure:  x_cp/c = 1/4 - c_m,c/4 / c_l
TA.xcp   = 1/4 - TA.Cm_c4 ./ TA.Cl;

end