function [cdL, cdT, cdM] = skinFriction(Re, ReCrit)
%SKINFRICTION  Skin-friction drag estimate: flat plate, both sides wetted.
%
%   Re     : chord Reynolds number (scalar or vector)
%   ReCrit : transition Reynolds number (default 5e5)
%   cdL, cdT, cdM : cd = 2*Cf for laminar, turbulent and mixed flow
%
%   Reference: Anderson, Fundamentals of Aerodynamics, 6th ed., Sec. 4.12
%              (Eqs. 4.86, 4.88, 4.89, 4.92; Examples 4.8-4.11).

if nargin < 2, ReCrit = 5e5; end

CfLam  = @(R) 1.328 ./ sqrt(R);      % laminar,   Eq. 4.86
CfTurb = @(R) 0.074 ./ R.^0.2;       % turbulent, Eq. 4.88

% Mixed: laminar up to x_cr, turbulent behind it (Eq. 4.92).
% x_cr/c = ReCrit/Re (Eq. 4.89). If Re < ReCrit, all laminar.

ReX = min(ReCrit, Re);               % Reynolds at transition point
xcr = ReX ./ Re;                     % x_cr / c
CfM = CfTurb(Re) - xcr.*CfTurb(ReX) + xcr.*CfLam(ReX);

% Two faces wetted -> "net Cf" = 2*Cf = cd

cdL = 2 * CfLam(Re);
cdT = 2 * CfTurb(Re);
cdM = 2 * CfM;

end