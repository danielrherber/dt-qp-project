%--------------------------------------------------------------------------
% SpaceShuttleReentry_helper.m
% Helper function for the symbolic equations in the Space Shuttle Reentry
% example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clear; clc;

% initialize
syms t
syms rad lon lat v fpa azi aoa bank
syms y1 y2 y3 y4 y5 y6 u1 u2
syms Re cd0 cd1 cd2 rho0 H cl0 cl1 S mass mu m2ft

% intermediate parameters (dynamics)
cbank = cos(bank);
sbank = sin(bank);
altitude = rad - Re;
CD = cd0 + cd1*aoa + cd2*aoa^2;
rho = rho0*exp(-altitude/H);
CL = cl0 + cl1*aoa;
q = 0.5*rho*v^2;
D = q*S*CD/mass;
L = q*S*CL/mass;
gravity  = mu/rad^2;

% intermediate parameters (heating)
vI = 3.280839895*v;
hI = 3.280839895*(rad - Re);
rho = 0.002378*exp(-hI/23800);
qr = 17700*sqrt(rho).*(0.0001*vI).^3.07;
a = rad2deg(aoa);
qa = 1.0672181 + -0.192137741e-1*a + 0.21286289e-3*a.^2 + -0.10117249e-5*a.^3;
q = qa.*qr;

% state derivatives
raddot = v*sin(fpa);
londot = v*cos(fpa)*sin(azi)/(rad.*cos(lat));
latdot = v*cos(fpa)*cos(azi)/rad;
vdot   = -D-gravity*sin(fpa);
fpadot = (L*cbank-cos(fpa)*(gravity-v^2/rad))/v;
azidot = (L*sbank/cos(fpa)+v^2*cos(fpa)*sin(azi)*tan(lat)/rad)/v;

% combine
eqs = [raddot;londot;latdot;vdot;fpadot;azidot;q];

% substitutions
eqs = subs(eqs,rad,y1);
eqs = subs(eqs,lon,y2);
eqs = subs(eqs,lat,y3);
eqs = subs(eqs,v,y4);
eqs = subs(eqs,fpa,y5);
eqs = subs(eqs,azi,y6);
eqs = subs(eqs,aoa,u1);
eqs = subs(eqs,bank,u2);

% display
disp(eqs)