%--------------------------------------------------------------------------
% DSuspensionNested_DesignConstraints.m
% (Nonlinear) plant design constraints
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [c,ceq] = DSuspensionNested_DesignConstraints(xd,param)

% extract plant variables
d = xd(1);  % spring wire diameter
D = xd(2);  % spring helix diameter
p = xd(3);  % spring pitch
Na = xd(4); % number of active spring coils
Dp = xd(6); % piston diameter
Ds = xd(7); % damper stroke

% extract problem parameters
G = param.G; L0max = param.L0max; DoMax = param.DoMax; nd = param.nd;
A = param.A; m = param.m; dsc = param.dsc; dwt = param.dwt;
ld1 = param.ld1; ld2 = param.ld2;

% intermediate parameters
C = D/d;
L0 = p*Na + 2*d;
Ls = d*(Na + 1.75 - 1);
ks = d^4*G/(8*D^3*Na*(1+1/2/C^2));
Fs = ks*(L0 - Ls);
KB = (4*C+2)/(4*C-3);
tau_s = KB*8*Fs*D/(pi*d^3);
d_bar = 1000*d;
Sut = A*1e6/(d_bar^m);
Ssy = 0.65*Sut;

%--- inequality constraints
% initialize
% c = zeros(9,1);
c = zeros(4,1);

% spring manufacturing constraint
% In DTCD: g(1) = 4 - C; (DRH checked)
% NOTE: make linear
% c(1) = (4*d - D)/DoMax;

% spring manufacturing constraint
% In DTCD: g(2) = C - 12; (DRH checked)
% NOTE: make linear
% c(2) = (D - 12*d)/DoMax;

% spring stability/buckling constraint for squared ground ends
% In DTCD: g(3) = L0 - 5.26*D; (DRH checked)
% c(3) = p*Na/(5.26*DoMax) + 2*d/(5.26*DoMax) - 5.26*D/(5.26*DoMax);
c(1) = p*Na/(5.26*DoMax) + 2*d/(5.26*DoMax) - 5.26*D/(5.26*DoMax);

% spring packaging constraint
% In DTCD: g(4) = L0 - L0max; (DRH checked)
% c(4) = (p*Na+2*d)/L0max - 1;
c(2) = (p*Na+2*d)/L0max - 1;

% spring packaging constraint
% In DTCD: g(5) = d + D - DoMax; (DRH checked)
% NOTE: make linear
% c(5) = d/DoMax + D/DoMax - 1;

% spring-damper interference constraint
% suspension stop constraint (ensures don't hit stops during test bump)
% In DTCD: g(6) = d - D + Dp + 2*(dsc+dwt); (DRH checked)
% NOTE: make linear
% c(6) = d/(2*(dsc+dwt)) - D/(2*(dsc+dwt)) + Dp/(2*(dsc+dwt)) + 1;

% shear stress constraint
% In DTCD: g(7) = (tau_s*nd - Ssy)/Ssy; (DRH checked)
% c(7) = (tau_s*nd/Ssy - 1);
c(3) = (tau_s*nd/Ssy - 1);

% ensures adequate damper stroke
% In DTCD: g(8) = L0 - Ls - Ds; (DRH checked)
% c(8) =  p*Na/L0max + (2*d - d*(Na + 1.75 - 1))/L0max - Ds/L0max;
c(4) =  p*Na/L0max + (2*d - d*(Na + 1.75 - 1))/L0max - Ds/L0max;

% ensures enough space for damper
% In DTCD: g(9) = 2*Ds + ld1 + ld2 - L0max;
% NOTE: make linear
% c(9) = 2*Ds + ld1 + ld2 - L0max;

% c = []; % remove all (for debugging)

%--- equality constraints
ceq = [];

end