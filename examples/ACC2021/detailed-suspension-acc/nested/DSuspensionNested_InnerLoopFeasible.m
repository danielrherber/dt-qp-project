%--------------------------------------------------------------------------
% DSuspensionNested_InnerLoopFeasible.m
% Simply determine if the inner-loop problem is feasible
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function feasible = DSuspensionNested_InnerLoopFeasible(xd,param)

% initialize as feasible
feasible = true;

% extract plant variables
d = xd(1);  % spring wire diameter
D = xd(2);  % spring helix diameter
p = xd(3);  % spring pitch
Na = xd(4); % number of active spring coils
Do = xd(5); % orifice diameter
Dp = xd(6); % piston diameter

% extract problem parameters
ms = param.ms; g = param.g; G = param.G; Q = param.Q; rho1 = param.rho1;
kv = param.kv; Afa = param.Afa; eta = param.eta; Cd = param.Cd;
Pallow = param.Pallow; x3d_allow = param.x3d_allow;
xv_allow = param.xv_allow; nd = param.nd; m = param.m; A = param.A;
xi = param.xi;

% constraint to ensure spring linearity (gi2)
ks = d^4*G/(8*D^3*Na*(1 + d^2/(2*D^2) ));
delta_g = ms*g/ks;
L0 = p*Na + 2*d;
Ls = d*(Na+Q-1);
b1 = (L0-Ls-1.15*delta_g)/(1.1*1.15);
if b1 <= 0
    feasible = false;
    return
end

% Soderberg fatigue criterion (gi3)
A2 = (G*d*nd*(1000*d)^m*((4*D)/d + 2))/(240000*A*D^2*Na*pi*(d^2/(2*D^2) + 1)*((4*D)/d - 3));
if A2 <= 0
    disp(1)
end
b2 = 1-(981*D*ms*(1000*d)^m*((4*D)/d + 2))/(8125000*A*d^3*pi*((4*D)/d - 3));
if b2/A2 <= 0
    feasible = false;
    return
end

% Zimmerli fatigue criterion (gi4)
A3 = zeros(1,4);
A3(3) = (3*G*d*((4*D)/d + 2))/(602500000*D^2*Na*pi*(d^2/(2*D^2) + 1)*((4*D)/d - 3));
if A3(3) <= 0
    feasible = false;
    return
end

% rattlespace constraint (gi1)
b1 = -(Na*p+2*d)+(d*(Na+Q-1))+0.02+(4*ms*g*D*Na*(2*D^2+d^2)/(G*d^4)); % same value as above
if -b1 <= 0
    feasible = false;
    return
end

end