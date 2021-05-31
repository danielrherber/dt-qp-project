%--------------------------------------------------------------------------
% DSuspensionNested_RoughStateConstraints.m
% State-dependent constraints for rough road load case
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function setup = DSuspensionNested_RoughStateConstraints(xd,param,setup)

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

% initialize index
idx = 1;

%--------------------------------------------------------------------------
% constraint to ensure spring linearity (gi2)
% In DTCD: g(1) = 0.15+1-(L0-Ls)/(delta_g+x(3)*1.1);
% NOTE: make UB/LB
ks = d^4*G/(8*D^3*Na*(1 + d^2/(2*D^2) ));
delta_g = ms*g/ks; % dg = ms*g/(4*ks);
L0 = p*Na + 2*d;
Ls = d*(Na+Q-1);
b1 = (L0-Ls-1.15*delta_g)/(1.1*1.15);
A1 = [0,0,1,0];

% Z(idx).linear(1).right = 2;
% Z(idx).linear(1).matrix = A1;
% Z(idx).b = b1;
% idx = idx+1;

% additional condition (delta_g+x(3)*1.1 >= 0)
A21 = [0,0,-1.1/delta_g,0];
% Z(idx).linear(1).right = 2;
% Z(idx).linear(1).matrix = A2;
% Z(idx).b = 1;
% idx = idx+1;
% not convinced the additional condition is correct
% Z(idx).linear(1).right = 2;
% Z(idx).linear(1).matrix = -A1;
% Z(idx).b = b1;
% idx = idx+1;

% validation
% x3 = 2*(rand(10000,1)-0.5);
% C = D/d;
% L0 = p*Na + 2*d;
% Ls = d*(Na + 1.75 - 1);
% ks = d^4 * G/ (8* D^3*Na*(1+1/2/C^2));
% delta_g = param.ms*9.81/ks;
% gDTCD = 0.15+1-(L0-Ls)./(delta_g+x3*1.1);
% gDTQP1 = A1(3)*x3-b1;
% gDTQP2 = A2(3)*x3-1;
% CDTCD = gDTCD<=0;
% CDTQP = (gDTQP1<=0) & (gDTQP2<=0); % combine the conditions
% % norm(CDTQP-CDTCD,inf) % should be the same (within FP errors)
% % disp(sum(CDTQP-CDTCD))
% disp(isequal(CDTQP,CDTCD))

%--------------------------------------------------------------------------
% Soderberg fatigue criterion (gi3)
% In DTCD: g(2) = tau_a/Se2 + tau_m/Ssy - 1;
% NOTE: make UB/LB
A2 = (G*d*nd*(1000*d)^m*((4*D)/d + 2))/(240000*A*D^2*Na*pi*(d^2/(2*D^2) + 1)*((4*D)/d - 3));
b2 = 1-(981*D*ms*(1000*d)^m*((4*D)/d + 2))/(8125000*A*d^3*pi*((4*D)/d - 3));

% Z(idx).linear(1).right = 2;
% Z(idx).linear(1).matrix = [0,0,A2,0];
% Z(idx).b = b2;
% idx = idx+1;

% validation
% x3 = 2*(rand(1000,1)-0.5);
% C = D/d;
% ks = d^4 * G/ (8* D^3*Na*(1+1/2/C^2));
% KB = (4*C+2)/(4*C-3);   % Bergstrausser factor
% delta_g = param.ms*9.81/ks;
% d_bar = 1000*d;
% Sut = A*1e6/(d_bar^m);
% Ssy = 0.65*Sut;
% Se2 = 0.24*Sut/nd;      % shear endurance strength (Stoicescu)
% Fmax = ks * (x3+ delta_g);        % maxium spring axial force
% Fmin = ks * (delta_g - x3);       % minium spring axial force
% Fm = (Fmax + Fmin)/2;               % mean spring axial force
% Fa = (Fmax - Fmin)/2;               % spring axial force amplitude
% tau_m = KB*8*Fm*D/(pi*d^3);
% tau_a = KB*8*Fa*D/(pi*d^3);
% gDTCD = tau_a/Se2 + tau_m/Ssy -1;    % Soderberg fatigue criterion
% gDTQP = A2*x3-b2;
% norm(gDTCD-gDTQP,inf) % should be the same (within FP errors)

% symbolic derivation
% syms D d G Na ms A m y3 nd real
% C = D/d;
% ks = d^4 * G/ (8* D^3*Na*(1+1/2/C^2));
% KB = (4*C+2)/(4*C-3);   % Bergstrausser factor
% delta_g = ms*9.81/ks;
% d_bar = 1000*d;
% Sut = A*1e6/(d_bar^m);
% Ssy = 0.65*Sut;
% Se2 = 0.24*Sut/nd;      % shear endurance strength (Stoicescu)
% Fmax = ks * (y3+ delta_g);        % maxium spring axial force
% Fmin = ks * (delta_g - y3);       % minium spring axial force
% Fm = (Fmax + Fmin)/2;               % mean spring axial force
% Fa = (Fmax - Fmin)/2;               % spring axial force amplitude
% tau_m = KB*8*Fm*D/(pi*d^3);
% tau_a = KB*8*Fa*D/(pi*d^3);
% [A,b] = equationsToMatrix(tau_a/Se2 + tau_m/Ssy,y3)

%--------------------------------------------------------------------------
% Zimmerli fatigue criterion (gi4)
% In DTCD: g(3) = (tau_a*nd-241e+6)/241e+6;
% NOTE: make UB/LB
A3 = zeros(1,4);
A3(3) = (3*G*d*((4*D)/d + 2))/(602500000*D^2*Na*pi*(d^2/(2*D^2) + 1)*((4*D)/d - 3));
b3 = 1;

% Z(idx).linear(1).right = 2;
% Z(idx).linear(1).matrix = A3;
% Z(idx).b = b3;
% idx = idx+1;

% validation (using DTCD implementation)
% x = 2*(rand(10,1)-0.5);
% C = D/d;
% ks = d^4 * G/ (8* D^3*Na*(1+1/2/C^2));
% delta_g = ms*9.81/ks;
% Fmax = ks * (x + delta_g);
% Fmin = ks * (delta_g - x);
% Fa = (Fmax - Fmin)/2;
% KB = (4*C+2)/(4*C-3);
% tau_a = KB*8*Fa*D/(pi*d^3);
% gDTCD = (tau_a*nd-241e+6)/241e+6;
% gDTQP = A3(3)*x-1;
% norm(gDTCD-gDTQP,inf) % should be the same (within FP errors)

% symbolic derivation
% syms ks dg y3 D d F G Na g ms
% ks = d^4*G/(8*D^3*Na*(1 + d^2/(2*D^2) ));
% dg = ms*g/ks;
% C = D/d;
% Fmax = ks*(y3+dg);
% Fmin = ks*(dg-y3);
% Fa = (Fmax-Fmin)/2;
% tau(F) = (4*C+2)/(4*C-3)*(8*F*D)/(pi*d^3);
% ex = 1.2*tau(Fa)/(241e6) - 1 == 0
% [A,b] = equationsToMatrix(ex,y3)

%--------------------------------------------------------------------------
% create single UB and LB constraints
[Mb,Ib] = min([b1,b2/A2,b3/A3(3)]);

setup.UB(end+1).right = 2;
setup.UB(end).matrix = [inf;inf;Mb;inf];
setup.LB(end+1).right = 2;
setup.LB(end).matrix = [-inf;-inf;A21(3);-inf];

%--------------------------------------------------------------------------
% damper pressure constraint (gi5)
% In DTCD: Pmax = g(4) = (Pmax - Pallow)/Pallow;
a4 = sqrt(rho1/(2*Pallow^3))*kv*Dp^2/(pi*eta*Afa*Cd*Do^3); % same value as above
A4 = [0,-a4,0,a4];

% Z(idx).linear(1).right = 2;
% Z(idx).linear(1).matrix = A4;
% Z(idx).b = 1;
% idx = idx+1;
%
% Z(idx).linear(1).right = 2;
% Z(idx).linear(1).matrix = -A4;
% Z(idx).b = 1;
% idx = idx+1;

% validation
% x2 = 2*(rand(10,1)-0.5); x4 = 2*(rand(10,1)-0.5);
% Ao = pi*Do^2/4;         % valve pressure surface area (M^2)
% xm = Ao*Pallow/kv;      % lift @Pallow (m)
% C2 = xi*Afa*sqrt(xm);   % damper valve parameter (m^0.5)
% KvC2 = sqrt(kv)/C2;     % valve parameter: sqrt(kv)/C2 (n^0.5/m)
% cs = Dp^4*KvC2*sqrt(pi*rho1/2)/(8*Cd*Do^2);   % damping ratio coefficient cs
% Pmax = 4*cs*(-x2+x4)/(pi*Dp^2); % maximum damper pressure
% gDTCD = (Pmax - Pallow)/Pallow;      % damper pressure constraint
% gDTQP = A4(2)*x2 + A4(4)*x4 - 1;
% norm(gDTCD-gDTQP,inf) % should be the same (within FP errors)
% gDTCD = (-Pmax - Pallow)/Pallow;
% gDTQP = -A4(2)*x2 + -A4(4)*x4 - 1;
% norm(gDTCD-gDTQP,inf) % should be the same (within FP errors)

%--------------------------------------------------------------------------
% damper velocity constraint (gi6)
% In DTCD: g(6) = -x(2) + x(4) - x3d_allow;
A6 = [0,-1/x3d_allow,0,1/x3d_allow]; % verified correct

% Z(idx).linear(1).right = 2;
% Z(idx).linear(1).matrix = A6;
% Z(idx).b = 1;
% idx = idx+1;
%
% Z(idx).linear(1).right = 2;
% Z(idx).linear(1).matrix = -A6;
% Z(idx).b = 1;
% idx = idx+1;

% validation
% x2 = 2*(rand(10,1)-0.5); x4 = 2*(rand(10,1)-0.5);
% gDTCD = -x2 + x4 - x3d_allow;
% gDTQP = (A6(2)*x2 + A6(4)*x4 - 1)*x3d_allow;
% norm(gDTCD-gDTQP,inf) % should be the same (within FP errors)
% gDTCD = x2 - x4 - x3d_allow;     % damper velocity constraint
% gDTQP = (-A6(2)*x2 + -A6(4)*x4 - 1)*x3d_allow;
% norm(gDTCD-gDTQP,inf) % should be the same (within FP errors)

%--------------------------------------------------------------------------
% damper valve lift constraint (gi7)
% In DTCD: g(8) = (xvmax - xv_allow)/xv_allow;
a8 = sqrt(rho1/(32*Pallow))*Dp^2/(eta*Afa*Cd*xv_allow*Do); % same value as above
A8 = [0,-a8,0,a8];

% Z(idx).linear(1).right = 2;
% Z(idx).linear(1).matrix = A8;
% Z(idx).b = 1;
% idx = idx+1;
%
% Z(idx).linear(1).right = 2;
% Z(idx).linear(1).matrix = -A8;
% Z(idx).b = 1;

% validation
% x2 = 2*(rand(10,1)-0.5); x4 = 2*(rand(10,1)-0.5);
% Ao = pi*Do^2/4;         % valve pressure surface area (M^2)
% xm = Ao*Pallow/kv;      % lift @Pallow (m)
% C2 = xi*Afa*sqrt(xm);   % damper valve parameter (m^0.5)
% KvC2 = sqrt(kv)/C2;     % valve parameter: sqrt(kv)/C2 (n^0.5/m)
% cs = Dp^4*KvC2*sqrt(pi*rho1/2)/(8*Cd*Do^2);   % damping ratio coefficient cs
% Pmax = 4*cs*(-x2+x4)/(pi*Dp^2); % maximum damper pressure
% xvmax = Ao*Pmax/kv;
% gDTCD = (xvmax - xv_allow)/xv_allow; % damper valve lift constraint
% gDTQP = A8(2)*x2 + A8(4)*x4 - 1;
% norm(gDTCD-gDTQP,inf) % should be the same (within FP errors)
% gDTCD = (-xvmax - xv_allow)/xv_allow;
% gDTQP = -A8(2)*x2 + -A8(4)*x4 - 1;
% norm(gDTCD-gDTQP,inf) % should be the same (within FP errors)

% only add the required x2 - x4 constraint
[Ms,Is] = min([A4(2),A6(2),A8(2)]);
As = [0,-Ms,0,Ms];
% disp(Is)

Z(idx).linear(1).right = 2;
Z(idx).linear(1).matrix = As;
Z(idx).b = 1;
idx = idx+1;

Z(idx).linear(1).right = 2;
Z(idx).linear(1).matrix = -As;
Z(idx).b = 1;

% assign
setup.Z = Z;

end