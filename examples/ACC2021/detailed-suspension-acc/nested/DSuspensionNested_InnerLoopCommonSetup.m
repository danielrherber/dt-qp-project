%--------------------------------------------------------------------------
% DSuspensionNested_InnerLoopCommonSetup.m
% Common inner-loop problem setup
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [setup,Bz] = DSuspensionNested_InnerLoopCommonSetup(x,p)

% number of states
ns = 4;

% update linear dynamic model
[A,Bu,Bz] = UpdateLinearModel(x,p);
setup.A = A; setup.B = Bu;

% matrices for each objective term
C1 = [1,0,0,0]; D1u = Bu(1,:);
C3 = A(4,:); D3u = Bu(4,:);

% objective weights
w1 = p.w1;
w2 = p.w2;
w3 = p.w3;

% Lagrange term (DRH verified, note we use CTR, DTCD uses CEF)
L(1).left = 2; L(1).right = 2; % states
L(1).matrix = w1*(C1'*C1) + w2*(C3'*C3);
L(2).left = 1; L(2).right = 1; % controls
L(2).matrix = w1*(D1u'*D1u) + w2*(D3u'*D3u) + w3;
L(3).left = 1; L(3).right = 2; %  controls-states
L(3).matrix = 2*w1*(D1u'*C1) + 2*w2*(D3u'*C3);
setup.L = L;

% initial value constraints
LB(1).right = 4; % initial states
LB(1).matrix = zeros(ns,1);
UB(1).right = 4; % initial states
UB(1).matrix = zeros(ns,1);
setup.UB = UB; setup.LB = LB;

end

% update linear dynamic model
function [A,Bu,Bz] = UpdateLinearModel(x,p)

% calculate intermediate plant variables
xp = UpdateIntermediatePlant(x,p);

% extract
mu = p.mu; ms = p.ms; bt = p.bt; kt = p.kt;
bs = xp(1); ks = xp(2);

% state
A = [0,1,0,0;
    -kt/mu,-(bs+bt)/mu,ks/mu,bs/mu;
    0,-1,0,1;
    0,bs/ms,-ks/ms,-bs/ms];

% control
Bu = [0;-1/mu;0;1/ms];

% disturbance
Bz = [-1;bt/mu;0;0];

end

% spring and damper constants
function xp = UpdateIntermediatePlant(xd,param)

% extract plant variables
d = xd(1);  % spring wire diameter
D = xd(2);  % spring helix diameter
Na = xd(4); % number of active spring coils
Do = xd(5); % orifice diameter
Dp = xd(6); % piston diameter

% extract problem parameters
G = param.G; rho1 = param.rho1; kv = param.kv; Cd = param.Cd;
Afa = param.Afa; Pallow = param.Pallow; xi = param.xi;

% intermediate parameters
C = D/d;                % Spring index
Ao = pi*Do^2/4;         % valve pressure surface area (M^2)
xm = Ao*Pallow/kv;      % lift @ Pallow (m)
C2 = xi*Afa*sqrt(xm);   % damper valve parameter (m^0.5)

% spring and damper constants
xp(1) = Dp^4*sqrt(kv* pi*rho1/2)/(8*Cd*Do^2*C2); % NEED COMMENT
xp(2) = d^4*G/(8* D^3*Na*(1+1/2/C^2)); % NEED COMMENT

end