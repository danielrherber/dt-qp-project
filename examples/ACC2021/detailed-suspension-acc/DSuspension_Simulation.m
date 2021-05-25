%--------------------------------------------------------------------------
% DSuspension_Simulation.m
% Run a simulation given a plant design for Detailed Suspension problem
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [T,U,Y] = DSuspension_Simulation(xd,x0,p)

% extract plant variables
d = xd(1);  % spring wire diameter
D = xd(2);  % spring helix diameter
Na = xd(4); % number of active spring coils
Do = xd(5); % orifice diameter
Dp = xd(6); % piston diameter

% extract problem parameters
G = p.G; rho1 = p.rho1; kv = p.kv; Cd = p.Cd;
Afa = p.Afa; Pallow = p.Pallow; xi = p.xi;

% intermediate parameters
C = D/d;                % Spring index
Ao = pi*Do^2/4;         % valve pressure surface area (M^2)
xm = Ao*Pallow/kv;      % lift @ Pallow (m)
C2 = xi*Afa*sqrt(xm);   % damper valve parameter (m^0.5)

% spring and damper constants
xp(1) = Dp^4*sqrt(kv* pi*rho1/2)/(8*Cd*Do^2*C2); % damper coefficient
xp(2) = d^4*G/(8* D^3*Na*(1+1/2/C^2)); % spring constant

% extract
mu = p.mu; ms = p.ms; bt = p.bt; kt = p.kt;
bs = xp(1); ks = xp(2);
z0dot = p.z0dot;
ramp_in = 25/100*10;

% state
A = [0,1,0,0;
    -kt/mu,-(bs+bt)/mu,ks/mu,bs/mu;
    0,-1,0,1;
    0,bs/ms,-ks/ms,-bs/ms];

% control
Bu = [0;-1/mu;0;1/ms];

% disturbance
Bz = [-1;bt/mu;0;0];

% duplicate matrices for two design load cases
A2 = blkdiag(A,A);
Bu2 = [Bu;Bu];
Bz2 = [Bz;Bz];

% derivative function
fd = @(t,y) A2*y + Bu2*0 + [Bz*z0dot(t);Bz*ramp_in];

% ode options
OPTIONS = odeset('RelTol',1e-8,'AbsTol',1e-10);

% run the simulation
[T,Y] = ode15s(@(t,y) fd(t,y),[0 p.tf],x0,OPTIONS);

% zero control
U = zeros(length(T),2);

end