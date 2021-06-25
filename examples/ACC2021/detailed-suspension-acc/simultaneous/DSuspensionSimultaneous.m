%--------------------------------------------------------------------------
% DSuspensionSimultaneous.m
% Simultaneous implementation of the control co-design problem in:
% J. T. Allison, T. Guo, and Z. Han, "Co-Design of an Active Suspension
% Using Simultaneous Dynamic Optimization," Journal of Mechanical Design,
% vol. 136, no. 8, Jun. 2014, doi: 10.1115/1.4027335
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function sol = DSuspensionSimultaneous(varargin)

% check if inputs were provided
if isempty(varargin)
    p.derivatives = 'symbolic'; % derivative method
    p.nt = 200; % number of time points
    p.OptimalityTolerance = 1e-7; % optimality tolerance
    p.FeasibilityTolerance = 1e-12; % feasibility tolerance
else
    in = varargin;
    p.derivatives = in{1}; % derivative method
    p.nt = in{2}; % number of time points
    p.OptimalityTolerance = in{3}; % optimality tolerance
    p.FeasibilityTolerance = in{4}; % feasibility tolerance
end

% get DTQP default options for this problem
opts = DSuspensionSimultaneous_opts;

% assign options
opts.method.derivatives = p.derivatives;
opts.dt.nt = p.nt;
opts.solver.Otolerance = p.OptimalityTolerance;
opts.solver.Ftolerance = p.FeasibilityTolerance;

% opts for sens study
if ~isempty(varargin)
    opts.general.displevel = 1;
end

% problem parameters
p = DSuspensionProblem_Parameters;

% time horizon
p.t0 = 0; p.tf = 2;

%% setup
% number of controls, states, and parameters
n.nu = 2; n.ny = 8; n.np = 7;

%--- system dynamics
% rough road input
str{1} =  '[';
str{end+1} = 'y2 - z0d;';
str{end+1} = '(- kt*y1 - ((281.5194*p6^4/(p5^3)) + ct)*y2 + (G*p1^4/(4*p2*p4*(2*p2^2+p1^2)))*y3 + (281.5194*p6^4/(p5^3))*y4 + ct*z0d - u1)/mus;';
str{end+1} = 'y4 - y2;';
str{end+1} = '((281.5194*p6^4/(p5^3))*y2 - (G*p1^4/(4*p2*p4*(2*p2^2+p1^2)))*y3 - (281.5194*p6^4/(p5^3))*y4 + u1)/ms;';

% ramp road input
str{end+1} = 'y6 - R;';
str{end+1} = '(- kt*y5 - ((281.5194*p6^4/(p5^3)) + ct)*y6 + (G*p1^4/(4*p2*p4*(2*p2^2+p1^2)))*y7 + (281.5194*p6^4/(p5^3))*y8 + ct*R - u2)/mus;';
str{end+1} = 'y8 - y6;';
str{end+1} = '((281.5194*p6^4/(p5^3))*y6 - (G*p1^4/(4*p2*p4*(2*p2^2+p1^2)))*y7 - (281.5194*p6^4/(p5^3))*y8 + u2)/ms';
str{end+1} = ']';

% combine
element.dynamics = horzcat(str{:});

%--- objective function
% rough road input
str0{1} = 'w1*y1^2 + ';
str0{2} = 'w2/ms^2*((281.5194*p6^4/(p5^3))*y2 - (G*p1^4/(4*p2*p4*(2*p2^2+p1^2)))*y3 - (281.5194*p6^4/(p5^3))*y4 + u1)^2 + ';
str0{3} = 'w3*u1^2+';

% ramp road objective
str0{4} = '(w1*y5^2 + ';
str0{5} = 'w2/ms^2*((281.5194*p6^4/(p5^3))*y6 - (G*p1^4/(4*p2*p4*(2*p2^2+p1^2)))*y7 - (281.5194*p6^4/(p5^3))*y8 + u2)^2 + ';
str0{6} = 'w3*u2^2)/100';

% combine
element.lagrange = horzcat(str0{:});

%--- symbolic parameters
element.parameter_list = 'ct kt mus ms w1 w2 w3 G z0d R m A nd';
element.parameter_values = {p.bt p.kt p.mu p.ms p.w1 p.w2 p.w3 p.G p.z0dot p.ramp_in 0.108 1974 1.2 };

%--- initial state values
LB(1).right = 4;
LB(1).matrix = zeros(8,1);
UB(1).right = 4;
UB(1).matrix = zeros(8,1);

%--- linear inequality constraints
% initialize
idz = 0;

% manufacturability g_{o,1}
idz = idz + 1; % increment
Z(idz).linear(1).right = 3;
Z(idz).linear(1).matrix = [4/0.25,-1/0.25,0,0,0,0,0]';
Z(idz).b = 0;

% tangling constraint g_{o,2}
idz = idz + 1; % increment
Z(idz).linear(1).right = 3; % parameters
Z(idz).linear(1).matrix = [-12/0.25,1/0.25,0,0,0,0,0]';
Z(idz).b = 0;

% interference constraint g_{o,5}
idz = idz + 1; % increment
Z(idz).linear(1).right = 3; % parameters
Z(idz).linear(1).matrix = [1/0.25,1/0.25,0,0,0,0,0]';
Z(idz).b = 1;

% damper size constraint g_{o,7}
idz = idz + 1; % increment
Z(idz).linear(1).right = 3; % parameters
Z(idz).linear(1).matrix = [1/0.022,-1/0.022,0,0,0,1/0.022,0]';
Z(idz).b = -1;

% damper range of motion g_{o,8}
idz = idz + 1; % increment
Z(idz).linear(1).right = 3; % parameters
Z(idz).linear(1).matrix = [0,0,0,0,0,0,2/0.34]';
Z(idz).b = 1;

% maximum velocity constraints g_{i,6}, ramp
idz = idz + 1; % increment
Z(idz).linear(1).right = 2; % states
Z(idz).linear(1).matrix = [0,-1/5,0,1/5,0,0,0,0]';
Z(idz).b = 1;

idz = idz + 1; % increment
Z(idz).linear(1).right = 2; % states
Z(idz).linear(1).matrix = [0,1/5,0,-1/5,0,0,0,0]';
Z(idz).b = 1;

% maximum velocity constraints g_{i,6}, rough road
idz = idz + 1; % increment
Z(idz).linear(1).right = 2; % states
Z(idz).linear(1).matrix = [0,0,0,0,0,-1/5,0,1/5]';
Z(idz).b = 1;

idz = idz + 1; % increment
Z(idz).linear(1).right = 2; % states
Z(idz).linear(1).matrix = [0,0,0,0,0,1/5,0,-1/5]';
Z(idz).b = 1;

%% rough road input
% initialize
idx = 1;
strC{idx} = '[';

% buckling constraint g_{o,3}
idx = idx + 1; % increment
strC{idx} = '(p3*p4+2*p1)/(5.26*.25)-(p2/0.25);';

% pocket length constraint g_{o,4}
idx = idx + 1; % increment
strC{idx} = '(p3*p4/0.4)+(2*p1/0.4)-1;';

% damper range of motion g_{o,9}
idx = idx + 1; % increment
strC{idx} = '(p3*p4/0.4)+(1.25*p1/0.4)-(p1*p4/0.4)-(p7/0.4);';

% scaled stress constraint g_{o,6}
idx = idx + 1; % increment
strC{idx} = '309.2441*p1^(4.108)*(4*p2+2*p1)*(p3*p4+1.25*p1-p1*p4)/(p2*p4*(2*p2^2+p1^2)*(4*p2-3*p1))-1;';

% spring linearity constraint g_{i,2}
idx = idx + 1; % increment
strC{idx} = 'y3 - (400*p1)/253 - (200*p4*p3)/253 + (200*p1*(p4 + 3/4))/253 + (3924*p2^3*p4*ms*(p1^2/(2*p2^2) + 1))/(55*G*p1^4);';

% Soderberg constraint g_{i,3}
idx = idx + 1; % increment
strC{idx} = '(9.81*p2*ms*(1000*p1)^m*((4*p2)/p1 + 2))/(81250*A*p1^3*pi*((4*p2)/p1 - 3)) + (G*p1*nd*y3*(1000*p1)^m*((4*p2)/p1 + 2))/(240000*A*p2^2*p4*pi*(p1^2/(2*p2^2) + 1)*((4*p2)/p1 - 3)) - 1;';

% Zimmeli constraint g_{i,4}
idx = idx + 1; % increment
strC{idx} = '(3*G*p1*y3*((4*p2)/p1 + 2))/(602500000*p2^2*p4*pi*(p1^2/(2*p2^2) + 1)*((4*p2)/p1 - 3)) - 1 ;';

% damper pressure constraints g_{i,5}
idx = idx + 1; % increment
strC{idx} = '7.5461*10^-5*(p6^2/p5^3)*(y4-y2)-1;';

idx = idx + 1; % increment
strC{idx} = '7.5461*10^-5*(p6^2/p5^3)*(-y4+y2)-1;';

% valve lift constraints g_{i,7}
idx = idx + 1; % increment
strC{idx} = '1.2512*(p6^2/p5)*(y4-y2)-1;';

idx = idx + 1; % increment
strC{idx} = '1.2512*(p6^2/p5)*(-y4+y2)-1;';

%% ramp input
% rattlespace constraint g_{i,1}
idx = idx + 1; % increment
strC{idx} = 'y7 - 2*p1 - p4*p3 + p1*(p4 + 3/4) + (9.81*4*p2*p4*ms*(2*p2^2 + p1^2))/(G*p1^4) + 1/50;';

% damper pressure constraints g_{i,5}
idx = idx + 1; % increment
strC{idx} = '7.5461*10^-5*(p6^2/p5^3)*(y8-y6)-1;';

idx = idx + 1; % increment
strC{idx} = '7.5461*10^-5*(p6^2/p5^3)*(-y8+y6)-1;';

% valve lift constraints g_{i,7}
idx = idx + 1; % increment
strC{idx} = '1.2512*(p6^2/p5)*(y8-y6)-1;';

idx = idx + 1; % increment
strC{idx} = '1.2512*(p6^2/p5)*(-y8+y6)-1;';

% combine
idx = idx + 1; % increment
strC{idx} = ']';
element.g.func = horzcat(strC{:});
element.g.pathboundary = [0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1];

% simple parameter bounds
LB(2).right = 3;
LB(2).matrix = [0.005 0.05 0.02 3 0.003 0.03 0.1];
UB(2).right = 3;
UB(2).matrix = [0.02 0.4 0.5 16 0.012 0.08 0.3];

% get initial guess
setup = DSuspensionSimultaneous_guess([]);

% combine structures
setup.element = element; setup.UB = UB; setup.LB = LB;setup.Z = Z;
setup.t0 = p.t0; setup.tf = p.tf; setup.p = p; setup.n = n;

%% solve
t1 = tic;
[T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);
tc = toc(t1);

%% output
% assign output structure
sol.T = T; sol.U = U; sol.Y = Y; sol.P = P; sol.F = F;
sol.in = in; sol.opts = opts;
sol.tc = tc - opts.timer.sym;
sol.xp = P;

%% plot
if isempty(varargin)
% assign ramp and rough road states and controls
Yrough = Y(:,1:4); Yramp = Y(:,5:end);
Urough = U(:,1); Uramp = U(:,2);

% plot ramp solution
DSuspension_plot(T,Uramp,Yramp,P,F,in,opts,sol)

% plot rough road solution
DSuspension_plot(T,Urough,Yrough,P,F,in,opts,sol)
end
end

% initial guess
function setup = DSuspensionSimultaneous_guess(setup)

% run a simulation with the initial plant design
simflag = false;

% initial plant number (see below)
InitialPlant = 2;

% set initial plant design
switch InitialPlant
    %----------------------------------------------------------------------
    case 1 % initial values from Allison2014b
    p0 = [0.01,0.12,0.05,6,0.0067,0.04,0.15];
    %----------------------------------------------------------------------
    case 2 % initial values from Allison2014b
    p0 = [0.01,0.129,0.106,3.57,0.006,0.035,0.17];
    %----------------------------------------------------------------------
    case 3 % (lb+ub)/2
    p0 = [0.0125,0.2250,0.2600,9.5000,0.0075,0.0550,0.2000];
    %----------------------------------------------------------------------
    case 4 % optimal values from Allison2014b
    p0 = [0.0097,0.0620,0.0201,15.3,0.0061,0.0303,0.170];
    %----------------------------------------------------------------------
    case 5 % optimal values using simultaneous
    p0 = [0.0200 0.1820 0.0335 10.7308 0.0097 0.0405 0.1700];
    %----------------------------------------------------------------------
    case 6 % optimal values using nested
    p0 = [0.0166 0.1686 0.0396 6.5015 0.0094 0.0395 0.1700];
    %----------------------------------------------------------------------
end

% determine if a simulation needs to be run
if simflag

    % initial states
    x0 = zeros(8,1);

    % run the simulation
    [T0,U0,Y0] = DSuspension_Simulation(p0,x0,p);

    % replicate constant plant design
    P0 = repmat(p0,length(T0),1);

    % assign time field
    p.Tguess = T0;

else

    % zero states and controls
    Y0 = [[0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0]];
    U0 = [[0,0];[0,0]];

    % replicate constant plant design
    P0 = [p0;p0];

end

% combine
setup.guess.X = [U0,Y0,P0];

end

% default DTQP options
function opts = DSuspensionSimultaneous_opts

% options
opts.general.displevel = 2;
opts.general.plotflag = 1;
opts.dt.defects = 'TR';
opts.dt.quadrature = 'CTR';
opts.method.form = 'nonlinearprogram';
opts.method.derivatives = 'real-central';
opts.method.olqflag = true; % need to fix
opts.solver.tolerance = 1e-12;
opts.solver.maxiters = 5000;
opts.solver.function = 'IPFMINCON';

% mesh
opts.dt.nt = 200;
opts.dt.mesh = 'ED';
% opts.dt.mesh = 'USER';
% Allison2014b.rough = load('RoughMeshAllison2014b','t');
% opts.dt.t = Allison2014b.rough.t;

end