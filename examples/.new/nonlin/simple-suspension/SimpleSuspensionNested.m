%--------------------------------------------------------------------------
% SimpleSuspensionNested.m
% D. R. Herber and A. K. Sundarrajan, "On the uses of linear-quadratic
% methods in solving nonlinear dynamic optimization problems with direct
% transcription", in ASME International Mechanical Engineering Congress &
% Exposition, 2020.
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clearvars -except externalInput

% auxiliary data
auxdata = SimpleSuspensionProblemParameters;
auxdata.t0 = 0; auxdata.tf = 3;

% limits for plant design
lb = [auxdata.bmin,auxdata.kmin];
ub = [auxdata.bmax,auxdata.kmax];

% solve the outer-loop, plant-only problem
t1 = tic;
[auxdata,~] = OuterLoop(auxdata,lb,ub);
toc(t1)

% final linear model update
[A,Bu,Bz,C1,C2,C3,D1u,D3u] = UpdateLinearModel(auxdata.xpopt,auxdata);

% final solution
[F,T,U,Y,P,in,opts] = InnerLoop(auxdata,A,Bu,Bz,C1,C2,C3,D1u,D3u);

% plots
SimpleSuspension_plot(T,U,Y,P,F,in,opts,[])

%--------------------------------------------------------------------------
% solve the outer-loop problem
%--------------------------------------------------------------------------
function [p,F] = OuterLoop(p,lb,ub,varargin)

% number of plant design variables
nxp = length(lb);

% ensure column vectors
lb = lb(:); ub = ub(:);

% scale
LB = Scaling(lb,lb,ub,1);
UB = Scaling(ub,lb,ub,1);

% initialize plant design
rng(8920543); % fix seed
X = LB + (UB-LB).*rand(nxp,1);

%--------------------------------------------------------------------------
% % patternsearch options
% options = psoptimset('Display','Iter','MaxIter',10,'Cache','on',...
%     'UseParallel',true,'CompletePoll','on','CompleteSearch','on',...
%     'SearchMethod', {@searchlhs,1,nxp*30});
%
% % global optimization
% X = patternsearch(@(x) Objective(x,p,lb,ub),X,[],[],[],[],LB,UB,[],options);

%--------------------------------------------------------------------------
% fmincon options
options = optimoptions('fmincon','Display','Iter','Algorithm','interior-point',...
    'UseParallel',true,'MaxIter',100,'OptimalityTolerance',1e-5,...
    'FiniteDifferenceStepSize',sqrt(eps),'FiniteDifferenceType','forward');

% local optimization
[X,F] = fmincon(@(x) Objective(x,p,lb,ub),X,[],[],[],[],LB,UB,[],options);

% unscale from linear 0 to 1 mapping
X = Scaling(X,lb,ub,2);

% assign unscaled plant variables
p.xpopt = X;

end

%--------------------------------------------------------------------------
% objective function for outer-loop problem
%--------------------------------------------------------------------------
function [F] = Objective(xp,p,lb,ub)

% unscale from linear 0 to 1 mapping
xp = Scaling(xp,lb,ub,2);

% update matrices
[A,Bu,Bz,C1,C2,C3,D1u,D3u] = UpdateLinearModel(xp,p);

% find inner-loop solution
F = InnerLoop(p,A,Bu,Bz,C1,C2,C3,D1u,D3u);

end

%--------------------------------------------------------------------------
% update the linear dynamic and output models
%--------------------------------------------------------------------------
function [A,Bu,Bz,C1,C2,C3,D1u,D3u] = UpdateLinearModel(xp,p)

% extract plant design
bs = xp(1); ks = xp(2);

% state dynamics
A = [0,1,0,0;
    -p.kt/p.mu,-(bs+p.bt)/p.mu,ks/p.mu,bs/p.mu;
    0,-1,0,1;
    0,bs/p.ms,-ks/p.ms,-bs/p.ms];
Bu = [0;-1/p.mu;0;1/p.ms];
Bz = [-1;p.bt/p.mu;0;0];

% output equation
C1 = [1,0,0,0];
C2 = [0,0,1,0];
C3 = A(4,:);
D1u = Bu(1,:);
D3u = Bu(4,:);

end

%--------------------------------------------------------------------------
% construct and solve the inner-loop problem
%--------------------------------------------------------------------------
function [F,varargout] = InnerLoop(p,A,Bu,Bz,C1,C2,C3,D1u,D3u)

% initialize struct
setup = DTQP_setup_initialize;

% extract
w1 = p.w1; w2 = p.w2; w3 = p.w3; % objective weights

% number of states
ns = 4;

% disturbance
d = cell(ns,1);
for idx = 1:ns
    d{idx,1} = @(t,p) Bz(idx)*p.z0dot(t)';
end

% Lagrange terms
L(1).left = 2; % states
L(1).right = 2; % states
L(1).matrix = w1*(C1'*C1) + w2*(C3'*C3);
L(2).left = 1; % controls
L(2).right = 1; % controls
L(2).matrix = w1*(D1u'*D1u) + w2*(D3u'*D3u) + w3;
L(3).left = 1; % controls
L(3).right = 2; % states
L(3).matrix = 2*w1*(D1u'*C1) + 2*w2*(D3u'*C3);

% initial state values
LB(1).right = 4; % initial states
LB(1).matrix = zeros(ns,1);
UB(1).right = 4; % initial states
UB(1).matrix = zeros(ns,1);

% rattlespace constraints
LB(2).right = 2; % states
LB(2).matrix = [-inf,-inf,-p.rmax,-inf];
UB(2).right = 2; % states
UB(2).matrix = [inf,inf,p.rmax,inf];

% combine
setup.auxdata = p;
setup.t0 = 0;
setup.tf = p.tf;
setup.lq.lagrange = L;
setup.lq.dynamics.A = A;
setup.lq.dynamics.B = Bu;
setup.lq.dynamics.Bv = d;
setup.lq.ub = UB;
setup.lq.lb = LB;

% DTQP options
opts = SimpleSuspensionNested_opts;

% form and solve the LQDO problem
[T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);

% optional additional outputs
if nargout > 1
    varargout = {T,U,Y,P,in,opts};
end

end

%--------------------------------------------------------------------------
% inner-loop, DTQP options
%--------------------------------------------------------------------------
function opts = SimpleSuspensionNested_opts

opts.dt.defects = 'TR';
opts.dt.quadrature = 'CTR';
opts.dt.mesh = 'ED';
opts.dt.nt = 1000;
opts.general.displevel = 1;
opts.solver.tolerance = 1e-16;

end

%--------------------------------------------------------------------------
% local scaling function
%--------------------------------------------------------------------------
function x = Scaling(x,l,u,type)

% reshape to column vector
x = reshape(x,[],1);

% scale or unscale
if type == 1 % scale

    % log scale
    % x = log10(x); l = log10(l); u = log10(u);

    % scale to linear 0 to 1 mapping
    x = (x-l)./(u-l);

elseif type == 2 % unscale

    % log scale
    % l = log10(l); u = log10(u);

    % unscale from linear 0 to 1 mapping
    x = l + x.*(u-l);

    % unscale from log mapping
    % x = 10.^(x);

end

end