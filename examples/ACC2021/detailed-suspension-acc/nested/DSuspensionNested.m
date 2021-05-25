%--------------------------------------------------------------------------
% DSuspensionNested.m
% Nested implementation of the control co-design problem in:
% J. T. Allison, T. Guo, and Z. Han, "Co-Design of an Active Suspension
% Using Simultaneous Dynamic Optimization,” Journal of Mechanical Design,
% vol. 136, no. 8, Jun. 2014, doi: 10.1115/1.4027335
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function sol = DSuspensionNested(varargin)

if ~isempty(varargin)
    p = varargin{1};
else
    p.GAflag = true;
    p.FiniteDifferenceType = 'central'; % derivative method
    p.nt = 100; % number of time points
    p.InnerLoopTolerance = 1e-14; % inner-loop tolerance
    p.OuterLoopTolerance = 1e-3; % outer-loop tolerance
end

% select mesh number for the inner-loop problems
p.meshnum = 1; % see below

% parameters
p = DSuspensionProblem_Parameters(p);

%% solve outer-loop problem
t1 = tic;
[xp,F] = OuterLoop(p);
tc = toc(t1);

%% output
[~,sol.ramp] = DSuspensionNested_InnerLoopRamp(xp,p);
[~,sol.rough] = DSuspensionNested_InnerLoopRough(xp,p);
sol.xp = xp;
sol.F = F;
sol.tc = tc;

%% plot
if isempty(varargin)
%close all open figures
close all

% extract
ramp = sol.ramp; rough = sol.rough;

% change plot flag
ramp.opts.general.plotflag = true;
rough.opts.general.plotflag = true;

% plot ramp solution
DSuspension_plot(ramp.T,ramp.U,ramp.Y,sol.xp,ramp.F,ramp.in,ramp.opts,[])

% plot rough road solution
DSuspension_plot(rough.T,rough.U,rough.Y,sol.xp,rough.F,rough.in,rough.opts,[])
end
end

% solve the outer-loop problem
function [xp,F] = OuterLoop(p)

% extract
UB = p.xpUpperBound; LB = p.xpLowerBound;

% number of plant design variables
nxp = length(LB);

% linear outer-loop (plant design) constraints
[A,b] = DSuspensionNested_LinearDesignConstraints(p);

%% initialize plant design

% initial plant number (see below)
InitialPlant = 7;

% set initial plant design
switch InitialPlant
    %----------------------------------------------------------------------
    case 1 % initial values from allison2014b
    X = [0.01,0.12,0.05,6,0.0067,0.04,0.15];
    %----------------------------------------------------------------------
    case 2 % initial values from allison2014b
    X = [0.01,0.129,0.106,3.57,0.006,0.035,0.17];
    %----------------------------------------------------------------------
    case 3 % (lb+ub)/2
    X = [0.0125,0.2250,0.2600,9.5000,0.0075,0.0550,0.2000];
    %----------------------------------------------------------------------
    case 4 % optimal values from allison2014b
    X = [0.0097,0.0620,0.0201,15.3,0.0061,0.0303,0.170];
    %----------------------------------------------------------------------
    case 5 % optimal values using simultaneous
    X = [0.0200 0.1820 0.0335 10.7308 0.0097 0.0405 0.1700];
    %----------------------------------------------------------------------
    case 6 % optimal values using nested
    X = [0.0166 0.1686 0.0396 6.5015 0.0094 0.0395 0.1700];
    %----------------------------------------------------------------------
    case 7 % feasible inner-loop
    X = [0.0199,0.1567,0.03835,7.95,0.006528,0.05099,0.16999];
    %----------------------------------------------------------------------
end

%% --- global search method
% if genetic algorithm global search
if p.GAflag

    % fix random seed
    rng(56438790)

    % ga options
    options = optimoptions('ga','Display','Iter','UseParallel',true,...
        'PopulationSize',nxp*5,'MaxGenerations',1,...
        'CreationFcn',@gacreationnonlinearfeasible);

    % run ga
    [X,~] = ga(@(x) InnerLoopCombined(x,p),length(X),A,b,[],[],LB,UB,...
        @(x) DSuspensionNested_DesignConstraints(x,p),[],options);

end

%% --- local search method
% fmincon options
options = optimoptions('fmincon','Display','Iter','Algorithm','interior-point',...
    'UseParallel',true,'MaxIter',1000,'MaxFunctionEvaluations',3000,...
    'FiniteDifferenceStepSize',10*sqrt(eps),'OptimalityTolerance',p.OuterLoopTolerance,...
    'FiniteDifferenceType',p.FiniteDifferenceType,...
    'ScaleProblem','obj-and-constr','HonorBounds',true);

% run fmincon
[X,F] = fmincon(@(x) InnerLoopCombined(x,p),X,A,b,[],[],LB,UB,@(x) DSuspensionNested_DesignConstraints(x,p),options);

%% return final plant variables
xp = X;

end

% compute inner-loop problems
function [F] = InnerLoopCombined(x,p)

% use simple inner-loop feasibility checks
feasible = DSuspensionNested_InnerLoopFeasible(x,p);

% stop if inner-loop is infeasible
if ~feasible
    F = nan;
    return
end

% ramp load case
[Framp] = DSuspensionNested_InnerLoopRamp(x,p);

% rough road load case
[Frough] = DSuspensionNested_InnerLoopRough(x,p);

% combine
F = p.wramp*Framp + p.wrough*Frough;

end