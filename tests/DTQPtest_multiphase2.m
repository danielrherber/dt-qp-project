%--------------------------------------------------------------------------
% DTQPtest_multiphase2.m
% Test multiple-interval BrysonDenham problem
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clear; clc

%% tunable parameters
p.ell = 1/9;

t0 = 0; tf = 1;
xT = linspace(t0,tf,20);

opts.dt.nt = 4;
opts.dt.defects = 'PS';
opts.dt.quadrature = 'G';
opts.dt.mesh = 'LGL';

%% setup
% system dynamics
A = [0 1;0 0]; B = [0;1];

% Lagrange term
L(1).left = 1; L(1).right = 1; L(1).matrix = 1/2; % 1/2*u^2

% simple bounds
UB(1).right = 2; UB(1).matrix = [p.ell;Inf]; % states

%% linkage equality constraints
idx = 0;
q = eye(length(A));
for k = 1:length(A)
    idx = idx + 1;
    LY(idx).left.linear.right = 5; % final states
    LY(idx).left.linear.matrix = q(:,k);
    LY(idx).right.linear.right = 4; % initial states
    LY(idx).right.linear.matrix = -q(:,k);
    LY(idx).b = 0;
end

%% phases
% combine structures
for phs = 1:length(xT)-1
    setup(phs).A = A;
    setup(phs).B = B;
    setup(phs).L = L;

    UBt = UB;
    LBt = [];

    if phs == 1
        UBt(2).right = 4; UBt(2).matrix = [0;1]; % initial states
        LBt(1).right = 4; LBt(1).matrix = [0;1];
    end

    if phs == length(xT)-1
        UBt(2).right = 5; UBt(2).matrix = [0;-1]; % final states
        LBt(1).right = 5; LBt(1).matrix = [0;-1];
    end

    setup(phs).UB = UBt;
    if ~isempty(LBt)
        setup(phs).LB = LBt; 
    end

    setup(phs).t0 = xT(phs);
    setup(phs).tf = xT(phs+1);
    setup(phs).p = p;

    if phs < length(xT)-1
        setup(phs).LY = LY;
    end

end

%% solve
[T,U,Y,P,F,p,opts] = DTQP_solve(setup,opts);

%% plot
figure
plot(T,U); hold on

figure
plot(T,Y); hold on