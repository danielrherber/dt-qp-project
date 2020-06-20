%--------------------------------------------------------------------------
% DTQP_QLIN_update_control_constraint.m
% Convert the nonlinear control constraint to a Linear constraint for the
% LQDO problem
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [setup] = DTQP_QLIN_update_control_constraint(setup,opts)

% get the nonlinear variables
symb = setup.symb;
Linc = symb.c;
oc = symb.oc;

% get the linearized matrix values
U = rand(opts.dt.nt,oc.nu);
T = linspace(setup.p.t0,setup.p.tf,opts.dt.nt)';
X = [U];
param = [];

B = TaylorConvert(Linc.B,T,X,param);
d = TaylorConvert(Linc.d,T,X,param);

% linear equality constraints
Y(1).linear(1).right = 1;Y(1).linear(1).matrix = B;
Y(1).b = d;

% return the values
setup.Y = Y;

end