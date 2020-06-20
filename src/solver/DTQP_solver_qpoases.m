%--------------------------------------------------------------------------
% DTQP_SOLVER_qpoases.m
% Interface to qpOASES software
%--------------------------------------------------------------------------
% See https://projects.coin-or.org/qpOASES
% NOTE: initial implementation
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [X,F,in,opts] = DTQP_SOLVER_qpoases(H,f,A,b,Aeq,beq,lb,ub,in,opts)

% qp options
qp = opts.qp;

% options
options = qpOASES_options('default');
options.printLevel = qp.printLevel;
options.maxIter = qp.maxIter;
options.maxCpuTime = qp.maxCpuTime;

% number of optimization variables
n = in.nx;

% require f to be the correct size
if isempty(f)
    f = zeros(n,1);
end

% arbitrary large number
Nlarge = 1e9;

% remove simple upper bound constraints with no bound (inf)
ub(isinf(ub)) = Nlarge;
lb(isinf(lb)) = -Nlarge;

% combine matrices/vectors
AL = [A;Aeq];
lbL = full([-Nlarge*ones(size(b));beq]); % require full matrix
ubL = full([b;beq]); % require full matrix

% solve the QP
[X, F, EXITFLAG] = qpOASES(H,f,AL,lb,ub,lbL,ubL,options);

end