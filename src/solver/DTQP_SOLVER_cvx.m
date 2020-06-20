%--------------------------------------------------------------------------
% DTQP_SOLVER_cvx.m
% Interface to cvx software
%--------------------------------------------------------------------------
% See http://cvxr.com/cvx/
% NOTE: initial implementation
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [X,F,in,opts] = DTQP_SOLVER_cvx(H,f,A,b,Aeq,beq,lb,ub,in,opts)

% number of optimization variables
n = in.nx;

% options
cvx_precision('best');

% find simple equality constraints in ub/lb
ieq = find(lb==ub);
Aeq = [Aeq;sparse(1:length(ieq),ieq,1,length(ieq),n)];
beq = [beq;ub(ieq)];
lb(ieq) = -inf;
ub(ieq) = inf;

% remove simple upper bound constraints with no bound (inf)
iub = find(~isinf(ub));
rub = ub(iub);
Aub = sparse(1:length(rub),iub,1,length(rub),n);

% remove simple lower bound constraints with no bound (-inf)
ilb = find(~isinf(lb));
rlb = lb(ilb);
Alb = sparse(1:length(rlb),ilb,1,length(rlb),n);

% combine
A = [A;Aub;Alb];
b = [b;rub;rlb];

% objective function string
strObj = '';
if ~isempty(H)
    strObj = [strObj,'0.5*quad_form(x,H) '];
end
if ~isempty(f)
    strObj = [strObj,'f''*x '];
end
if ~isempty(strObj)
	strObj = ['minimize ( ',strObj,')'];
end

% linear inequality constraint string
strA = '';
if size(A,1) > 0
    strA = 'A*x <= b;';
end

% linear equality constraint string
strAeq = '';
if size(Aeq,1) > 0
    strAeq = 'Aeq*x == beq;';
end

% construct the cvx problem
cvx_begin
    variable x(n)
    eval(strObj)
    eval(strA)
    eval(strAeq)
cvx_end

% assign the outputs
X = x;
F = cvx_optval;

end