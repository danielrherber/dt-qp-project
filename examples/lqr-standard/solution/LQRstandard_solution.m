%--------------------------------------------------------------------------
% LQRstandard_solution.m
% Solution function for LQRstandard example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function D = LQRstandard_solution(A,B,R,Q,M,p,opts)

% sqrt of the number of costates
p.snp = p.ns;

% copy the matrices
p.A = full(A); p.B = full(B); p.iR = inv(R); p.Q = full(Q); p.M = full(M);

% find indices of diagonal and lower triangular entries
p.NS = sum(sum(tril(ones(p.snp,p.snp))));
p.Ilower = find(tril(ones(p.snp,p.snp)));
p.Idiag = find(eye(p.snp,p.snp));

%--------------------------------------------------------------------------
% START: ode solution
%--------------------------------------------------------------------------
% copy p to output structure
pode = p;

% ode options
options = odeset('RelTol',opts.tolode,'AbsTol',opts.tolode*1e-3,'InitialStep',1e-6);

% initial costates
P0 = reshape(full(M),[],1);
P0 = P0(p.Ilower);

% backward integration for the costates
[tode,Pode] = ode15s(@(x,y) odefun_dP(x,y,p),[p.t(end) p.t(1)],P0,options);

% flip the solution to be forward in time
tode = flipud(tode);
Pode = flipud(Pode);

% same time vector for ode/bvp
pode.t = tode;

% forward integration for the states
[~,Yode] = ode15s(@(x,y) odefun_dX(x,y,pode,Pode),pode.t,p.x0,options);

% calculate the optimal control
Uode = calcU(Yode,Pode,p,tode);

%--------------------------------------------------------------------------
% END: ode solution
%--------------------------------------------------------------------------

% obtain bvp solution or report the ode solution
if strcmp(opts.solmethod,'bvp')
    %----------------------------------------------------------------------
    % START: bvp solution
    %----------------------------------------------------------------------
    % ordered nodes of the initial mesh
    T = p.t;

    % initial guess for the solution as intepolation of previous solutions
    yinit = @(t) [interp1(tode,Yode,t,'pchip'),interp1(tode,Pode,t,'pchip')];
    
    % initialize solution
    solinit = bvpinit(T,yinit);
    
    % bvp options
    options = bvpset('RelTol',opts.tolbvp,'AbsTol',opts.tolbvp*1e-3,...
        'NMax',10000,'Stats','on');
    
    % solve the bvp
    sol = bvp4c(@(x,y) odefun(x,y,p),@(ya,yb) bcfun(ya,yb,p),solinit,options);
    
    % time mesh
    D.T = sol.x';
    
    % states
    D.Y = sol.y(1:p.ns,:)';
    
    % costates
    Psol = sol.y(p.ns+1:end,:)';
    
    % calculate the optimal control
    D.U = calcU(D.Y,Psol,p,D.T);
    %----------------------------------------------------------------------
    % END: bvp solution
    %----------------------------------------------------------------------
else
    % use the ode-based solution
    D.T = tode;
    D.Y = Yode;
    D.U = Uode;
    
end

%--------------------------------------------------------------------------
% START: objective function
%--------------------------------------------------------------------------
% Mayer term
Yf = D.Y(end,:)'; % final states
FM = Yf'*M*Yf;

% Lagrange term
FQ = integral(@(t) quadIntegrand(t,D.T,D.Y,Q),D.T(1),D.T(end),...
    'RelTol',opts.tolode,'AbsTol',opts.tolode*1e-3); % Y'*Q*Y
FR = integral(@(t) quadIntegrand(t,D.T,D.U,R),D.T(1),D.T(end),...
    'RelTol',opts.tolode,'AbsTol',opts.tolode*1e-3); % U'*R*U
FL = FQ + FR; % combine

% objective function value
D.F = FM + FL; % combine
%--------------------------------------------------------------------------
% END: objective function
%--------------------------------------------------------------------------

end
% costate ordinary differential equation function
function dP = odefun_dP(~,P,p)
% reshape the costates
q = P;
P = zeros(p.snp,p.snp);
P(p.Ilower) = q;
Pdiag = diag(P);
P = P+P';
P(p.Idiag) = Pdiag;
% co-state equation
dP = -P*p.A - p.A'*P - p.Q + P*p.B*p.iR*p.B'*P;
% reshape
dP = dP(p.Ilower);
end
% state ordinary differential equation function for ode option
function dX = odefun_dX(x,X,p,P)
P = interp1(p.t,P,x,'pchip'); % costates
% reshape the costates
q = P;
P = zeros(p.snp,p.snp);
P(p.Ilower) = q;
Pdiag = diag(P);
P = P+P';
P(p.Idiag) = Pdiag;
% state equation
dX = (p.A - p.B*p.iR*p.B'*P)*X;
end
% ordinary differential equation function for bvp option
function dY = odefun(~,Y,p)
% get matrices
A = p.A; B = p.B; Q = p.Q; iR = p.iR;
% extract
X = Y(1:p.ns); % states
P = Y(p.ns+1:end); % costates
% reshape the states
X = reshape(X,[],1);
% reshape the costates
q = P;
P = zeros(p.snp,p.snp);
P(p.Ilower) = q;
Pdiag = diag(P);
P = P+P';
P(p.Idiag) = Pdiag;
% co-state equation
dP = -(P*A + A'*P + Q - P*B*iR*B'*P);
% control equation
U = -iR*B'*P*X;
% state equation
dX = A*X + B*U;
% reshape
dX = reshape(dX,[],1);
dP = dP(p.Ilower);
% combine
dY = [dX;dP];
end
% boundary value function for bvp option
function res = bcfun(Y0,Yf,p)
X0 = Y0(1:p.ns); % initial state 
Pf = Yf(p.ns+1:end); % costate final conditions
% reshape
X0 = reshape(X0,[],1);
% Pf = reshape(Pf,p.snp,p.snp);
% residual equations
res1 = X0 - p.x0;
res2 = Pf - p.M(p.Ilower);
% reshape
res1 = reshape(res1,[],1);
res2 = reshape(res2,[],1);
% combine
res = [res1; res2];
end
% calculate the optimal control from states and costates
function U = calcU(Y,P,p,T)
    % intialize
    U = zeros(size(p.B,2),length(T));
    for k = 1:length(T)
        % states
        YY = reshape(Y(k,:),[],1);
        % costates
        q = P(k,:);
        PP = zeros(p.snp,p.snp);
        PP(p.Ilower) = q;
        Pdiag = diag(PP);
        PP = PP+PP';
        PP(p.Idiag) = Pdiag;
        % control equation
        U(:,k) = -p.iR*p.B'*PP*YY;
    end
    % transpose
    U = U';
end
% calculate quadratic integrand (vectorized)
function I = quadIntegrand(t,T,A,B)
    I = zeros(size(t));
    for k = 1:length(t)
        X = interp1(T,A,t(k),'pchip')';
        I(k) = X'*B*X;
    end
end