%--------------------------------------------------------------------------
% TavallaeiTousi1_solution.m
% Solution function for TavallaeiTousi1 example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function D = TavallaeiTousi1_solution(in,opts)

% sqrt of the number of costates (or the number of states)
in.ny = in.phase_info(1).ny;
in.snp = in.ny;

% copy the matrices
auxdata = in.auxdata; A = full(auxdata.A); B = full(auxdata.B);
iR = full(inv(auxdata.r)); R = full(auxdata.r);
Q = full(auxdata.Q); x0 = full(auxdata.y0); K = auxdata.k;
M = zeros(in.ny);

% find indices of diagonal and lower triangular entries
in.NS = sum(sum(tril(ones(in.snp,in.snp))));
in.Ilower = find(tril(ones(in.snp,in.snp)));
in.Idiag = find(eye(in.snp,in.snp));

%--------------------------------------------------------------------------
% START: ode solution
%--------------------------------------------------------------------------
% copy in to output structure
pode = in;

% ode options
options = odeset('RelTol',opts.tolode,'AbsTol',opts.tolode*1e-3,'InitialStep',1e-6);

% initial costates
P0 = reshape(full(M),[],1);
P0 = P0(in.Ilower);

% backward integration for the costates
[tode,Pode] = ode15s(@(x,y) odefun_dP(x,y,A,B,iR,Q,K,in),[in.tf in.t0],P0,options);

% flip the solution to be forward in time
tode = flipud(tode);
Pode = flipud(Pode);

% same time vector for ode/bvp
pode.t = tode;

% create interpolation functions
Pode_interp = griddedInterpolant(tode,Pode,'spline');

% forward integration for the states
[~,Yode] = ode15s(@(x,y) odefun_dX(x,y,A,B,iR,pode,Pode_interp),pode.t,x0,options);

% calculate the optimal control
Uode = calcU(Yode,Pode,B,iR,in,tode);

%--------------------------------------------------------------------------
% END: ode solution
%--------------------------------------------------------------------------

% obtain bvp solution or report the ode solution
if strcmp(opts.solmethod,'bvp')
    %----------------------------------------------------------------------
    % START: bvp solution
    %----------------------------------------------------------------------
    % ordered nodes of the initial mesh
    T = in.t;

    % initial guess for the solution as interpolation of previous solutions
    yinit = @(t) [interp1(tode,Yode,t,'pchip'),interp1(tode,Pode,t,'pchip')];

    % initialize solution
    solinit = bvpinit(T,yinit);

    % bvp options
    options = bvpset('RelTol',opts.tolbvp,'AbsTol',opts.tolbvp*1e-3,...
        'NMax',10000,'Stats','on');

    % solve the bvp
    sol = bvp4c(@(x,y) odefun(x,y,A,B,Q,K,iR,in),@(ya,yb) bcfun(ya,yb,x0,M,in),...
            solinit,options);

    % time mesh
    D.T = sol.x';

    % states
    D.Y = sol.y(1:in.ny,:)';

    % costates
    Psol = sol.y(in.ny+1:end,:)';

    % calculate the optimal control
    D.U = calcU(D.Y,Psol,B,iR,in,D.T);
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
FQ = integral(@(t) quadIntegrandTavealleaiTousi1(t,D.T,D.Y,Q,K),D.T(1),D.T(end),...
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
function dP = odefun_dP(t,P,A,B,iR,Q,K,in)
% reshape the costates
q = P;
P = zeros(in.snp,in.snp);
P(in.Ilower) = q;
Pdiag = diag(P);
P = P+P';
P(in.Idiag) = Pdiag;
% co-state equation
dP = -P*A - A'*P - t.^K.*Q + P*B*iR*B'*P;
% reshape
dP = dP(in.Ilower);
end
% state ordinary differential equation function for ode option
function dX = odefun_dX(t,X,A,B,iR,in,Pode_interp)
P = Pode_interp(t); % costates
% reshape the costates
q = P;
P = zeros(in.snp,in.snp);
P(in.Ilower) = q;
Pdiag = diag(P);
P = P+P';
P(in.Idiag) = Pdiag;
% state equation
dX = (A - B*iR*B'*P)*X;
disp(t)
end
% ordinary differential equation function for bvp option
function dY = odefun(t,Y,A,B,Q,K,iR,in)
% extract
X = Y(1:in.ny); % states
P = Y(in.ny+1:end); % costates
% reshape the states
X = reshape(X,[],1);
% reshape the costates
q = P;
P = zeros(in.snp,in.snp);
P(in.Ilower) = q;
Pdiag = diag(P);
P = P+P';
P(in.Idiag) = Pdiag;
% co-state equation
dP = -(P*A + A'*P + t.^K.*Q - P*B*iR*B'*P);
% control equation
U = -iR*B'*P*X;
% state equation
dX = A*X + B*U;
% reshape
dX = reshape(dX,[],1);
dP = dP(in.Ilower);
% combine
dY = [dX;dP];
end
% boundary value function for bvp option
function res = bcfun(Y0,Yf,x0,M,in)
X0 = Y0(1:in.ny); % initial state
Pf = Yf(in.ny+1:end); % costate final conditions
% reshape
X0 = reshape(X0,[],1);
% Pf = reshape(Pf,in.snp,in.snp);
% residual equations
res1 = X0 - x0;
res2 = Pf - M(in.Ilower);
% reshape
res1 = reshape(res1,[],1);
res2 = reshape(res2,[],1);
% combine
res = [res1; res2];
end
% calculate the optimal control from states and costates
function U = calcU(Y,P,B,iR,in,T)
    % initialize
    U = zeros(size(B,2),length(T));
    for k = 1:length(T)
        % states
        YY = reshape(Y(k,:),[],1);
        % costates
        q = P(k,:);
        PP = zeros(in.snp,in.snp);
        PP(in.Ilower) = q;
        Pdiag = diag(PP);
        PP = PP+PP';
        PP(in.Idiag) = Pdiag;
        % control equation
        U(:,k) = -iR*B'*PP*YY;
    end
    % transpose
    U = U';
end
% calculate quadratic integrand (vectorized)
function I = quadIntegrand(t,T,A,B)
    A_interp = griddedInterpolant(T,A,'spline')';
    I = zeros(size(t));
    for k = 1:length(t)
        X = A_interp(t(k))';
        I(k) = X'*B*X;
    end
end
% special quadratic integrand (vectorized) for this example
function I = quadIntegrandTavealleaiTousi1(t,T,A,B,K)
    A_interp = griddedInterpolant(T,A,'spline')';
    I = zeros(size(t));
    for k = 1:length(t)
        X = A_interp(t(k))';
        I(k) = t(k)^K*X'*B*X;
    end
end