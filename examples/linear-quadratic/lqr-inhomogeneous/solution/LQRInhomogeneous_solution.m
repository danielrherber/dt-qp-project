%--------------------------------------------------------------------------
% LQRInhomogeneous_solution.m
% Solution function for LQRInhomogeneous example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function D = LQRInhomogeneous_solution(in,opts)

% sqrt of the number of costates (or the number of states)
in.ny = in.phase_info(1).ny;
in.snp = in.ny;

% copy the matrices
auxdata = in.auxdata; A = auxdata.A; B = auxdata.B; d = auxdata.d;
R = auxdata.R; Q = auxdata.Q; M = auxdata.M; x0 = auxdata.x0;

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
options = odeset('RelTol',opts.tolode,'AbsTol',opts.tolode*1e-3,'InitialStep',1e-12);

% initial costates
P0 = reshape(full(M),[],1);
P0 = P0(in.Ilower);
k0 = zeros(in.ny,1);
V0 = [P0;k0];

% backward integration for the costates
[tode,Pkode] = ode15s(@(x,y) odefun_dPk(x,y,A,B,d,Q,R,in),[in.tf in.t0],V0,options);

% extract
Pode = Pkode(:,1:length(in.Ilower));
kode = Pkode(:,length(in.Ilower)+1:end);

% flip the solution to be forward in time
tode = flipud(tode);
Pode = flipud(Pode);
kode = flipud(kode);

% same time vector for ode/bvp
pode.t = tode;

% create interpolation functions
Pode_interp = griddedInterpolant(tode,Pode,'spline');
kode_interp = griddedInterpolant(tode,kode,'spline');

% forward integration for the states
[~,Yode] = ode15s(@(x,y) odefun_dX(x,y,A,B,d,R,pode,Pode_interp,kode_interp),pode.t,x0,options);

% calculate the optimal control
Uode = calcU(Yode,Pode,kode,B,R,in,tode);

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
    yinit = @(t) [interp1(tode,Yode,t,'spline'),interp1(tode,Pode,t,'spline'),...
        interp1(tode,kode,t,'spline')];

    % initialize solution
    solinit = bvpinit(T,yinit);

    % bvp options
    options = bvpset('RelTol',opts.tolbvp,'AbsTol',opts.tolbvp*1e-3,...
        'NMax',10000,'Stats','on');

    % solve the bvp
    sol = bvp4c(@(x,y) odefun(x,y,A,B,d,Q,R,in),@(ya,yb) bcfun(ya,yb,x0,M,in),solinit,options);

    % time mesh
    D.T = sol.x';

    % states
    D.Y = sol.y(1:in.ny,:)';

    % costates
    Psol = sol.y(in.ny+1:in.ny+length(in.Ilower),:)';

    % additional costates
    ksol = sol.y(in.ny+length(in.Ilower)+1:end,:)';

    % calculate the optimal control
    D.U = calcU(D.Y,Psol,ksol,B,R,in,D.T);
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
FQ = integral(@(t) quadIntegrand(t,D.T,D.Y,Q,in),D.T(1),D.T(end),...
    'RelTol',opts.tolode,'AbsTol',opts.tolode*1e-3); % Y'*Q*Y
FR = integral(@(t) quadIntegrand(t,D.T,D.U,R,in),D.T(1),D.T(end),...
    'RelTol',opts.tolode,'AbsTol',opts.tolode*1e-3); % U'*R*U
FL = FQ + FR; % combine

% objective function value
D.F = FM/2 + FL/2; % combine
%--------------------------------------------------------------------------
% END: objective function
%--------------------------------------------------------------------------

end
% costate ordinary differential equation function
function dPk = odefun_dPk(t,Pk,A,B,d,Q,R,in)
% matrix values at current time
auxdata = in.auxdata;
A = DTQP_tmatrix(A,auxdata,t); A = squeeze(A);
B = DTQP_tmatrix(B,auxdata,t); B = squeeze(B);
d = DTQP_tmatrix(d,auxdata,t); d = squeeze(d); d = d(:);
Q = DTQP_tmatrix(Q,auxdata,t); Q = squeeze(Q);
R = DTQP_tmatrix(R,auxdata,t); R = squeeze(R);
% extract
P = Pk(1:length(in.Ilower));
k = Pk(length(in.Ilower)+1:end);
% reshape the costates
q = P;
P = zeros(in.snp,in.snp);
P(in.Ilower) = q;
Pdiag = diag(P);
P = P+P';
P(in.Idiag) = Pdiag;
% costate equation
dP = -P*A - A'*P - Q + P*B*(R\(B'))*P;
% additional costate equation
dk = (P*B*(R\(B')) - A')*k - P*d;
% reshape
dP = dP(in.Ilower);
% combine
dPk = [dP;dk];
end
% state ordinary differential equation function for ode option
function dX = odefun_dX(t,X,A,B,d,R,in,Pode_interp,kode_interp)

% interpolate
P = Pode_interp(t);
k = kode_interp(t); % additional costates
k = k(:);
% matrix values at current time
auxdata = in.auxdata;
A = DTQP_tmatrix(A,auxdata,t); A = squeeze(A);
B = DTQP_tmatrix(B,auxdata,t); B = squeeze(B);
d = DTQP_tmatrix(d,auxdata,t); d = squeeze(d); d = d(:);
R = DTQP_tmatrix(R,auxdata,t); R = squeeze(R);
% reshape the costates
q = P;
P = zeros(in.snp,in.snp);
P(in.Ilower) = q;
Pdiag = diag(P);
P = P+P';
P(in.Idiag) = Pdiag;
% control
U = -(R\(B'))*(P*X+k);
% state equation
dX = A*X + B*U + d;
disp(t)
end
% ordinary differential equation function for bvp option
function dY = odefun(t,Y,A,B,d,Q,R,in)
% matrix values at current time
auxdata = in.auxdata;
A = DTQP_tmatrix(A,auxdata,t); A = squeeze(A);
B = DTQP_tmatrix(B,auxdata,t); B = squeeze(B);
d = DTQP_tmatrix(d,auxdata,t); d = squeeze(d); d = d(:);
Q = DTQP_tmatrix(Q,auxdata,t); Q = squeeze(Q);
R = DTQP_tmatrix(R,auxdata,t); R = squeeze(R);
% extract
X = Y(1:in.ny); % states
P = Y(in.ny+1:in.ny+length(in.Ilower)); % costates
k = Y(in.ny+length(in.Ilower)+1:end); % additional costates
% reshape the states and additional costates
X = reshape(X,[],1);
k = reshape(k,[],1);
% reshape the costates
q = P;
P = zeros(in.snp,in.snp);
P(in.Ilower) = q;
Pdiag = diag(P);
P = P+P';
P(in.Idiag) = Pdiag;
% co-state equation
dP = -P*A - A'*P - Q + P*B*(R\(B'))*P;
% additional costate equation
dk = (P*B*(R\(B')) - A')*k - P*d;
% control equation
U = -(R\(B'))*(P*X+k);
% state equation
dX = A*X + B*U + d;
% reshape
dX = reshape(dX,[],1);
dP = dP(in.Ilower);
% combine
dY = [dX;dP;dk];
end
% boundary value function for bvp option
function res = bcfun(Y0,Yf,x0,M,in)
X0 = Y0(1:in.ny); % initial state
Pf = Yf(in.ny+1:end); % costate and additional costate final conditions
% residual equations
res1 = X0 - x0;
res2 = Pf - [M(in.Ilower);zeros(in.ny,1)];
% reshape
res1 = reshape(res1,[],1);
res2 = reshape(res2,[],1);
% combine
res = [res1; res2];
end
% calculate the optimal control from states and costates
function U = calcU(Y,P,K,B,R,in,T)
    % intialize
    U = zeros(size(B,2),length(T));
    for k = 1:length(T)
        % matrix values at current time
        auxdata = in.auxdata;
        B = DTQP_tmatrix(B,auxdata,T(k)); B = squeeze(B);
        R = DTQP_tmatrix(R,auxdata,T(k)); R = squeeze(R);
        % states
        YY = reshape(Y(k,:),[],1);
        % costates
        q = P(k,:);
        PP = zeros(in.snp,in.snp);
        PP(in.Ilower) = q;
        Pdiag = diag(PP);
        PP = PP+PP';
        PP(in.Idiag) = Pdiag;
        % additional costates
        KK = reshape(K(k,:),[],1);
        % control equation
        U(:,k) = -(R\(B'))*(PP*YY+KK);
    end
    % transpose
    U = U';
end
% calculate quadratic integrand (vectorized)
function I = quadIntegrand(t,T,A,B,in)
    A_interp = griddedInterpolant(T,A,'spline')';
    I = zeros(size(t));
    for k = 1:length(t)
        auxdata = in.auxdata;
        H = DTQP_tmatrix(B,auxdata,t(k)); H = squeeze(H);
        X = A_interp(t(k))';
        I(k) = X'*H*X;
    end
end