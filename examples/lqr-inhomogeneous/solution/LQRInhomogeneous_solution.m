%--------------------------------------------------------------------------
% LQRInhomogeneous_solution.m
% Solution function for LQRInhomogeneous example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function D = LQRInhomogeneous_solution(A,B,d,R,Q,M,p,opts)

% sqrt of the number of costates
p.snp = p.ns;

% copy the matrices
p.A = A; p.B = B; p.d = d; p.R = R; p.Q = Q; p.M = M;

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
k0 = zeros(p.ns,1);
V0 = [P0;k0];

% backward integration for the costates
[tode,Pkode] = ode15s(@(x,y) odefun_dPk(x,y,p),[p.t(end) p.t(1)],V0,options);

% extract
Pode = Pkode(:,1:length(p.Ilower));
kode = Pkode(:,length(p.Ilower)+1:end);

% flip the solution to be forward in time
tode = flipud(tode);
Pode = flipud(Pode);
kode = flipud(kode);

% same time vector for ode/bvp
pode.t = tode;

% forward integration for the states
[~,Yode] = ode15s(@(x,y) odefun_dX(x,y,pode,Pode,kode),pode.t,p.x0,options);

% calculate the optimal control
Uode = calcU(Yode,Pode,kode,p,tode);

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
    yinit = @(t) [interp1(tode,Yode,t,'pchip'),interp1(tode,Pode,t,'pchip'),...
        interp1(tode,kode,t,'pchip')];
    
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
    Psol = sol.y(p.ns+1:p.ns+length(p.Ilower),:)';
    
    % additional costates
    ksol = sol.y(p.ns+length(p.Ilower)+1:end,:)';
    
    % calculate the optimal control
    D.U = calcU(D.Y,Psol,ksol,p,D.T);
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
D.F = FM/2 + FL/2; % combine
%--------------------------------------------------------------------------
% END: objective function
%--------------------------------------------------------------------------

end
% costate ordinary differential equation function
function dPk = odefun_dPk(t,Pk,p)
% matrix values at current time
p.t = t;
A = DTQP_tmatrix(p.A,p);
B = DTQP_tmatrix(p.B,p);
d = DTQP_tmatrix(p.d,p);
Q = DTQP_tmatrix(p.Q,p);
R = DTQP_tmatrix(p.R,p);
% extract
P = Pk(1:length(p.Ilower));
k = Pk(length(p.Ilower)+1:end);
% reshape the costates
q = P;
P = zeros(p.snp,p.snp);
P(p.Ilower) = q;
Pdiag = diag(P);
P = P+P';
P(p.Idiag) = Pdiag;
% costate equation
dP = -P*A - A'*P - Q + P*B*(R\(B'))*P;
% additional costate equation
dk = (P*B*(R\(B')) - A')*k - P*d;
% reshape
dP = dP(p.Ilower);
% combine
dPk = [dP;dk];
end
% state ordinary differential equation function for ode option
function dX = odefun_dX(t,X,p,P,k)
% interpolate
P = interp1(p.t,P,t,'pchip'); % costates
k = interp1(p.t,k,t,'pchip'); % additional costates
k = k(:);
% matrix values at current time
p.t = t;
A = DTQP_tmatrix(p.A,p);
B = DTQP_tmatrix(p.B,p);
d = DTQP_tmatrix(p.d,p);
R = DTQP_tmatrix(p.R,p);
% reshape the costates
q = P;
P = zeros(p.snp,p.snp);
P(p.Ilower) = q;
Pdiag = diag(P);
P = P+P';
P(p.Idiag) = Pdiag;
% control
U = -(R\(B'))*(P*X+k);
% state equation
dX = A*X + B*U + d;
disp(t)
end
% ordinary differential equation function for bvp option
function dY = odefun(t,Y,p)
% matrix values at current time
p.t = t;
A = DTQP_tmatrix(p.A,p);
B = DTQP_tmatrix(p.B,p);
d = DTQP_tmatrix(p.d,p);
Q = DTQP_tmatrix(p.Q,p);
R = DTQP_tmatrix(p.R,p);
% extract
X = Y(1:p.ns); % states
P = Y(p.ns+1:p.ns+length(p.Ilower)); % costates
k = Y(p.ns+length(p.Ilower)+1:end); % additional costates
% reshape the states and additional costates
X = reshape(X,[],1);
k = reshape(k,[],1);
% reshape the costates
q = P;
P = zeros(p.snp,p.snp);
P(p.Ilower) = q;
Pdiag = diag(P);
P = P+P';
P(p.Idiag) = Pdiag;
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
dP = dP(p.Ilower);
% combine
dY = [dX;dP;dk];
end
% boundary value function for bvp option
function res = bcfun(Y0,Yf,p)
X0 = Y0(1:p.ns); % initial state 
Pf = Yf(p.ns+1:end); % costate and additional costate final conditions
% residual equations
res1 = X0 - p.x0;
res2 = Pf - [p.M(p.Ilower);zeros(p.ns,1)];
% reshape
res1 = reshape(res1,[],1);
res2 = reshape(res2,[],1);
% combine
res = [res1; res2];
end
% calculate the optimal control from states and costates
function U = calcU(Y,P,K,p,T)
    % intialize
    U = zeros(size(p.B,2),length(T));
    for k = 1:length(T)
        % matrix values at current time
        p.t = T(k);
        B = DTQP_tmatrix(p.B,p);
        R = DTQP_tmatrix(p.R,p);
        % states
        YY = reshape(Y(k,:),[],1);
        % costates
        q = P(k,:);
        PP = zeros(p.snp,p.snp);
        PP(p.Ilower) = q;
        Pdiag = diag(PP);
        PP = PP+PP';
        PP(p.Idiag) = Pdiag;
        % additional costates
        KK = reshape(K(k,:),[],1);
        % control equation
        U(:,k) = -(R\(B'))*(PP*YY+KK);
    end
    % transpose
    U = U';
end
% calculate quadratic integrand (vectorized)
function I = quadIntegrand(t,T,A,B)
    I = zeros(size(t));
    for k = 1:length(t)
        p.t = t(k);
        H = DTQP_tmatrix(B,p);
        X = interp1(T,A,t(k),'pchip')';
        I(k) = X'*H*X;
    end
end