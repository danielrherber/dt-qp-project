%--------------------------------------------------------------------------
% OutputTracking_solution.m
% Solution function for OutputTracking example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function D = OutputTracking_solution(p,opts)

% sqrt of the number of costates
p.snp = p.ns;

% copy the matrices
A = p.A; B = p.B; C = p.C; Q = p.Q; R = p.R; o = p.o;

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
P0 = reshape(zeros(p.snp),[],1);
P0 = P0(p.Ilower);
k0 = zeros(p.ns,1);
V0 = [P0;k0];

% backward integration for the costates
[tode,Pkode] = ode15s(@(x,y) odefun_dPk(x,y,p,A,B,C,Q,R,o),[p.t(end) p.t(1)],V0,options);

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
[~,Yode] = ode15s(@(x,y) odefun_dX(x,y,pode,Pode,kode,A,B,C,Q,R,o),pode.t,p.x0,options);

% calculate the optimal control
Uode = calcU(Yode,Pode,kode,p,tode,B,R);

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

    % initial guess for the solution as interpolation of previous solutions
    yinit = @(t) [interp1(tode,Yode,t,'pchip'),interp1(tode,Pode,t,'pchip'),...
        interp1(tode,kode,t,'pchip')];
    
    % initialize solution
    solinit = bvpinit(T,yinit);
    
    % bvp options
    options = bvpset('RelTol',opts.tolbvp,'AbsTol',opts.tolbvp*1e-3,...
        'NMax',10000,'Stats','on');
    
    % solve the bvp
    sol = bvp4c(@(x,y) odefun(x,y,p,A,B,C,Q,R,o),@(ya,yb) bcfun(ya,yb,p),solinit,options);
    
    % time mesh
    D.T = sol.x';
    
    % states
    D.Y = sol.y(1:p.ns,:)';
    
    % costates
    Psol = sol.y(p.ns+1:p.ns+length(p.Ilower),:)';
    
    % additional costates
    ksol = sol.y(p.ns+length(p.Ilower)+1:end,:)';
    
    % calculate the optimal control
    D.U = calcU(D.Y,Psol,ksol,p,D.T,B,R);
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
% actual output values
for k = 1:size(C,1)
    o_dt(:,k) = sum(bsxfun(@times,C(k,:),D.Y),2);
end
% desired output values
o_actual = DTQP_tmatrix(o,p,D.T);

% output errors
e = o_actual - o_dt;

% Lagrange term
FQ = integral(@(t) quadIntegrand(t,D.T,e,Q,p),D.T(1),D.T(end),...
    'RelTol',opts.tolode,'AbsTol',opts.tolode*1e-3); % e'*Q*e
FR = integral(@(t) quadIntegrand(t,D.T,D.U,R,p),D.T(1),D.T(end),...
    'RelTol',opts.tolode,'AbsTol',opts.tolode*1e-3); % U'*R*U
FL = FQ + FR; % combine

% objective function value
D.F = FL; % combine
%--------------------------------------------------------------------------
% END: objective function
%--------------------------------------------------------------------------

end
% costate ordinary differential equation function
function dPk = odefun_dPk(t,Pk,p,A,B,C,Q,R,o)
% matrix values at current time
A = DTQP_tmatrix(A,p,t); A = squeeze(A);
B = DTQP_tmatrix(B,p,t); B = squeeze(B);
C = DTQP_tmatrix(C,p,t); C = squeeze(C);
Q = DTQP_tmatrix(Q,p,t); Q = squeeze(Q);
R = DTQP_tmatrix(R,p,t); R = squeeze(R);
o = DTQP_tmatrix(o,p,t); o = squeeze(o)';
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
dP = -P*A - A'*P - C'*Q*C + P*B*(R\(B'))*P;
% additional costate equation
dk = (P*B*(R\(B')) - A')*k + C'*Q*o;
% reshape
dP = dP(p.Ilower);
% combine
dPk = [dP;dk];
end
% state ordinary differential equation function for ode option
function dX = odefun_dX(t,X,p,P,k,A,B,C,Q,R,o)
% interpolate
P = interp1(p.t,P,t,'pchip'); % costates
k = interp1(p.t,k,t,'pchip'); % additional costates
k = k(:);
% matrix values at current time
A = DTQP_tmatrix(A,p,t); A = squeeze(A);
B = DTQP_tmatrix(B,p,t); B = squeeze(B);
R = DTQP_tmatrix(R,p,t); R = squeeze(R);
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
dX = A*X + B*U;
disp(t)
end
% ordinary differential equation function for bvp option
function dY = odefun(t,Y,p,A,B,C,Q,R,o)
% matrix values at current time
A = DTQP_tmatrix(A,p,t); A = squeeze(A);
B = DTQP_tmatrix(B,p,t); B = squeeze(B);
C = DTQP_tmatrix(C,p,t); C = squeeze(C);
Q = DTQP_tmatrix(Q,p,t); Q = squeeze(Q);
R = DTQP_tmatrix(R,p,t); R = squeeze(R);
o = DTQP_tmatrix(o,p,t); o = squeeze(o)';
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
dP = -P*A - A'*P - C'*Q*C + P*B*(R\(B'))*P;
% additional costate equation
dk = (P*B*(R\(B')) - A')*k + C'*Q*o;
% control equation
U = -(R\(B'))*(P*X+k);
% state equation
dX = A*X + B*U;
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
res2 = Pf - [zeros(numel(p.Ilower),1);zeros(p.ns,1)];
% reshape
res1 = reshape(res1,[],1);
res2 = reshape(res2,[],1);
% combine
res = [res1; res2];
end
% calculate the optimal control from states and costates
function U = calcU(Y,P,K,p,T,B,R)
    % intialize
    U = zeros(size(B,2),length(T));
    for k = 1:length(T)
        % matrix values at current time
        B = DTQP_tmatrix(B,p,T(k)); B = squeeze(B);
        R = DTQP_tmatrix(R,p,T(k)); R = squeeze(R);
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
function I = quadIntegrand(t,T,A,B,p)
    I = zeros(size(t));
    for k = 1:length(t)
        H = DTQP_tmatrix(B,p,t(k)); H = squeeze(H);
        X = interp1(T,A,t(k),'pchip')';
        I(k) = X'*H*X;
    end
end