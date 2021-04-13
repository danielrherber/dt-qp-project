%--------------------------------------------------------------------------
% DTQP_TEST_FO.m
% Test new approaches and functions for the first-order hold method
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clear; clc

% random system
% rng(235930) % nonsingular A
rng(45354) % singular A
ny = 3; % number of states
nu = 2; % number of inputs
sys = rss(ny,1,nu); % random state-space model

% extract matrices
A = sys.A;
B = sys.B;

% needed data
p = [];
nt = 20;
t = linspace(0,100,nt)';
h = diff(t);
np = 0;
nd = 0;

% assign
in.p = p;
in.nt= nt;
in.t = t;
in.h = h;
in.ny = ny;
in.nu = nu;
in.np = np;
in.nd = nd;
in.nx = (ny+nu)*nt + np;
in.tm = t(1:end-1) + in.h/2; % midpoints

% c1
C1 = DTQP_DEFECTS_expm(A,in,[]);

% c2
C2 = DTQP_DEFECTS_convolution_integral_type0(A,B,in,[]);

% c3
C3 = DTQP_DEFECTS_convolution_integral_type1(A,B,in,[]);

% select input
u_num = 4;

switch u_num
    %----------------------------------------------------------------------
    case 1
    u_fun = @(t) [sin(t);exp(-t)];
    %----------------------------------------------------------------------
    case 2
    u_fun = @(t) [t;-t];
    %----------------------------------------------------------------------
    case 3
	u_fun = @(t) [sin(t/10);0];
    %----------------------------------------------------------------------
    case 4
    t_in = t;
    U = rand(length(t),nu);
    u_fun = @(t) interp1(t_in,U,t,'linear','extrap')';
    %----------------------------------------------------------------------
end

% simulation options
OPTIONS = odeset('RelTol',1e-13,'AbsTol',1e-20);

% time horizon
TSPAN = [t(1) t(end)];

% initial conditions
X0 = rand(ny,1);

% forward simulation
[T_ode2,X_ode2] = ode45(@(t,x) deriv(t,x,u_fun,A,B),TSPAN,X0,OPTIONS);
[T_ode,X_ode] = ode45(@(t,x) deriv(t,x,u_fun,A,B),t,X0,OPTIONS);

% initialize
T_foh = t;
X_foh = nan(nt,ny);
X_foh(1,:) = X0(:)'; % initial value

% shift for easy indexing
C1 = shiftdim(C1,1);
C2 = shiftdim(C2,1);
C3 = shiftdim(C3,1);

% initial input value
U_foh(:,1) = u_fun(t(1));

% go through each interval
for k = 1:nt-1

    % end point input value
    U_foh(:,k+1) = u_fun(t(k+1));

    % extract
    c1 = C1(:,:,k);
    c2 = C2(:,:,k);
    c3 = C3(:,:,k);

    % compute next states
    X_foh(k+1,:) = c1*X_foh(k,:)' + c2*U_foh(:,k) + c3*(U_foh(:,k+1) - U_foh(:,k));

end

% compute states using some other approaches
X_tr = compare_with_TR(A,B,X0,in,[],U_foh);
X_hs = compare_with_HS(A,B,X0,in,[],U_foh);
X_fo = compare_with_FO(A,B,X0,in,[],U_foh);

% defect constraint matrices
[Aeq_HS,beq_HS] = DTQP_DEFECTS_HS(A,B,[],[],in,[]);
[Aeq_FO,beq_FO] = DTQP_DEFECTS_FO(A,B,[],[],in,[]);

%% plot states
hf = figure; hf.Color = 'w'; hold on
plot(T_ode2,X_ode2,'linewidth',2)
plot(T_foh,X_foh,'.','markersize',20)
plot(T_foh,X_tr,'.-','markersize',6)
plot(T_foh,X_hs,'.-','markersize',6)
plot(T_foh,X_fo,'.-','markersize',6)
ha = gca; ha.LineWidth = 1; ha.FontSize = 14;

%% plot controls
hf = figure; hf.Color = 'w'; hold on
U_ode = zeros(nu,length(T_ode2));
for k = 1:length(T_ode2)
   U_ode(:,k) =  u_fun(T_ode2(k));
end
plot(T_ode2,U_ode,'linewidth',2)
plot(T_foh,U_foh,'.','markersize',20)
ha = gca; ha.LineWidth = 1; ha.FontSize = 14;

%% plot errors
hf = figure; hf.Color = 'w'; hold on
plot(t,abs(X_foh-X_ode),'k.-','markersize',10)
plot(t,abs(X_hs-X_ode),'r.-','markersize',10)
plot(t,abs(X_tr-X_ode),'b.-','markersize',10)
plot(t,abs(X_fo-X_ode),'g.-','markersize',10)
ha = gca; ha.LineWidth = 1; ha.FontSize = 14;
ha.YScale = 'log';

% state derivative function
function Dx = deriv(t,x,u_fun,A,B)

% input
u = u_fun(t);

% compute state derivative
Dx = A*x + B*u;

end

% compute states using simultaneous HS defect constraints
function X_tr = compare_with_TR(A,B,X0,in,opts,U_foh)

% defect constraint matrices
[Aeq,beq] = DTQP_DEFECTS_TR(A,B,[],[],in,opts);

% add initial condition
Aeq(end+1,(in.nu+0)*in.nt+1) = 1;
Aeq(end+1,(in.nu+1)*in.nt+1) = 1;
Aeq(end+1,(in.nu+2)*in.nt+1) = 1;
beq(end+1,1) = X0(1);
beq(end+1,1) = X0(2);
beq(end+1,1) = X0(3);

% transpose inputs
Ut = U_foh';

% solve for states
X = Aeq(:,in.nu*in.nt+1:end)\(beq-Aeq(:,1:(in.nu*in.nt))*Ut(:));

% reshape
X_tr = reshape(X,in.nt,[]);

end

% compute states using simultaneous HS defect constraints
function X_hs = compare_with_HS(A,B,X0,in,opts,U_foh)

% defect constraint matrices
[Aeq,beq] = DTQP_DEFECTS_HS(A,B,[],[],in,opts);

% add initial condition
Aeq(end+1,(in.nu+0)*in.nt+1) = 1;
Aeq(end+1,(in.nu+1)*in.nt+1) = 1;
Aeq(end+1,(in.nu+2)*in.nt+1) = 1;
beq(end+1,1) = X0(1);
beq(end+1,1) = X0(2);
beq(end+1,1) = X0(3);

% transpose inputs
Ut = U_foh';

% solve for states
X = Aeq(:,in.nu*in.nt+1:end)\(beq-Aeq(:,1:(in.nu*in.nt))*Ut(:));

% reshape-
X_hs = reshape(X,in.nt,[]);

end

% compute states using simultaneous FO defect constraints
function X_fo = compare_with_FO(A,B,X0,in,opts,U_foh)

% defect constraint matrices
[Aeq,beq] = DTQP_DEFECTS_FO(A,B,[],[],in,opts);

% add initial condition
Aeq(end+1,(in.nu+0)*in.nt+1) = 1;
Aeq(end+1,(in.nu+1)*in.nt+1) = 1;
Aeq(end+1,(in.nu+2)*in.nt+1) = 1;
beq(end+1,1) = X0(1);
beq(end+1,1) = X0(2);
beq(end+1,1) = X0(3);

% transpose inputs
Ut = U_foh';

% solve for states
X = Aeq(:,in.nu*in.nt+1:end)\(beq-Aeq(:,1:(in.nu*in.nt))*Ut(:));

% reshape
X_fo = reshape(X,in.nt,[]);

end