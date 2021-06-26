%--------------------------------------------------------------------------
% MinimumEnergyTransfer_solution.m
% Solution function for MinimumEnergyTransfer example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [U,Y,F] = MinimumEnergyTransfer_solution(T,in)

% extract
auxdata = in.auxdata;
A = full(auxdata.A);
B = full(auxdata.B);
y0 = auxdata.y0;
yf = auxdata.yf;
t0 = in.t0;
tf = in.tf;

% precompute
BB = B*B';
c1 = (expm(A*(tf-t0))*y0 - yf);

% integral options
opts = {'ArrayValued',true,'RelTol',1e-14};

% controllability gramian integrand
f = @(tau) expm(A*tau)*BB*expm(A'*tau);

% controllability gramian integral
W = @(t) integral(f, t0, t,opts{:});

% compute final value
Wf = W(tf);

% precompute
c2 = Wf\c1;

% compute control
U = zeros(size(B,2),length(T));
for k = 1:length(T)
    U(:,k) = -B'*expm(A'*(tf-T(k)))*c2;
end
U = U';

% simulation options
OPTIONS = odeset('RelTol',1e-13,'AbsTol',1e-13);

% compute states
[~,Y] = ode15s(@(t,y) ODEFUN(t,y,A,B,c2,tf),T,y0,OPTIONS);

% objective function
c3 = (yf - expm(A*(tf-t0))*y0);
F = c3'*(Wf\c3);

end
% state dynamics
function yp = ODEFUN(t,y,A,B,c2,tf)
% derivative function
yp = A*y - B*B'*expm(A'*(tf-t))*c2;

end