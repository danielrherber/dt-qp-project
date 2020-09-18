%--------------------------------------------------------------------------
% LinearPendulum_F.m
% Calculate the closed-form solution for the objective function
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function F = LinearPendulum_F(m,k,umax,x0,v0,tf)

% scaling constants
at = sqrt(m/k); % time
aq = umax/k; % position

% reshape and scale
x0 = x0/aq;
v0 = v0/aq*at;

% number of transitions
n = ceil(tf/at/pi);

% phase shift
p = mod(tf/at,pi);

% go through each phase
for idx = 1:n
    % initial time
    i = max(p + (idx-2)*pi,0);

    % final time
    f = p + (idx-1)*pi;

    % determine current sign of the control
    s = -sign(sin(tf/at))*(-1)^idx;

    % store current initial states
    X0 = x0; V0 = v0;

    % update end phase state values
    x0 = LinearPendulum_Y1(s,i,X0,V0,f);
    v0 = LinearPendulum_Y2(s,i,X0,V0,f);
end

% unscale
x0 = x0*aq;

% objective function value
F = -x0;

end