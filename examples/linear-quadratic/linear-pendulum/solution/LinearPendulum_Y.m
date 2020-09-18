%--------------------------------------------------------------------------
% LinearPendulum_Y.m
% Calculate the closed-form solution for the states
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function Y = LinearPendulum_Y(m,k,umax,x0,v0,tf,T)

% scaling constants
at = sqrt(m/k); % time
aq = umax/k; % position

% reshape and scale
T = T(:)/at;
x0 = x0/aq;
v0 = v0/aq*at;

% number of transitions
n = ceil(tf/at/pi);

% phase shift
p = mod(tf/at,pi);

% initialize
Y1c = cell(n,1); Y2c = Y1c; Tc = Y1c;

% go through each phase
for idx = 1:n
    % initial time
    i = max(p + (idx-2)*pi,0);

    % final time
    f = p + (idx-1)*pi;

    % current time vector
    t = T( i <= T & T <= f );

    % determine current sign of the control
    s = -sign(sin(tf/at))*(-1)^idx;

    % calculate states
    Y1c{idx} = LinearPendulum_Y1(s,i,x0,v0,t);
    Y2c{idx} = LinearPendulum_Y2(s,i,x0,v0,t);
    Tc{idx} = t;

    % store current initial states
    X0 = x0; V0 = v0;

    % update end phase state values
    x0 = LinearPendulum_Y1(s,i,X0,V0,f);
    v0 = LinearPendulum_Y2(s,i,X0,V0,f);
end

% combine
Yl = [vertcat(Y1c{:}),vertcat(Y2c{:})];
Tl = vertcat(Tc{:});

% output requested state values (because there might be duplicates)
Y = interp1(Tl,Yl,T);

% unscale
Y = Y*aq;
Y(:,2) = Y(:,2)/at;

end