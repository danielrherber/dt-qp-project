%--------------------------------------------------------------------------
% NeuenhofenKerriganX2_solution.m
% Creates solution for NeuenhofenKerriganX2 example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [T,U,Y,F] = NeuenhofenKerriganX2_solution(varargin)

% check if mesh was provided
if isempty(varargin)
    T = linspace(0,1,1e4)';
else
    T = varargin{1};
end

% initialize
U = nan(length(T),1);
Y = nan(length(T),2);

% transition times
t1 = 1-sqrt(41)/10;
t2 = t1 + log(2) - log(sqrt(41)-5)/2;

% ode options
options = odeset('RelTol',1e-12,'AbsTol',eps);

%--- phase 1
% indices for this phase
I1 = T<t1;

% user time mesh for this phase
T1 = T(I1);

% control
U(I1) = -1;

% ode mesh should include the endpoints
T1ode = unique([T1;t1]);

% run accurate simulation (phase 1)
[T1ode,Y1ode] = ode15s(@(t,y) [-0.5/y(1); 4*y(1)^4 + 1],T1ode,[1;0],options);

% interpolate the solution (should be exact to ode solution)
Y(I1,:) = interp1(T1ode,Y1ode,T1);

%--- phase 2
% indices for this phase
I2 = (T>=t1) & (T<t2);

% user time mesh for this phase
T2 = T(I2);

% control
U(I2) = 0.8*sinh(2*(T2-t2));

% ode mesh should include the endpoints
T2ode = unique([t1;T2;t2]);

% run accurate simulation (phase 1)
[T2ode,Y2ode] = ode15s(@(t,y) [0.5*0.8*sinh(2*(t-t2))/y(1); 4*y(1)^4 + (0.8*sinh(2*(t-t2)))^2],T2ode,Y1ode(end,:),options);

% interpolate the solution (should be exact to ode solution)
Y(I2,:) = interp1(T2ode,Y2ode,T2);

%--- phase 3
% indices for this phase
I3 = T>=t2;

% user time mesh for this phase
T3 = T(I3);

% control
U(I3) = 0;

% ode mesh should include the endpoints
T3ode = unique([t2;T3]);

% run accurate simulation (phase 1)
[T3ode,Y3ode] = ode15s(@(t,y) [0; 4*y(1)^4],T3ode,Y2ode(end,:),options);

% interpolate the solution (should be exact to ode solution)
Y(I3,:) = interp1(T3ode,Y3ode,T3);

%--- objective function
F = Y3ode(end,2);

%--- plot solution
% figure; hold on
% plot(T,U);
% plot(T,Y)

end