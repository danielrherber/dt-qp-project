%--------------------------------------------------------------------------
% NeuenhofenKerriganX1_solution.m
% Creates solution for NeuenhofenKerriganX1 example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [T,U,Y,F] = NeuenhofenKerriganX1_solution(varargin)

% check if mesh was provided
if isempty(varargin)
    T = linspace(0,pi/2,1e4);
else
    T = varargin{1};
end

% control
U = 1/2 - 1.5./((cos(T)-2).^2);

% run accurate simulation
options = odeset('RelTol',1e-12,'AbsTol',eps);
[~,Yode] = ode15s(@(t,y) [y(1)^2/2 + 1/2 - 1.5/((cos(t)-2)^2); y(1)^2],T,[0;0],options);

% plot simulation
% plot(Tode,Yode)

% state
Y = Yode(:,1);

% control component of the objective function
Iu = -(pi*3^(1/2))/9;

% state component of the objective function
Iy = Yode(end,2);

% objective function
F = Iy + Iu;

end