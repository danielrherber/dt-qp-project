%--------------------------------------------------------------------------
% BrysonHo63_solution.m
% Create solution for BrysonHo63 example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [F,E] = BrysonHo63_solution(tf,V,w,T)

% eccentricity
e = w/V;

% determine semi-major axis length
a = tf*sqrt(V^2 - w^2)/(2*pi);

% determine semi-minor axis length
b = sqrt(1 - e^2)*a;

% ellipse
tau = linspace(pi,-pi,length(T))';
E(:,1) = -b*sin(tau);
E(:,2) = a*cos(tau) + a;

% compute F
F = V^2*tf^2/(4*pi)*(1 - w^2/V^2)^(3/2);

end