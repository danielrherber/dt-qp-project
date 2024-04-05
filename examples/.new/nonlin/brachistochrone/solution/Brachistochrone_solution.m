%--------------------------------------------------------------------------
% Brachistochrone_solution.m
% Create solution for Brachistochrone example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [tf,Y,U] = Brachistochrone_solution(t,casenum,G,XF,varargin)

switch casenum
%--------------------------------------------------------------------------
case {1,2}
% extract
YF = varargin{1};

% number of time points
Nt = 10000;

% initialize symbolic variables
syms theta k xf yf

% define symbolic functions
eqnX = symfun( k^2/2*(theta - sin(theta)) == xf , xf);
eqnY = symfun( k^2/2*(1 - cos(theta)) == yf , yf );

% solve for thetaf and k
S = vpasolve([eqnX(XF),eqnY(YF)],[theta,k],[0 2*pi; 0 Inf]);

% extract
thetaf = S.theta; k = S.k;

% scaling constant
at = sqrt(k^2/(2*G));

% final time
tf = double(at*thetaf);

% convert to double
thetaf = double(thetaf); k = double(k);

% theta and time horizon
theta = linspace(0,thetaf,Nt)';
T = linspace(0,tf,Nt)';

% states
Xs = k^2/2*(theta - sin(theta));
Ys = k^2/2*(1 - cos(theta));
Vs = k*sqrt(2*G)*sin(theta/2);

% control
Us = theta/2;

% interpolate states
Y = interp1(T,[Xs,Ys,Vs],t,'spline');

% interpolate controls
U = interp1(T,Us,t,'spline');
%--------------------------------------------------------------------------
case 3
% final time
tf = sqrt(pi*XF/G);

% intermediate constant
w = sqrt(pi/4*G/XF);

% control
U = pi/2 - w*t;

% states
Xs = (2*XF/pi)*(w*t - sin(2*w*t)/2);
Ys = (2*XF/pi)*sin(w*t).^2;
Vs = (G*sin(t*w))/w;

% combine
Y = [Xs,Ys,Vs];
%--------------------------------------------------------------------------
case 4
% extract
theta = varargin{1};
h = varargin{2};

% final time
tf = sqrt(2/G*(XF + h*cot(theta))*(theta + cot(theta))) - sqrt(2*h/G*cot(theta)*(theta - pi/2 + cot(theta)));

% intermediate constants
w1 = sqrt(G/2*(theta - (pi/2) + cot(theta))/(h*cot(theta)));
w2 = sqrt(G/2*(theta+cot(theta))/(XF + h*cot(theta)));
t1 = (pi/2 - theta)/w1;
t2 = tf - theta/w2;

U = Brachistochrone_U(tf,t,t1,t2,theta,w1,w2);
Y = Brachistochrone_Y(tf,G,t,t1,t2,theta,w1,w2);

end