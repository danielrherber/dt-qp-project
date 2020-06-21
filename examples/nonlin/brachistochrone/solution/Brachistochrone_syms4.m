%--------------------------------------------------------------------------
% Brachistochrone_syms4.m
% Creates solution for Brachistochrone example for case number 4
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clear; clc

syms t f w1 w2 t1 t2 theta g a real
syms x(t) y(t) v(t)

% control
U1 = pi/2 - w1*t;
U2 = theta;
U3 = w2*(f - t);

% states
V1 = dsolve( diff(v,t) == g*sin(U1), v(0) == 0);
V2 = dsolve( diff(v,t) == g*sin(U2), v(t1) == subs(V1,t,t1) );
V3 = dsolve( diff(v,t) == g*sin(U3), v(t2) == subs(V2,t,t2) );

X1 = dsolve( diff(x,t) == V1*cos(U1), x(0) == 0);
X2 = dsolve( diff(x,t) == V2*cos(U2), x(t1) == subs(X1,t,t1) );
X3 = dsolve( diff(x,t) == V3*cos(U3), x(t2) == subs(X2,t,t2) );

Y1 = dsolve( diff(y,t) == V1*sin(U1), y(0) == 0);
Y2 = dsolve( diff(y,t) == V2*sin(U2), y(t1) == subs(Y1,t,t1) );
Y3 = dsolve( diff(y,t) == V3*sin(U3), y(t2) == subs(Y2,t,t2) );

% piecewise functions
U1 = heaviside(t1 - t)*U1;
U2 = heaviside(t2 - t)*heaviside(t - t1)*U2;
U3 = heaviside(t - t2)*U3;

V1 = heaviside(t1 - t)*V1;
V2 = heaviside(t2 - t)*heaviside(t - t1)*V2;
V3 = heaviside(t - t2)*V3;

X1 = heaviside(t1 - t)*X1;
X2 = heaviside(t2 - t)*heaviside(t - t1)*X2;
X3 = heaviside(t - t2)*X3;

Y1 = heaviside(t1 - t)*Y1;
Y2 = heaviside(t2 - t)*heaviside(t - t1)*Y2;
Y3 = heaviside(t - t2)*Y3;

% current file path
path = msavename(mfilename('fullpath'),'');

% problem name
prob = 'Brachistochrone';

% create matlab functions
aU = matlabFunction(U1+U2+U3,'file',[path,prob,'_U']);
aY = matlabFunction([X1+X2+X3,Y1+Y2+Y3,V1+V2+V3],'file',[path,prob,'_Y']);

return

%% example
THETA = atan(0.5);
G = 32.1740;
H = 0.1;
XF = 1;

TF = sqrt(2/G*(XF + H*cot(THETA))*(THETA + cot(THETA))) - sqrt(2*H/G*cot(THETA)*(THETA - pi/2 + cot(THETA)));
W1 = sqrt(G/2*(THETA - (pi/2) + cot(THETA))/(H*cot(THETA)));
W2 = sqrt(G/2*(THETA+cot(THETA))/(XF + H*cot(THETA)));
T1 = (pi/2 - THETA)/W1;
T2 = TF - THETA/W2;

T = linspace(0,TF,10000)';

U = Brachistochrone_U(TF,T,T1,T2,THETA,W1,W2);
Y = Brachistochrone_Y(TF,G,T,T1,T2,THETA,W1,W2);

figure; hold on
plot(T,U,'.-')

figure; hold on
plot(T,Y,'.-')