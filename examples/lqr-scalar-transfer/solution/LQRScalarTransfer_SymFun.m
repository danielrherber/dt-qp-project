%--------------------------------------------------------------------------
% LQRScalarTransfer_SymFun.m
% High-precision solution for LQRScalarTransfer problem
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [U,Y,F] = LQRScalarTransfer_SymFun(ain,bin,cin,din,fin,qin,rin,tin)

% % see LQRScalarTransfer_syms.m
% syms a b c d f q r t
% Y = str2sym('c*exp(-(t*(a^2*r^2 + q*b^2*r)^(1/2))/r) - (exp(-(t*(a^2*r^2 + q*b^2*r)^(1/2))/r)*(exp((2*t*(a^2*r^2 + q*b^2*r)^(1/2))/r) - 1)*(c - d*exp((f*(a^2*r^2 + q*b^2*r)^(1/2))/r)))/(exp((2*f*(a^2*r^2 + q*b^2*r)^(1/2))/r) - 1)');
% U = str2sym('- (exp(-(t*(r*(r*a^2 + q*b^2))^(1/2))/r)*(a*d*exp(((r*(r*a^2 + q*b^2))^(1/2)*(f + 2*t))/r) + a*c*exp((2*f*(r*(r*a^2 + q*b^2))^(1/2))/r) - a*d*exp((f*(r*(r*a^2 + q*b^2))^(1/2))/r) - a*c*exp((2*t*(r*(r*a^2 + q*b^2))^(1/2))/r)))/(b*(exp((2*f*(r*(r*a^2 + q*b^2))^(1/2))/r) - 1)) - (exp(-(t*(r*(r*a^2 + q*b^2))^(1/2))/r)*(c*exp((2*f*(r*(r*a^2 + q*b^2))^(1/2))/r)*(r*(r*a^2 + q*b^2))^(1/2) - d*exp((f*(r*(r*a^2 + q*b^2))^(1/2))/r)*(r*(r*a^2 + q*b^2))^(1/2) + c*exp((2*t*(r*(r*a^2 + q*b^2))^(1/2))/r)*(r*(r*a^2 + q*b^2))^(1/2) - d*exp(((r*(r*a^2 + q*b^2))^(1/2)*(f + 2*t))/r)*(r*(r*a^2 + q*b^2))^(1/2)))/(b*r*(exp((2*f*(r*(r*a^2 + q*b^2))^(1/2))/r) - 1))');
% F = str2sym('((a^2*r^2 + q*b^2*r)^(1/2)*(c^2 + d^2) + a*r*(c^2 - d^2))/b^2 + ((a^2*r^2 + q*b^2*r)^(1/2)*(2*c^2 - 4*exp((f*(a^2*r^2 + q*b^2*r)^(1/2))/r)*c*d + 2*d^2))/(b^2*(exp((2*f*(a^2*r^2 + q*b^2*r)^(1/2))/r) - 1))');

% load the symbolic variables and functions
load('LQRScalarTransferSymFun.mat');

% evaluate Y
Y = subs(Y,a,ain);
Y = subs(Y,b,bin);
Y = subs(Y,c,cin);
Y = subs(Y,d,din);
Y = subs(Y,f,fin);
Y = subs(Y,q,qin);
Y = subs(Y,r,rin);
Y = subs(Y,t,tin);
Y = double(vpa(Y,20));

% evaluate U
U = subs(U,a,ain);
U = subs(U,b,bin);
U = subs(U,c,cin);
U = subs(U,d,din);
U = subs(U,f,fin);
U = subs(U,q,qin);
U = subs(U,r,rin);
U = subs(U,t,tin);
U = double(vpa(U,20));

% evaluate F
F = subs(F,a,ain);
F = subs(F,b,bin);
F = subs(F,c,cin);
F = subs(F,d,din);
F = subs(F,f,fin);
F = subs(F,q,qin);
F = subs(F,r,rin);
F = double(vpa(F,20));

end