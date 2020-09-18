%-------------------------------------------------------------------------------
% This file provides some examples illustrating how to implement different
% linear-quadratic dynamic optimization (LQDO) problem elements using the
% structure-based syntax of the DT QP Project
%-------------------------------------------------------------------------------
% Examples Link: https://github.com/danielrherber/dt-qp-project/examples
% Project Link: https://github.com/danielrherber/dt-qp-project
%-------------------------------------------------------------------------------
% Comments
% All items should be put into a setup structure
% k is an appropriate index in the examples below, need to be incremental
%-------------------------------------------------------------------------------
% Time horizon and auxiliary data (setup.p = p)
% p.t0 = 0; p.tf = 20;
% p.a = 1;
%-------------------------------------------------------------------------------
% Dynamics
% setup.A = [0,1;-3,-2]; % ns = 2
% setup.A = {@(t) sin(t)}; % ns = 1
% setup.B = [1,2]; % ns = 1, nu = 2
% setup.B = {@(t) cos(t); @(t) sin(t)}; % ns = 2, nu = 1
% setup.d = {1;0; @(t) sin(t)}; % ns = 3
% setup.d = -1; % ns = 1
%-------------------------------------------------------------------------------
% Lagrange terms (setup.L = L)
% L(k).left = 1; L(k).right = 1; L(k).matrix = {@(t) sin(t)}; % nu = 1
% L(k).left = 4; L(k).right = 2; L(k).matrix = [0,0;1,0]; % ns = 2
% L(k).left = 0; L(k).right = 2; L(k).matrix = {@(t) exp(-t), 0, 1}; % ns = 3
% L(k).left = 0; L(k).right = 0; L(k).matrix = {@(t) sin(t)};
%-------------------------------------------------------------------------------
% Mayer terms (setup.M = M)
% M(k).left = 3; M(k).right = 3; M(k).matrix = 1; % np = 1
% M(k).left = 5; M(k).right = 5; M(k).matrix = eye(2); % ns = 2
% M(k).left = 0; M(k).right = 5; M(k).matrix = [1;-2]; % ns = 2
%-------------------------------------------------------------------------------
% Additional linear equality constraints (setup.Y = Y)
% # Y(k).linear(1).right = 1; Y(k).linear(1).matrix = 1;
% # Y(k).linear(2).right = 3; Y(k).linear(2).matrix = -1;
% # Y(k).b = @(t) sin(t); % nu = 1, np = 1
% & Y(k).linear(1).right = 5; Y(k).linear(1).matrix =[0;1]; % ns = 2
% & Y(k).b = 1;
%-------------------------------------------------------------------------------
% Additional linear inequality constraints (setup.Z = Z)
% # Z(k).linear(1).right = 1; Z(k).linear(1).matrix = 1;
% # Z(k).linear(2).right = 3; Z(k).linear(2).matrix = -1;
% # Z(k).b = @(t) sin(t); % nu = 1, np = 1
% & Z(k).linear(1).right = 2; Z(k).linear(1).matrix = [1;-2];
% & Z(k).b = 0; % ns = 2
%-------------------------------------------------------------------------------
% Simple upper bounds (setup.UB = UB)
% UB(k).right = 2; UB(k).matrix = [inf;pi]; % ns = 2
% UB(k).right = 4; UB(k).matrix = [inf;0;inf]; % ns = 3
%-------------------------------------------------------------------------------
% Simple lower bounds (setup.LB = LB)
% LB(k).right = 1; LB(k).matrix = {@(t) sin(t);-inf}; % nu = 2
% LB(k).right = 4; LB(k).matrix = [-inf;0;-inf]; % ns = 3
%-------------------------------------------------------------------------------
disp('type "help DTQP_template" in the command window')