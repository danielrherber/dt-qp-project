%--------------------------------------------------------------------------
% MultiphaseParameter_solution.m
% Solution function for MultiphaseParameter example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [Y,P,F] = MultiphaseParameter_solution(in)

% extract
p = in.p; Y0 = p.Y0; E1 = p.E1; E2 = p.E2; t1 = p.t1; t2 = p.t2; t3 = p.t3;

% fminbnd options
options = optimset('Display','none','TolX',1e-15);
pmin = -100;
pmax = 100;

% find the minimum
[P,F] = fminbnd(@(P) objective(P,t1,t2,t3,Y0,E1,E2),pmin,pmax,options);

% initialize
Y2 = []; Y3 = [];

% phase 1
T1 = in(1).t;
Y1 = MultiphaseParameter_Y1(P,T1',Y0);

% phase 2
if length(in) > 1
	T2 = in(2).t;
    Y1T1 = Y1(:,end) - E1;
    Y2 = MultiphaseParameter_Y2(P,T2',Y1T1,t1);
end

% phase 3
if length(in) > 2
	T3 = in(3).t;
    Y2T2 = Y2(:,end) - E2;
    Y3 = MultiphaseParameter_Y3(0,T3',Y2T2,t2);
end

% combine
Y = [Y1,Y2,Y3]';

end
% objective function
function F = objective(P,t1,t2,t3,Y0,E1,E2)
% states at the end of phase 1
Y1T1 = MultiphaseParameter_Y1(P,t1,Y0) - E1;

% states at the end of phase 2
Y2T2 = MultiphaseParameter_Y2(P,t2,Y1T1,t1) - E2;

% objective function
F = MultiphaseParameter_F1(P,t1,Y0) + ...
    MultiphaseParameter_F2(P,t1,t2,Y1T1) + ...
	MultiphaseParameter_F3(0,t2,t3,Y2T2);

end