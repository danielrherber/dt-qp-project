function [c,ceq] = JadduShimemura2_constraints(x)

% phase time knots
% t1 = x(1);
% y2_t1 = x(2);

% y2_t1 = x(1);
% OPTIONS = [];
% % t1 = fzero(@(t1) JadduShimemura2_U_t1_error(t1,y2_t1),0.99,OPTIONS);
% t1 = fminbnd(@(t1) abs(JadduShimemura2_U_t1_error(t1,y2_t1)),0.01,0.99,OPTIONS);

t1 = x(1);
OPTIONS = [];
y2_t1 = fzero(@(y2_t1) JadduShimemura2_U_t1_error(t1,y2_t1),0.99,OPTIONS);
% y2_t1 = fminbnd(@(y2_t1) abs(JadduShimemura2_U_t1_error(t1,y2_t1)),0.01,0.99,OPTIONS);

% phase 1 states
T1 = linspace(0,t1,1e6)';
Y_t1 = JadduShimemura2_Y_t1(T1,t1,y2_t1); 

% phase 2 states 
T2 = linspace(t1,1,1e6)';
Y_t2 = JadduShimemura2_Y_t2(T2,t1,y2_t1);

% extract state 1
y2_t1 = Y_t1(:,1);
y2_t2 = Y_t2(:,1);

% path constraint
c1 = y2_t1 - 8*(T1-0.5).^2 + 0.5;
c2 = y2_t2 - 8*(T2-0.5).^2 + 0.5;

c = max([c1;c2]);
ceq = [];

end