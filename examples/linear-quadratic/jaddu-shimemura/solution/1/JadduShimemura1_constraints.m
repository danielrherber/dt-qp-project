function [c,ceq] = JadduShimemura1_constraints(x)

% phase time knots
% t1 = x(1);
% t2 = x(2);
% y1_t1 = x(3);

y1_t1 = x(1);
OPTIONS = [];
t1 = fzero(@(t1) JadduShimemura1_U_t1_error(t1,y1_t1),0.5,OPTIONS);
t2 = fzero(@(t2) JadduShimemura1_U_t2_error(t1,t2,y1_t1),0.75,OPTIONS);

%
T1 = linspace(0,t1,1e4)';
Y_t1 = JadduShimemura1_Y_t1(T1,t1,y1_t1); 

% 
T3 = linspace(t2,1,1e4)';
Y_t3 = JadduShimemura1_Y_t3(T3,t1,t2,y1_t1);

% 
y2_t1 = Y_t1(:,2);
y2_t3 = Y_t3(:,2);


c1 = y2_t1 - 8*(T1-0.5).^2 + 0.5;
c3 = y2_t3 - 8*(T3-0.5).^2 + 0.5;

c = [c1;c3];
ceq = [];

end