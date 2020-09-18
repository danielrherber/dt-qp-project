function F = JadduShimemura1_objective(x)

% phase time knots
% t1 = x(1);
% t2 = x(2);
% y1_t1 = x(3);

y1_t1 = x(1);
OPTIONS = [];
t1 = fzero(@(t1) JadduShimemura1_U_t1_error(t1,y1_t1),0.5,OPTIONS);
t2 = fzero(@(t2) JadduShimemura1_U_t2_error(t1,t2,y1_t1),0.75,OPTIONS);

% calculate the integrals
I1 = integral(@(t) JadduShimemura1_I_t1(t,t1,y1_t1),0,t1,'AbsTol',eps,'RelTol',eps);
I2 = integral(@(t) JadduShimemura1_I_t2(t,t1,y1_t1),t1,t2,'AbsTol',eps,'RelTol',eps);
I3 = integral(@(t) JadduShimemura1_I_t3(t,t1,t2,y1_t1),t2,1,'AbsTol',eps,'RelTol',eps);

% sum
F = I1 + I2 + I3;

end

