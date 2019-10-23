function F = JadduShimemura2_objective(x)

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

% calculate the integrals
I1 = integral(@(t) JadduShimemura2_I_t1(t,t1,y2_t1),0,t1,'AbsTol',eps,'RelTol',eps);
I2 = integral(@(t) JadduShimemura2_I_t2(t,t1,y2_t1),t1,1,'AbsTol',eps,'RelTol',eps);

% sum
F = I1 + I2;

end

