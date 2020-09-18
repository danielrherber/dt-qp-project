function F = JadduShimemura3_objective(x)

y2_t1 = x(1);

% calculate the integrals
I1 = integral(@(t) JadduShimemura3_I_t1(t,y2_t1),0,0.5,'AbsTol',eps,'RelTol',eps);
I2 = integral(@(t) JadduShimemura3_I_t2(t,y2_t1),0.5,1,'AbsTol',eps,'RelTol',eps);

% sum
F = I1 + I2;

end

