function out1 = BrysonHo248_Y(a,b,c,t,t1,t2,tf)
%BRYSONHO248_Y
%    OUT1 = BRYSONHO248_Y(A,B,C,T,T1,T2,TF)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    11-Apr-2018 21:31:56

t4 = c.*t2;
t5 = -t1+t2;
t6 = exp(t5);
t7 = t1-t2;
t8 = exp(t7);
t9 = t1.^2;
t10 = t+2.0;
t11 = c.*t.*t10.*(1.0./2.0);
t12 = t-t1;
t13 = heaviside(t12);
t14 = t-t2;
t15 = heaviside(t14);
t16 = a.*2.0;
t17 = b.*t1.*2.0;
t18 = c.*t1.*2.0;
t19 = exp(-t);
t20 = exp(t1);
t21 = t16+t17+t18-c.*t9;
t22 = t13-t15;
t23 = t-tf;
t24 = heaviside(t23);
t25 = t15-t24;
t26 = a.*t8.*(1.0./2.0);
t27 = b.*t1.*t8.*(1.0./2.0);
t28 = c.*t1.*t8.*(1.0./2.0);
t29 = c.*t1.*t6.*(1.0./2.0);
t30 = c.*t6.*t9.*(1.0./4.0);
t33 = c.*t;
t31 = t4+t26+t27+t28+t29+t30-t33-a.*t6.*(1.0./2.0)-b.*t6-b.*t1.*t6.*(1.0./2.0)-c.*t8.*t9.*(1.0./4.0);
t32 = t13-1.0;
out1 = [t25.*(t4-t11+a.*t8-t.*t31+c.*t2.^2.*(1.0./2.0)-a.*t2.*t6.*(1.0./2.0)+a.*t2.*t8.*(1.0./2.0)-b.*t2.*t6+b.*t1.*t8+c.*t1.*t8-c.*t8.*t9.*(1.0./2.0)-b.*t1.*t2.*t6.*(1.0./2.0)+b.*t1.*t2.*t8.*(1.0./2.0)+c.*t1.*t2.*t6.*(1.0./2.0)+c.*t1.*t2.*t8.*(1.0./2.0)+c.*t2.*t6.*t9.*(1.0./4.0)-c.*t2.*t8.*t9.*(1.0./4.0))-t32.*(a+t11+t.*(b-c.*t))+t19.*t20.*t21.*t22.*(1.0./2.0),-t25.*t31-t32.*(b-t33)+t22.*(exp(-t1).*exp(t).*(b.*4.0+t16+t17-t18-c.*t9).*(1.0./4.0)-t19.*t20.*t21.*(1.0./4.0))];