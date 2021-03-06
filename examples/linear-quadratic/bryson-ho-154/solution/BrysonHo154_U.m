function u = BrysonHo154_U(t,c1,c2,t0,tf,v0,x0)
%BRYSONHO154_U
%    U = BRYSONHO154_U(T,C1,C2,T0,TF,V0,X0)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    19-Dec-2017 12:06:00

t2 = c1.*c2.*x0.*6.0;
t3 = c1.*c2.*t.*v0.*6.0;
t4 = c2.*v0.*1.2e1;
t5 = c2.*t.*v0.*1.2e1;
t6 = c1.*c2.*t.*x0.*1.2e1;
t7 = tf.^2;
t8 = t0.^2;
u = -(-t8.*(t2+t3)+c1.*v0.*1.2e1-tf.*(t5+t6-c2.*x0.*1.2e1)+t0.*(t5+t6-tf.*(t4-c1.*c2.*t.*v0.*1.2e1)-c1.*c2.*t7.*v0.*6.0)+t7.*(t2-t3+t4)-c2.*t.*x0.*1.2e1+c1.*c2.*t0.*t8.*v0.*2.0+c1.*c2.*t7.*tf.*v0.*4.0)./(c1.*t0.*-1.2e1+c1.*tf.*1.2e1-c2.*t0.*t7.*1.2e1-c2.*t0.*t8.*4.0+c2.*t7.*tf.*4.0+c2.*t8.*tf.*1.2e1+c1.*c2.*t7.^2+c1.*c2.*t8.^2+c1.*c2.*t7.*t8.*6.0-c1.*c2.*t0.*t7.*tf.*4.0-c1.*c2.*t0.*t8.*tf.*4.0+1.2e1);
