function u = BrysonHo116_U(t,tf,v0,x0)
%BRYSONHO116_U
%    U = BRYSONHO116_U(T,TF,V0,X0)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    14-Dec-2017 14:51:08

t2 = tf.^2;
t3 = v0.^2;
t4 = t2-t3-x0.*4.0-tf.*v0.*2.0;
t5 = sqrt(t4);
u = -heaviside(t-tf)-heaviside(t+1.0)+heaviside(t-t5.*(1.0./2.0)-tf.*(1.0./2.0)-v0.*(1.0./2.0))+heaviside(t+t5.*(1.0./2.0)-tf.*(1.0./2.0)-v0.*(1.0./2.0));
