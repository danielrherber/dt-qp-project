function Y3t = MultiphaseParameter_Y3(P,t,in3,tf2)
%MULTIPHASEPARAMETER_Y3
%    Y3T = MULTIPHASEPARAMETER_Y3(P,T,IN3,TF2)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    10-May-2018 23:44:59

y01 = in3(1,:);
y02 = in3(2,:);
y03 = in3(3,:);
y04 = in3(4,:);
y05 = in3(5,:);
y06 = in3(6,:);
t2 = -t+tf2;
t3 = exp(t2);
t4 = t3-1.0;
t5 = P.*t4;
Y3t = [-P.*t4+t3.*y01;t3.*y02;-P.*t4+t3.*y03;t5+t3.*y04;-t5+t3.*y05;-t5+t3.*y06];