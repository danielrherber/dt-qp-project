clear; clc; close all

% initialize
syms t real
syms em um a S f t0 y0 real
syms eta real

% intermediate variables
Q = em^(-2);
R = um^(-2);
Pss = R*(-1/a + (Q/R + 1/a^2)^(1/2));
Ab = -(Q/R+1/a^2)^(1/2);
Zss = -(2*R*(Q/R+1/a^2)^(1/2))^(-1);
T = log(-(((S-Pss)*Zss)^(-1)-1))/2;
Z = -2*Zss*exp(Ab*(t-f)+T)*sinh(Ab*(t-f)+T);
P = Pss + 1/Z;
cf = t-f;
c0 = t0-f;

% state
y = y0*sinh(Ab*cf+T)/sinh(Ab*c0+T) ...
    +(eta*exp(-T)/(1-Pss/S))...
    *(cosh(Ab*cf+T) - sinh(Ab*cf+T)*coth(Ab*c0+T)) ...
    +(Q*eta/R/Ab^2)*(1 - cosh(Ab*cf+T)*cosh(T) ...
    -sinh(Ab*cf+T)/sinh(Ab*c0+T) ...
    +sinh(Ab*cf+T)*coth(Ab*c0+T)*cosh(T));

% control
u = diff(y,t)+y/a;

%%
xpath = msavename(mfilename('fullpath'),'');
prob = 'TurnerChunJuang1';

aY = matlabFunction(y,'file',[xpath,prob,'_Y'],'Optimize',true);
aU = matlabFunction(u,'file',[xpath,prob,'_U'],'Optimize',true);

return

%%
p.t0 = 0; 
p.f = 10;
p.y0 = -5;
p.em = 1;
p.um = 1;
p.a = 1;
p.S = 5;
p.eta = 10;

tx = linspace(p.t0,p.f,100);

Y = TurnerChunJuang1_Y(p.S,p.a,p.em,p.eta,p.f,tx,p.t0,p.um,p.y0);
U = TurnerChunJuang1_U(p.S,p.a,p.em,p.eta,p.f,tx,p.t0,p.um,p.y0);
F = TurnerChunJuang1_F(p.S,p.a,p.em,p.eta,p.f,p.t0,p.um,p.y0);

Q = p.em^(-2);
R = p.um^(-2);
Pss = R*(-1/p.a + (Q/R + 1/p.a^2)^(1/2));
Ab = -(Q/R+1/p.a^2)^(1/2);
Zss = -(2*R*(Q/R+1/p.a^2)^(1/2))^(-1);
T = log(-(((p.S-Pss)*Zss)^(-1)-1))/2;

% conditions
disp(p.S-Pss)
disp(0)
disp(Zss)

% symbolic method
Ysyms = y;
Ysyms = subs(Ysyms,t0,p.t0);
Ysyms = subs(Ysyms,f,p.f);
Ysyms = subs(Ysyms,y0,p.y0);
Ysyms = subs(Ysyms,em,p.em);
Ysyms = subs(Ysyms,um,p.um);
Ysyms = subs(Ysyms,a,p.a);
Ysyms = subs(Ysyms,S,p.S);
Ysyms = subs(Ysyms,a,p.a);
Ysyms = subs(Ysyms,eta,p.eta);
Ysyms = subs(Ysyms,t,tx);

% evaluate symbolic
tic
Ysyms = double(vpa(Ysyms));
toc

% plot
figure; hold on
plot(tx,Ysyms);
plot(tx,Y);
plot(tx,U);