syms t real

u = -12*(t - 4)^2;
v = -4*(t - 4)^3;
x = 15 - (t - 4)^4;

L = x^2 + 1e-3*u^2;
I = int(L,t,34/15,4);

%%
path = msavename(mfilename('fullpath'),'');

prob = 'BettsBiehnCampbell1';

aU = matlabFunction(u,'file',[path,prob,'_U']);
aX = matlabFunction([x,v],'file',[path,prob,'_Y']);
aPSI = matlabFunction(I,'file',[path,prob,'_F']);