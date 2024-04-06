clear
clc

syms t real

p = sym(50);
x1 = sin(pi*t);
x2 = pi*cos(pi*t);
I = sym(0);

%%
path = msavename(mfilename('fullpath'),'');

prob = 'Biegler9p5';

aP = matlabFunction(p,'file',[path,prob,'_P']);
aX = matlabFunction([x1,x2],'file',[path,prob,'_Y']);
aPSI = matlabFunction(I,'file',[path,prob,'_F']);