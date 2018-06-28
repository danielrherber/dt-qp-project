close all; clear; clc

rng(20569180,'twister')

n = 6;

syms s t real
syms P real
syms tf0 tf1 tf2 tf3 real

y0 = sym('y0',[n 1]);

G = sym(randi([-1 1],n,1));

%% phase 1
tf0 = 0;

% A = sym(gallery('hanowa',n,-0.1));
A = sym(gallery('jordbloc',n,-1));

Y1t = expm(A*(t-tf0))*y0 + P*int(expm(A*(t-s))*G,s,tf0,t);
Y1t = simplify(Y1t,'steps',100);

F1 = int(sum(Y1t.^2),t,tf0,tf1);
F1 = simplify(F1,'steps',100);

funF1 = matlabFunction(F1,'Vars',{P,tf1,y0});
funY1t = matlabFunction(Y1t,'Vars',{P,t,y0});

%% phase 2
A = sym(gallery('hanowa',n,-0.05));

Y2t = expm(A*(t-tf1))*y0 + P*int(expm(A*(t-s))*G,s,tf1,t);
Y2t = simplify(Y2t,'steps',100);

F2 = int(sum(Y2t.^2),t,tf1,tf2);
F2 = simplify(F2,'steps',100);

funF2 = matlabFunction(F2,'Vars',{P,tf1,tf2,y0});
funY2t = matlabFunction(Y2t,'Vars',{P,t,y0,tf1});

%% phase 3
A = sym(-eye(n));

Y3t = expm(A*(t-tf2))*y0 + P*int(expm(A*(t-s))*G,s,tf2,t);
Y3t = simplify(Y3t,'steps',100);

F3 = int(sum(Y3t.^2),t,tf2,tf3);
F3 = simplify(F3,'steps',100);

funF3 = matlabFunction(F3,'Vars',{P,tf2,tf3,y0});
funY3t = matlabFunction(Y3t,'Vars',{P,t,y0,tf2});

%%
path = msavename(mfilename('fullpath'),'');

prob = 'MultiphaseParameter';

aF1 = matlabFunction(F1,'file',[path,prob,'_F1'],'Vars',{P,tf1,y0});
aF2 = matlabFunction(F2,'file',[path,prob,'_F2'],'Vars',{P,tf1,tf2,y0});
aF3 = matlabFunction(F3,'file',[path,prob,'_F3'],'Vars',{P,tf2,tf3,y0});

aY1t = matlabFunction(Y1t,'file',[path,prob,'_Y1'],'Vars',{P,t,y0});
aY2t = matlabFunction(Y2t,'file',[path,prob,'_Y2'],'Vars',{P,t,y0,tf1});
aY3t = matlabFunction(Y3t,'file',[path,prob,'_Y3'],'Vars',{P,t,y0,tf2});

return

%% plots
TF0 = 0;
TF1 = 5;
TF2 = 25;
TF3 = 30;

N = 10;

T1 = linspace(TF0,TF1,100*N);
T2 = linspace(TF1,TF2,100*N);
T3 = linspace(TF2,TF3,100*N);

% E1 = randi([-1 1],n,1);
% E2 = randi([-1 1],n,1);

E1 = [4;zeros(n-1,1)];
E2 = [0;zeros(n-1,1)];

% Pv = -0.5;
% Pv = 0;
Pv = -1.6;

Y0 = 10*ones(n,1);
Y0 = Y0(:);
Y0 = Y0(1:n);

figure
c = lines(n);

Y1 = funY1t(Pv,T1,Y0);
for k = 1:n
    plot(T1,Y1(k,:),'color',c(k,:)); hold on
end
Y1T1 = Y1(:,end) + E1;

Y2 = funY2t(Pv,T2,Y1T1,TF1);
for k = 1:n
    plot(T2,Y2(k,:),'color',c(k,:)); hold on
end
Y2T2 = Y2(:,end) + E2;

Y3 = funY3t(Pv,T3,Y2T2,TF2);
for k = 1:n
    plot(T3,Y3(k,:),'color',c(k,:)); hold on
end


Pv2 = linspace(-10,10,1000);
for k = 1:length(Pv2)
	Y1 = funY1t(Pv2(k),T1(end),Y0) + E1;
    Y2 = funY2t(Pv2(k),T3(1),Y1,T2(1)) + E2;
    F(k) = funI1(Pv2(k),TF1,Y0) + funI2(Pv2(k),TF1,TF2,Y1) + funI2(0,TF2,TF3,Y2);
end

figure
plot(Pv2,F)