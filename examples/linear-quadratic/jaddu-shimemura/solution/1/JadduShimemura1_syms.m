close all; clear; clc

% initial options
path = msavename(mfilename('fullpath'),'');
prob = 'JadduShimemura1';
simplifyopts = {'Steps',100,'criterion','preferReal','Seconds',10};
pooldata = gcp('nocreate');
% phase = [1,2,3];
phase = [];

% initialize
syms t t1 t2 positive
syms Y1(t) Y2(t) 
syms L1(t) L2(t)
syms y1_t1 y2_t1 real
syms y1_t2 y2_t2 real
Y = [Y1;Y2];
dYdt = diff(Y,t);
L  = [L1;L2];
dLdt = diff(L,t);

% problem parameters
A = [0,1;0,-1];
B = [0;1];
Q = eye(2);
R = 0.005;
y0 = [0;-1];

% equations
u = -inv(R)*B'*L;
eq1 = dYdt == A*Y + B*u;
eq2 = -dLdt == Q*Y + A'*L;

%% phase 1: constraint inactive
if any(phase==1)
tic

% solve the boundary value problem
sol_1 = dsolve( eq1, eq2, ...
    Y1(0) == 0, Y2(0) == -1, ...
    Y1(t1) == y1_t1, Y2(t1) == 8*(t1-0.5)^2 - 0.5 );

% extract and simplify
if ~isempty(pooldata)
    spmd
        switch labindex
            case 1
                Y1_t1 = sol_1.Y1; % Y1_t1 = simplify(Y1_t1,simplifyopts{:});
            case 2
                Y2_t1 = sol_1.Y2; % Y2_t1 = simplify(Y2_t1,simplifyopts{:});
            case 3
                L1_t1 = sol_1.L1; % L1_t1 = simplify(L1_t1,simplifyopts{:});
            case 4
                L2_t1 = sol_1.L2; % L2_t1 = simplify(L2_t1,simplifyopts{:});
            otherwise
                % do nothing
        end
    end
    
    % extract results
    Y1_t1 = Y1_t1{1};
    Y2_t1 = Y2_t1{2};
    L1_t1 = L1_t1{3};
    L2_t1 = L2_t1{4};
    
else
    Y1_t1 = sol_1.Y1; % Y1_t1 = simplify(Y1_t1,simplifyopts{:});
    Y2_t1 = sol_1.Y2; % Y2_t1 = simplify(Y2_t1,simplifyopts{:});
    L1_t1 = sol_1.L1; % L1_t1 = simplify(L1_t1,simplifyopts{:});
    L2_t1 = sol_1.L2; % L2_t1 = simplify(L2_t1,simplifyopts{:});
end

% calculate the control
U_t1 = -inv(R)*B'*[L1_t1;L2_t1];
% U_t1 = simplify(U_t1,simplifyopts{:});

% integrand
I_t1 = Y1_t1^2 + Y2_t1^2 + R*U_t1^2;
% I_t1 = simplify(I_t1,simplifyopts{:});

% save the functions
% aI_t1 = matlabFunction(I_t1,'file',[path,prob,'_I_t1']);
% aY_t1 = matlabFunction([Y1_t1,Y2_t1],'file',[path,prob,'_Y_t1']);
% aU_t1 = matlabFunction(U_t1,'file',[path,prob,'_U_t1']);
aI_t1 = matlabFunction(I_t1,'file',[path,prob,'_I_t1'],'Optimize',false);
aY_t1 = matlabFunction([Y1_t1,Y2_t1],'file',[path,prob,'_Y_t1'],'Optimize',false);
aU_t1 = matlabFunction(U_t1,'file',[path,prob,'_U_t1'],'Optimize',false);

disp('phase 1 done')
toc

end
%% phase 2: constraint active
if any(phase==2)
tic

% controls and states
Y2_t2 = 8*(t-0.5)^2 - 0.5;
Y1_t2 = int(Y2_t2,t,t1,t) + y1_t1;
U_t2 = diff(Y2_t2,t) + Y2_t2;

% integrand
I_t2 = Y1_t2^2 + Y2_t2^2 + R*U_t2^2;
% I_t2 = simplify(I_t2,simplifyopts{:});

% save the functions
% aI_t2 = matlabFunction(I_t2,'file',[path,prob,'_I_t2']);
% aY_t2 = matlabFunction([Y1_t2,Y2_t2],'file',[path,prob,'_Y_t2']);
% aU_t2 = matlabFunction(U_t2,'file',[path,prob,'_U_t2']);
aI_t2 = matlabFunction(I_t2,'file',[path,prob,'_I_t2'],'Optimize',false);
aY_t2 = matlabFunction([Y1_t2,Y2_t2],'file',[path,prob,'_Y_t2'],'Optimize',false);
aU_t2 = matlabFunction(U_t2,'file',[path,prob,'_U_t2'],'Optimize',false);

disp('phase 2 done')
toc

end

%% phase 3: constraint inactive
if any(phase==3)
tic

y1_t2 = subs(Y1_t2,t,t2);

% solve the boundary value problem
sol_t3 = dsolve( eq1, eq2, ...
    L1(1) == 0, L2(1) == 0, ...
    Y1(t2) == y1_t2, Y2(t2) == 8*(t2-0.5)^2 - 0.5);

% extract and simplify
if ~isempty(pooldata)
    spmd
        switch labindex
            case 1
                Y1_t3 = sol_t3.Y1; % Y1_t3 = simplify(Y1_t3,simplifyopts{:});
            case 2
                Y2_t3 = sol_t3.Y2; % Y2_t3 = simplify(Y2_t3,simplifyopts{:});
            case 3
                L1_t3 = sol_t3.L1; % L1_t3 = simplify(L1_t3,simplifyopts{:});
            case 4
                L2_t3 = sol_t3.L2; % L2_t3 = simplify(L2_t3,simplifyopts{:});
            otherwise
                % do nothing
        end
    end
    
    % extract results
    Y1_t3 = Y1_t3{1};
    Y2_t3 = Y2_t3{2};
    L1_t3 = L1_t3{3};
    L2_t3 = L2_t3{4};
    
else
    Y1_t3 = sol_t3.Y1; % Y1_t3 = simplify(Y1_t3,simplifyopts{:});
    Y2_t3 = sol_t3.Y2; % Y2_t3 = simplify(Y2_t3,simplifyopts{:});
    L1_t3 = sol_t3.L1; % L1_t3 = simplify(L1_t3,simplifyopts{:});
    L2_t3 = sol_t3.L2; % L2_t3 = simplify(L2_t3,simplifyopts{:});
end

% calculate the control
U_t3 = -inv(R)*B'*[L1_t3;L2_t3];
% U_t3 = simplify(U_t3,simplifyopts{:});

% integrand
I_t3 = Y1_t3^2 + Y2_t3^2 + R*U_t3^2;
% I_t3 = simplify(I_t3,simplifyopts{:});

% save the functions
% aI_t3 = matlabFunction(I_t3,'file',[path,prob,'_I_t3']);
% aY_t3 = matlabFunction([Y1_t3,Y2_t3],'file',[path,prob,'_Y_t3']);
% aU_t3 = matlabFunction(U_t3,'file',[path,prob,'_U_t3']);
aI_t3 = matlabFunction(I_t3,'file',[path,prob,'_I_t3'],'Optimize',false);
aY_t3 = matlabFunction([Y1_t3,Y2_t3],'file',[path,prob,'_Y_t3'],'Optimize',false);
aU_t3 = matlabFunction(U_t3,'file',[path,prob,'_U_t3'],'Optimize',false);

disp('phase 3 done')
toc

end

% return

%% solve the three variable optimization problem

% A = []; b = [];
% LB = [0.01,0.5,-2];
% UB = [0.5,0.99,2];
% 
% % X0 = [0.49,0.51,0];
% X0 = [0.3,0.7,-0.0583];

A = []; b = [];
LB = [-2];
UB = [2];

% X0 = [0.51,0];
X0 = [-0.075];


% options = optimoptions('ga','Display','Iter','PopulationSize',10,...
%     'UseParallel',true,'FunctionTolerance',1e-12);
% X = ga(@(x) JadduShimemura1_objective(x),3,A,b,[],[],LB,UB,...
%     @(x) JadduShimemura1_constraints(x),options)

options = optimoptions('patternsearch','Display','Iter','UseParallel',false,...
    'MeshTolerance',100*eps,'UseCompletePoll',true,'UseCompleteSearch',true);
X = patternsearch(@(x) JadduShimemura1_objective(x),X0,A,b,[],[],LB,UB,...
    @(x) JadduShimemura1_constraints(x),options) 

% X0 = X;
% options = optimoptions('fmincon','Display','Iter','UseParallel',true,...
%     'Algorithm','interior-point');
% X = fmincon(@(x) JadduShimemura1_objective(x),X0,A,b,[],[],LB,UB,...
%     @(x) JadduShimemura1_constraints(x),options)

F = JadduShimemura1_objective(X);

%%
close all

% t1 = 0.3238;
% t2 = 0.7022;
% y1_t1 = -0.0743;

% t1 = X(1);
% t2 = X(2);
% y1_t1 = X(3);

y1_t1 = X(1);
OPTIONS = [];
t1 = fzero(@(t1) JadduShimemura1_U_t1_error(t1,y1_t1),0.5,OPTIONS);
t2 = fzero(@(t2) JadduShimemura1_U_t2_error(t1,t2,y1_t1),0.75,OPTIONS);

aP = matlabFunction(sym([t1,t2,y1_t1]),'file',[path,prob,'_parameters'],'Optimize',true);

T1 = linspace(0,t1,100)';
T2 = linspace(t1,t2,100)';
T3 = linspace(t2,1,100)';

Y_1_T1 = JadduShimemura1_Y_t1(T1,t1,y1_t1);
Y_2_T2 = JadduShimemura1_Y_t2(T2,t1,y1_t1);
Y_3_T3 = JadduShimemura1_Y_t3(T3,t1,t2,y1_t1);

U_T1 = JadduShimemura1_U_t1(T1,t1,y1_t1);
U_T2 = JadduShimemura1_U_t2(T2);
U_T3 = JadduShimemura1_U_t3(T3,t1,t2,y1_t1);

figure; hold on
t22 = linspace(0.15,0.85,1000);
plot(t22,8*(t22-0.5).^2 - 0.5,'linewidth',1,'color','k'); hold on
plot(T1,Y_1_T1,'linewidth',1);
plot(T2,Y_2_T2,'linewidth',1);
plot(T3,Y_3_T3,'linewidth',1);


legend('Y2')

figure; hold on
plot(T1,U_T1);
plot(T2,U_T2);
plot(T3,U_T3);

F = JadduShimemura1_objective([t1,t2,y1_t1]);
disp(F)

return


%%
aDY2_1 = matlabFunction(diff(Y2_t1,t));

y1_t1 = -0.4;

options = optimset('Display','iter');
t1 = fzero( @(x) aDY2_1(x,x,y1_t1) - (16*x-8), 0+0.2, options);



function E = zerofun(t1,y1_t1)

Y_t1 = JadduShimemura1_Y_t2(t1,t1,y1_t1);
y2_t1 = Y_t1(2);

Y1 = JadduShimemura1_Y_t1(t1,t1,y1_t1,y2_t1);
Y2 = JadduShimemura1_Y_t2(t1,t1,y1_t1);

E = Y1(1)-Y2(2);

end
