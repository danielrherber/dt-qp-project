close all; clear; clc

% initial options
path = msavename(mfilename('fullpath'),'');
prob = 'JadduShimemura3';
simplifyopts = {'Steps',100,'criterion','preferReal','Seconds',10};
pooldata = gcp('nocreate');
phase = [1,2];
% phase = [];

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

%% phase 1: constrain active at the end
if any(phase==1)
tic

% solve the boundary value problem
sol_1 = dsolve( eq1, eq2, ...
    Y1(0) == 0, Y2(0) == -1, ...
    Y1(0.5) == 0.5, Y2(0.5) == y2_t1);

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

%% phase 2: constraint inactive
if any(phase==2)
tic

% solve the boundary value problem
sol_t2 = dsolve( eq1, eq2, ...
    L1(1) == 0, L2(1) == 0, ...
    Y1(0.5) == 0.5, Y2(0.5) == y2_t1);

% extract and simplify
if ~isempty(pooldata)
    spmd
        switch labindex
            case 1
                Y1_t2 = sol_t2.Y1; % Y1_t2 = simplify(Y1_t2,simplifyopts{:});
            case 2
                Y2_t2 = sol_t2.Y2; % Y2_t2 = simplify(Y2_t2,simplifyopts{:});
            case 3
                L1_t2 = sol_t2.L1; % L1_t2 = simplify(L1_t2,simplifyopts{:});
            case 4
                L2_t2 = sol_t2.L2; % L2_t2 = simplify(L2_t2,simplifyopts{:});
            otherwise
                % do nothing
        end
    end
    
    % extract results
    Y1_t2 = Y1_t2{1};
    Y2_t2 = Y2_t2{2};
    L1_t2 = L1_t2{3};
    L2_t2 = L2_t2{4};
    
else
    Y1_t2 = sol_t2.Y1; % Y1_t2 = simplify(Y1_t2,simplifyopts{:});
    Y2_t2 = sol_t2.Y2; % Y2_t2 = simplify(Y2_t2,simplifyopts{:});
    L1_t2 = sol_t2.L1; % L1_t2 = simplify(L1_t2,simplifyopts{:});
    L2_t2 = sol_t2.L2; % L2_t2 = simplify(L2_t2,simplifyopts{:});
end

% calculate the control
U_t2 = -inv(R)*B'*[L1_t2;L2_t2];
% U_t2 = simplify(U_t2,simplifyopts{:});

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

return

%% solve the three variable optimization problem

% A = []; b = [];
% LB = [0.25,-2];
% UB = [0.75,2];
% X0 = [0.45,-1];

% A = []; b = [];
% LB = [-2];
% UB = [2];
% X0 = [-0.5];

A = []; b = [];
LB = [-2];
UB = [2];
X0 = [1];

% options = optimoptions('ga','Display','Iter','PopulationSize',10,...
%     'UseParallel',true,'FunctionTolerance',1e-12);
% X = ga(@(x) JadduShimemura2_objective(x),1,A,b,[],[],LB,UB,...
%     @(x) JadduShimemura2_constraints(x),options)

options = optimoptions('patternsearch','Display','Iter','UseParallel',false,...
    'MeshTolerance',eps,'UseCompletePoll',true,'UseCompleteSearch',true,...
    'MaxFunctionEvaluations',1e6);
X = patternsearch(@(x) JadduShimemura3_objective(x),X0,A,b,[],[],LB,UB,...
    [],options) 

X0 = X;
options = optimoptions('fmincon','Display','Iter','UseParallel',true,...
    'Algorithm','interior-point','OptimalityTolerance',1e-15,...
    'FunctionTolerance',1e-15,'StepTolerance',1e-15);
X = fmincon(@(x) JadduShimemura3_objective(x),X0,A,b,[],[],LB,UB,...
    [],options)

F = JadduShimemura2_objective(X);

%%
close all

t1 = 0.5;
y2_t1 = X(1);
aP = matlabFunction(sym(y2_t1),'file',[path,prob,'_parameters'],'Optimize',true);

T1 = linspace(0,t1,1e5)';
T2 = linspace(t1,1,1e5)';

Y_1_T1 = JadduShimemura3_Y_t1(T1,y2_t1);
Y_2_T2 = JadduShimemura3_Y_t2(T2,y2_t1);

U_T1 = JadduShimemura3_U_t1(T1,y2_t1);
U_T2 = JadduShimemura3_U_t2(T2,y2_t1);

figure; hold on
t22 = linspace(0.15,0.85,1000);
plot(T1,Y_1_T1,'.','linewidth',1);
plot(T2,Y_2_T2,'.','linewidth',1);

legend('constraint','Y1','Y2')

figure; hold on
plot(T1,U_T1);
plot(T2,U_T2);

% F = JadduShimemura2_objective([t1,t2,y2_t1]);
% disp(F)

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
