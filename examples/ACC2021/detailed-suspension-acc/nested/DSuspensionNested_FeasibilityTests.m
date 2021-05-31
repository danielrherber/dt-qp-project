%--------------------------------------------------------------------------
% DSuspensionNested_FeasibilityTests.m
% Feasibility tests for the nested method
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clear; clc

% number of random plants to test
Nrand = 1e4;

% QP to test feasible inner-loop problem?
QPflag = false;

% suspension problem parameters
p = DSuspensionProblem_Parameters;

% inner-loop options
p.nt = 200; % number of time points
p.InnerLoopTolerance = 1e-6; % inner-loop tolerance

% extract
UB = p.xpUpperBound; LB = p.xpLowerBound;

% number of plant design variables
nxp = length(LB);

% create linear plant constraints
[A,b] = DSuspensionNested_LinearDesignConstraints(p);

%% run the tests
t1 = tic;

% generate the random points
Xrand = (UB-LB).*rand(Nrand,nxp) + LB;

% initialize
I1 = false(Nrand,1); I2 = I1; I0 = I1;

% go through the random plant designs
parfor k = 1:Nrand

    % check if linear plant constraints are satisfied
    I0(k) = all(A*(Xrand(k,:)') - b <= 0);

    % check if nonlinear plant constraints are satisfied
    [c,ceq] = DSuspensionNested_DesignConstraints(Xrand(k,:),p);
    I1(k) = all(c<=0);

    % determine inner-loop feasibility
    if QPflag

        % solve inner-loop problem with ramp
        [~,out1] = DSuspensionNested_InnerLoopRamp(Xrand(k,:)',p);

        % solve inner-loop problem with rough road
        [~,out2] = DSuspensionNested_InnerLoopRough(Xrand(k,:)',p);

        % determine if inner-loop problem was feasible
        feasible = ~(isnan(out1.F) || isnan(out2.F));
        I2(k) = feasible;

    else

        % alternative test
        feasible = DSuspensionNested_InnerLoopFeasible(Xrand(k,:),p);
        I2(k) = feasible;

    end

    % display every 1e4
    if mod(k,1e4) == 0
       disp(string(k))
    end

end
toc(t1)

% display to command window
disp(strcat(string(sum(I0)/Nrand*100),"% feasible"))
disp(strcat(string(sum(I1)/Nrand*100),"% feasible"))
disp(strcat(string(sum(I0.*I1)/Nrand*100),"% feasible"))
disp(strcat(string(sum(I2)/Nrand*100),"% feasible"))
disp(strcat(string(sum(I0.*I1.*I2)/Nrand*100),"% feasible"))

%% feasibility plots
close all

% check if we have too many points
if Nrand >= 1e5
    return
end

% feasibility region to plot
I = I2;

 % colors
cFeasible =[1 0 0];
cInfeasible = [0.8 0.8 0.8];

% assign infeasible/feasible colors
c = repmat(cInfeasible,Nrand,1);
c(I,1) = cFeasible(1); c(I,2) = cFeasible(2); c(I,3) = cFeasible(3);

% initialize figure
hf = figure; hold on
hf.Color = 'w';

% go through each variable combination
for ix = 1:nxp
    for jx = 1:nxp

        % assign subplot
        subplot(nxp,nxp,sub2ind([nxp nxp],ix,jx))

        % create scatter plot
        scatter(Xrand(:,ix),Xrand(:,jx),6,c)

        axis tight

    end
end