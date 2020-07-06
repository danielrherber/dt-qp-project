%--------------------------------------------------------------------------
% DTQP_TEST_qlin_parameters.m
% A collection of test cases using qlin with parameters
% Comparisons made to SQP with fmincon
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
clc; clear; close all

tests = 1:4;
% tests = 4;

% go through the tests
for k = 1:length(tests)

    % fmincon options
    options = optimoptions('fmincon','Algorithm','sqp','Display','none',...
        'ScaleProblem',false);

    % test setup
    switch tests(k)
        %------------------------------------------------------------------
        case 1
        symb.Ob = '(p1-1)^4 + p1^2 + (p2-pi)^2';
        [P2,F2] = fmincon(@(x) (x(1)-1)^4 + x(1)^2 + (x(2)-pi)^2,[1,1],[],[],[],[],[],[],[],options);
        %------------------------------------------------------------------
        case 2
        symb.Ob = '(p1-1)^4 + log(p1+4)^2 + (p2-pi)^2';
        [P2,F2] = fmincon(@(x) (x(1)-1)^4 + log(x(1)+4)^2 + (x(2)-pi)^2,[1,1],[],[],[],[],[],[],[],options);
        %------------------------------------------------------------------
        case 3
        symb.Ob = '(p1-1)^4 + (p2-3)^2 + p1^2*p2';
        [P2,F2] = fmincon(@(x) (x(1)-1)^4 + (x(2)-3)^2 + x(1)^2*x(2),[1,1],[],[],[],[],[],[],[],options);
        %------------------------------------------------------------------
        case 4
        symb.Ob = '(p1-1)^2 + (p2-3)^2 + p1*p2';
        [P2,F2] = fmincon(@(x) (x(1)-1)^2 + (x(2)-3)^2 + x(1)*x(2),[1,1],[],[],[],[],[],[],[],options);
    end

    % problem structure
    [setup,opts] = problem(symb);

    % run the test and time
    [T,U,Y,P1,F1,in,opts] = DTQP_solve(setup,opts);

    % test analysis
    disp(strcat("Case ",string(tests(k))))
    disp("DTQP | F"); disp(F1); disp("P1"); disp(P1(:)')
    disp("FMINCON | F"); disp(F2); disp("P2"); disp(P2(:)')
end

% problem structure
function [setup,opts] = problem(symb)

n.ny = 1; % number of states
n.np = 2; % number of parameters

% system dynamics
symb.D = '[0]';

% time horizon
t0 = 0; tf = 1;

% combine structures
setup.symb = symb; setup.n = n;
setup.t0 = t0; setup.tf = tf; setup.p = [];

% options
opts.general.displevel = 0;
opts.general.plotflag = 0;
opts.dt.nt = 2;
opts.method.form = 'qlin';

end