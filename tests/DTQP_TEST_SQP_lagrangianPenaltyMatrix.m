%--------------------------------------------------------------------------
% DTQP_TEST_SQP_lagrangianPenaltyMatrix.m
% Test function for DTQP_SQP_lagrangianPenaltyMatrix.m
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
clc; clear; close all

tests = [1:4];
% tests = 4;

% go through each test
for k = 1:length(tests)

    opts = []; setup = []; param = []; % initialize

    % determine which test to run
    switch tests(k)
        %----------------------------------------------------------------------
        case 1
            f = 'y1^3'; % state derivative function
            ny = 1; % number of states
            nu = 0; % number of inputs
            np = 0; % number of parameters
            opts.dt.nt = 3; % number of time points
            setup.t0 = 0; % time initial
            setup.tf = 1; % final time
        %----------------------------------------------------------------------
        case 2
            f = '[u1^2*y1^2]';
            ny = 1; % number of states
            nu = 1; % number of inputs
            np = 0; % number of parameters
            opts.dt.nt = 3; % number of time points
            setup.t0 = 0; % time initial
            setup.tf = 1; % final time
        %----------------------------------------------------------------------
        case 3
            f = '[y1^3;y2^2*y1;y3 + u1]';
            ny = 3; % number of states
            nu = 1; % number of inputs
            np = 0; % number of parameters
            opts.dt.nt = 3; % number of time points
            setup.t0 = 0; % time initial
            setup.tf = 1; % final time
        %----------------------------------------------------------------------
        case 4
            f = '[y1*y2*u1*u2;y1^2*y2^2*u1^2*u2^2]';
            ny = 2; % number of states
            nu = 2; % number of inputs
            np = 0; % number of parameters
            opts.dt.nt = 3; % number of time points
            setup.t0 = 0; % time initial
            setup.tf = 1; % final time
    end

    % dummy state variables
    setup.A = zeros(ny);
    setup.B = zeros(ny,nu);
    setup.G = zeros(ny,np);

    % get default options
    [setup,opts] = DTQP_default_opts(setup,opts);

    % initialize problem
    [setup,in] = DTQP_initialize(setup,opts.dt);

    % initialize some elements
    U = rand(opts.dt.nt,nu);
    Y = rand(opts.dt.nt,ny);
    P = repmat(rand(1,np),opts.dt.nt,1);
    lambda = ones(opts.dt.nt-1,ny);
    X = [U Y P];
    opts.lambda = lambda;

    % evaluate the elements of the SQP problem
    form = 3; D2flag = true;
    E = DTQP_qlin_taylor(f,3,in,D2flag);

    % go through each state derivative function
    D2matrix = cell(size(E.D2));
    for i = 1:length(E.D2)
        D2matrix{i} = DTQP_qlin_update4tmatrix(E.D2{i},in.t,X,param);
    end

    % create indices for Lagrangian penalty matrix
    [I,J,V] = DTQP_SQP_lagrangianPenaltyMatrix(D2matrix,in,opts);

    % display the sparsity pattern
    spy(sparse(I,J,1));

    % display sequences
    disp([I,J,V]);

end