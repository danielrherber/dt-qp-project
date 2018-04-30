%--------------------------------------------------------------------------
% LQRstandard_output.m
% Output function for LQRstandard example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [O,sol] = LQRstandard_output(T,U,Y,P,F,p,opts)

% compute matrices for the original example
matrices = SolutionMatrices;

% check if the solved problem is equivalent to the original example
if (isequal(matrices.A,p.A) && isequal(matrices.B,p.B) ...
        && isequal(matrices.R,p.R) && isequal(matrices.Q,p.Q) ...
        && isequal(matrices.M,p.M) && isequal(matrices.x0,p.x0)...
        && isequal(matrices.t0,p.t0) && isequal(matrices.tf,p.tf))

    % check if the solution exists
    if exist('LQRstandard_solution1.mat','file')
        
        % load the solution
        load('LQRstandard_solution1.mat','D');
        
    else % construct the solution for later use

        % display warning
        disp('warning: computing high accuracy solution for later use')
        disp('warning: this may take a few minutes but only needs to be done ONCE')
        
        % high accuracy solution
        opts.solmethod = 'bvp';
        opts.tolode = 1e-9;
        opts.tolbvp = 1e-7;
        D = LQRstandard_solution(p.A,p.B,p.R,p.Q,p.M,p,opts);
        
        % path for saving
        savepath = mfoldername('LQRstandard_solution','');
        
        % save the solution
        save(fullfile(savepath,'LQRstandard_solution1.mat'),'D');
        
    end
    
else % create the solution
    opts.solmethod = 'ode'; % 'bvp'
    opts.tolode = 1e-7; % 1e-11
    opts.tolbvp = 1e-3; % 1e-6
    D = LQRstandard_solution(p.A,p.B,p.R,p.Q,p.M,p,opts);
    
end

% interpolate
sol(1).T = T;
sol(1).U = interp1(D.T,D.U,T,'PCHIP');
sol(1).Y = interp1(D.T,D.Y,T,'PCHIP');
sol(1).F = D.F;

% solution on high resolution T
if opts.general.plotflag
    sol(2).T = linspace(p.t0,p.tf,1e4)';
    sol(2).U = interp1(D.T,D.U,sol(2).T,'PCHIP');
    sol(2).Y = interp1(D.T,D.Y,sol(2).T,'PCHIP');
    sol(2).F = D.F;
end

% errors
errorU = abs(U-sol(1).U);
errorY = abs(Y-sol(1).Y);
errorF = abs(F-sol(1).F);

% outputs
O(1).value = max(errorY(:));
O(1).label = 'Ymax';
O(2).value = max(errorU(:));
O(2).label = 'Umax';
O(3).value = max(errorF);
O(3).label = 'F';
O(4).value = max(opts.QPcreatetime);
O(4).label = 'QPcreatetime';
O(5).value = max(opts.QPsolvetime);
O(5).label = 'QPsolvetime';
end
% compute matrices for the original example
function matrices = SolutionMatrices
    matrices.ns = 20; % number of states
    matrices.nu = 10; % number of controls
    matrices.t0 = 0; % time horizon
    matrices.tf = 10; % time horizon
    matrices.x0 = linspace(-5,5,matrices.ns)'; % initial states
    rng(393872382) % specific random seed
    matrices.A = sprand(matrices.ns,matrices.ns,0.5,1);
    matrices.B = sprand(matrices.ns,matrices.nu,1,1);
    matrices.R = eye(matrices.nu);
    matrices.Q = sprand(matrices.ns,matrices.ns,0.2);
    matrices.Q = ((matrices.Q)*((matrices.Q)'))/100;
    matrices.M = 10*eye(matrices.ns); % objective
end