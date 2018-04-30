%--------------------------------------------------------------------------
% OutputTracking_output.m
% Output function for OutputTracking example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [O,sol] = OutputTracking_output(T,U,Y,P,F,p,setup,opts)

% initialize
sol = [];

% combine
prob.setup = setup; prob.opts = opts;
prob.opts = rmfield(prob.opts,'QPcreatetime'); prob.opts = rmfield(prob.opts,'QPsolvetime');

% hash the input file to determine if the solution is available
filehash = DataHash(prob);

% hashed file name
filename = [opts.general.mname,'_solution_',filehash,'.mat'];
fullname = fullfile(opts.general.mpath,'solution',filename);

% check if the solution to the problem is already available
if exist(fullname,'file')

    % load the solution
    load(fullname,'D');

else % construct the solution for later use

    % display warning
    disp('warning: computing high accuracy solution for later use')
    disp('warning: this may take a few minutes but only needs to be done ONCE')

    % high accuracy solution
    opts.solmethod = 'bvp';
    opts.tolode = 1e-8;
    opts.tolbvp = 1e-8;

    D = OutputTracking_solution(p,opts);

    % save the solution
    save(fullname,'D');
        
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