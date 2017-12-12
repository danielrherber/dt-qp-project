%--------------------------------------------------------------------------
% BrysonHo109_output.m
% Output function for BrysonHo109 example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [O,sol] = BrysonHo109_output(T,U,Y,P,F,p,opts)

%--------------------------------------------------------------------------
if (p.a == 1) && (p.x0 == 1) && (p.tf == 1) && strcmp(func2str(p.g),'@(t)t.*cos(20*pi*t)-1/4')

    % check if the solution exists
    if exist('BrysonHo109_solution1.mat','file')
        
        % load the solution
        load('BrysonHo109_solution1.mat','D');
        
    else % construct the solution for later use
        % use high resolution time grid
        D = BrysonHo109_solution(linspace(0,p.tf,1e5)',Y,p);
        
        % path for saving
        savepath = mfoldername('BrysonHo109_solution','');
        
        % save the solution
        save(fullfile(savepath,'BrysonHo109_solution1.mat'),'D');
    end
    
    % interpolate
    sol(1).T = T;
    sol(1).U = interp1(D.T,D.U,T,'PCHIP');
    sol(1).Y = interp1(D.T,D.Y,T,'PCHIP');
    sol(1).F = D.F;
    
    % solution on high resolution T
    if opts.plotflag
        sol(2).T = linspace(p.t0,p.tf,1e4)';
        sol(2).U = interp1(D.T,D.U,sol(2).T,'PCHIP');
        sol(2).Y = interp1(D.T,D.Y,sol(2).T,'PCHIP');
        sol(2).F = D.F;
    end
%--------------------------------------------------------------------------
elseif (p.a == 2) && (p.x0 == 1) && (p.tf == 1) && strcmp(func2str(p.g),'@(t)t.*cos(20*pi*t)-1/4')

    % check if the solution exists
    if exist('BrysonHo109_solution2.mat','file')
        
        % load the solution
        load('BrysonHo109_solution2.mat','D');
        
    else % construct the solution for later use
        % use high resolution time grid
        D = BrysonHo109_solution(linspace(0,p.tf,1e5)',Y,p);
        
        % path for saving
        savepath = mfoldername('BrysonHo109_solution','');
        
        % save the solution
        save(fullfile(savepath,'BrysonHo109_solution2.mat'),'D');
    end
    
    % interpolate
    sol(1).T = T;
    sol(1).U = interp1(D.T,D.U,T,'PCHIP');
    sol(1).Y = interp1(D.T,D.Y,T,'PCHIP');
    sol(1).F = D.F;
    
    % solution on high resolution T
    if opts.plotflag
        sol(2).T = linspace(p.t0,p.tf,1e4)';
        sol(2).U = interp1(D.T,D.U,sol(2).T,'PCHIP');
        sol(2).Y = interp1(D.T,D.Y,sol(2).T,'PCHIP');
        sol(2).F = D.F;
    end
%--------------------------------------------------------------------------
else
    % solution on T
    sol = BrysonHo109_solution(T,Y,p);

    % solution on high resolution T
    if opts.plotflag
        D = BrysonHo109_solution(linspace(p.t0,p.tf,1e4)',Y,p);
        sol(2).T = D.T;
        sol(2).U = D.U;
        sol(2).Y = D.Y;
        sol(2).F = D.F;
    end

end
%--------------------------------------------------------------------------

% errors
errorU = abs(U-sol(1).U);
errorY = abs(Y-sol(1).Y);
errorF = abs(F-sol(1).F);

% outputs
O(1).value = max(errorY(:,1));
O(1).label = 'Ymax';
O(2).value = max(errorU(:,1));
O(2).label = 'Umax';
O(3).value = max(errorF);
O(3).label = 'F';
O(4).value = max(opts.QPcreatetime);
O(4).label = 'QPcreatetime';
O(5).value = max(opts.QPsolvetime);
O(5).label = 'QPsolvetime';