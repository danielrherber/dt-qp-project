%--------------------------------------------------------------------------
% BrysonHo59_solution.m
% Create solution for BrysonHo59 example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [U,Y,F] = BrysonHo59_solution(tf,a,h,T)

% fsolve options
options = optimoptions('fsolve','OptimalityTolerance',eps,...
    'FunctionTolerance',eps,'StepTolerance',eps);

% implicit function for B0
FUN = @(B0) -(4*h)/(a*tf^2) + 1/sin(B0) - (log((sec(B0)+tan(B0))/(sec(B0)-tan(B0))))/(2*tan(B0)^2);

% determine B0
[B0,FVAL] = fsolve(@(B0) FUN(B0),1,options); %#ok<ASGLU>

% compute intermediate value c
c = 2*tan(B0)/tf;

% compute controls
U = atan(tan(B0)*(1-2*T/tf));

% compute states
Y(:,1) = a/c*log( (tan(B0) + sec(B0))./(tan(U) + sec(U)) );
Y(:,2) = a/c*(sec(B0) - sec(U));
Y(:,3) = a/c^2*(sec(B0) - sec(U) - ...
    tan(U).*log((tan(B0) + sec(B0))./(tan(U) + sec(U))));
Y(:,4) = a/(2*c^2)*((tan(B0) - tan(U)).*sec(B0) - ...
    (sec(B0)-sec(U)).*tan(U) - log((tan(B0) + sec(B0))./(tan(U) + sec(U))));

% compute objective
F = a*tf*log((sec(B0)+tan(B0))/(sec(B0)-tan(B0)))./(2*tan(B0));

end