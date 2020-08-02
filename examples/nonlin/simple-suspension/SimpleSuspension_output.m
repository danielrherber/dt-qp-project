%--------------------------------------------------------------------------
% SimpleSuspension_output.m
% Output function for SimpleSuspension example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [O,sol] = SimpleSuspension_output(T,U,Y,P,F,in,opts)

% solution
sol = []; % no exact solution

% outputs
O(1).value = max(in.QPcreatetime);
O(1).label = 'QPcreatetime';
O(2).value = max(in.QPsolvetime);
O(2).label = 'QPsolvetime';

% determine which method
if strcmpi(opts.solver.function,'quadprog')
    method = "nested";
    P = in.p.xpopt;
else
    method = "simultaneous";
end

% display to command window
disp(strcat("--- ",method," method ---"))
disp(strcat("bs: ",string(P(1))))
disp(strcat("ks: ",string(P(2))))
disp(strcat(" F: ",string(F)))

end